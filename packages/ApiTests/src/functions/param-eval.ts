import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

grok.functions.register({
  signature: 'List<String> paramEvalCities(String region)',
  run: (region: string) => [`${region}-1`, `${region}-2`],
});

grok.functions.register({
  signature: 'List<String> paramEvalFruits()',
  run: () => ['apple', 'banana', 'cherry'],
});

grok.functions.register({
  signature: 'List<String> paramEvalSuggest(String s)',
  run: (s: string) => [`${s}-hit`],
});

const host = grok.functions.register({
  signature: 'string paramEvalHost(string region, string city, string fruit, string name, int x, int broken)',
  run: () => '',
});

for (const input of host.inputs) {
  if (input.name === 'city')
    input.options['choices'] = 'paramEvalCities';
  else if (input.name === 'fruit')
    input.options['choices'] = 'paramEvalFruits()';
  else if (input.name === 'name')
    input.options['suggestions'] = 'paramEvalSuggest';
  else if (input.name === 'x')
    input.options['default'] = '2 + 2';
  else if (input.name === 'broken')
    input.options['default'] = 'paramEvalNoSuchFunc(1)';
}

grok.functions.register({
  signature: 'dataframe paramEvalCars()',
  run: () => DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'model', ['Mazda', 'Volvo']),
    DG.Column.fromList('int', 'mpg', [21, 30]),
    DG.Column.fromList('int', 'cyl', [6, 4]),
  ]),
});

const carHost = grok.functions.register({
  signature: 'string paramEvalCarHost(string model, string speed)',
  run: () => '',
});

for (const input of carHost.inputs) {
  if (input.name === 'model') {
    input.options['choices'] = 'paramEvalCars()';
    input.options['propagateChoice'] = 'all';
  } else if (input.name === 'speed')
    input.options['choices'] = '["slow", "fast"]';
}

async function expectRejection(action: () => Promise<any>): Promise<void> {
  try {
    await action();
  } catch (e) {
    return;
  }
  throw 'Expected a rejection';
}

category('Functions: ParamEval', () => {
  test('choices from client-registered provider', async () => {
    const call = host.prepare({region: 'EU'});
    const r = await call.evalParamChoices('city');
    expectArray(r.items, ['EU-1', 'EU-2']);
    expect(r.values['EU-1'], 'EU-1');
    expect(r.lookup == null, true);
    expectArray(r.dependsOn, ['region']);
  });

  test('choices follow the dependent param', async () => {
    const call = host.prepare({region: 'EU'});
    call.setParamValue('region', 'US');
    const r = await call.evalParamChoices('city');
    expectArray(r.items, ['US-1', 'US-2']);
  });

  test('choices without dependencies', async () => {
    const call = host.prepare();
    const r = await call.evalParamChoices('fruit');
    expectArray(r.items, ['apple', 'banana', 'cherry']);
    expect(r.values['banana'], 'banana');
    expectArray(r.dependsOn, []);
  });

  test('propagate lookup from a dataframe provider', async () => {
    const call = carHost.prepare();
    const r = await call.evalParamChoices('model');
    expectArray(r.items, ['Mazda', 'Volvo']);
    expect(typeof r.items[0], 'string');
    expect(r.values['Mazda'], 'Mazda');
    expect(r.lookup!['Mazda']['mpg'], 21);
    expect(r.lookup!['Mazda']['cyl'], 6);
    expect(r.lookup!['Volvo']['mpg'], 30);
    expect('model' in r.lookup!['Mazda'], false);
    expectArray(r.dependsOn, []);
  });

  test('static list-literal choices', async () => {
    const call = carHost.prepare();
    const r = await call.evalParamChoices('speed');
    expectArray(r.items, ['slow', 'fast']);
    expectArray(r.dependsOn, []);
  });

  test('suggestions receive the typed text', async () => {
    const call = host.prepare({name: 'stale'});
    const r = await call.evalParamSuggestions('name', 'typed');
    expectArray(r.items, ['typed-hit']);
    expect(r.tooltips['typed-hit'], 'typed-hit');
  });

  test('default from literal command', async () => {
    const call = host.prepare();
    expect(await call.evalParamDefault('x'), 4);
  });

  test('default from broken command rejects', async () => {
    const call = host.prepare();
    await expectRejection(() => call.evalParamDefault('broken'));
  });
}, {owner: 'askalkin@datagrok.ai'});

category('Functions: ScriptSync', () => {
  test('one-argument call', async () => {
    expect(grok.functions.scriptSync('2 + 3'), 5);
  });

  test('variables map', async () => {
    expect(grok.functions.scriptSync('a + b', {a: 2, b: 3}), 5);
  });

  test('variables map evaluates over a fresh context', async () => {
    grok.shell.setVar('paramEvalShellVar', 41 as any);
    expect(grok.functions.scriptSync('paramEvalShellVar + 1'), 42);
    let thrown = false;
    try {
      grok.functions.scriptSync('paramEvalShellVar + 1', {a: 1});
    } catch (e) {
      thrown = true;
    }
    expect(thrown, true);
  });
}, {owner: 'askalkin@datagrok.ai'});

category('Functions: Property options', () => {
  const propHost = grok.functions.register({
    signature: 'string paramEvalPropHost(string flavor, double dose)',
    run: () => '',
  });

  test('write from JS reads back through a fresh wrapper', async () => {
    const dose = propHost.inputs.find((p) => p.name === 'dose')!;
    dose.options['units'] = 'mg';
    dose.options['suggestions'] = 'paramEvalSuggest';
    const reread = propHost.inputs.find((p) => p.name === 'dose')!;
    expect(reread.options['units'], 'mg');
    expect(reread.options['suggestions'], 'paramEvalSuggest');
  });

  test('typed choices setter reads back', async () => {
    const flavor = propHost.inputs.find((p) => p.name === 'flavor')!;
    flavor.choices = ['vanilla', 'mint'];
    expectArray(propHost.inputs.find((p) => p.name === 'flavor')!.choices, ['vanilla', 'mint']);
  });

  test('viewer property options include tags', async () => {
    const viewer = DG.Viewer.scatterPlot(grok.data.demo.demog(20));
    const color = viewer.getProperties().find((p) => p.name === 'colorColumnName')!;
    try {
      expect(color.options['.is-legend-property'] != null, true);
      // viewer Look property descriptors are class-level and their wrapper (with its merged
      // options map) is cached, so a write survives re-reads session-wide — treat as readonly
      color.options['paramEvalProbe'] = 'x';
      const reread = viewer.getProperties().find((p) => p.name === 'colorColumnName')!;
      expect(reread.options['paramEvalProbe'], 'x');
    } finally {
      delete color.options['paramEvalProbe'];
      viewer.detach();
    }
  });
}, {owner: 'askalkin@datagrok.ai'});

category('Functions: TableParams', () => {
  const tableHost = grok.functions.register({
    signature: 'string paramEvalTableHost(dataframe df, column age, column_list cols)',
    run: () => '',
  });

  const regHost = grok.functions.register({
    signature: 'string paramEvalTableReg(dataframe data, column onCol)',
    run: () => '',
  });

  const metaScript = () => DG.Script.create([
    '//name: paramEvalTableMeta',
    '//language: javascript',
    '//input: dataframe df',
    '//input: column age {type: numerical}',
    '//input: column_list cols {type: categorical; table: df}',
    'return;',
  ].join('\n'));

  test('dataframe write keeps identity, same-reference write fires onChanged once', async () => {
    const df = DG.DataFrame.fromCsv('age,sex,height\n25,M,180\n31,F,165');
    df.name = 'paramEvalTableDemog';
    grok.shell.addTable(df);
    try {
      const call = tableHost.prepare();
      let fired = 0;
      const sub = call.inputParams.df.onChanged.subscribe(() => fired++);
      try {
        call.setParamValue('df', df);
        call.setParamValue('df', df);
      } finally {
        sub.unsubscribe();
      }
      expect(call.inputs.df.dart === df.dart, true);
      expect(fired, 1);
    } finally {
      grok.shell.closeTable(df);
    }
  });

  test('string into a column param reads back as the resolver FuncCall', async () => {
    const df = DG.DataFrame.fromCsv('age,sex,height\n25,M,180\n31,F,165');
    df.name = 'paramEvalTableResolve';
    grok.shell.addTable(df);
    try {
      const call = tableHost.prepare();
      call.setParamValue('df', df);
      call.setParamValue('age', 'height');
      expect(call.inputs.age instanceof DG.FuncCall, true);
      expect(call.inputs.age instanceof DG.Column, false);
    } finally {
      grok.shell.closeTable(df);
    }
  });

  test('column array into a column_list reads back as a ColumnList', async () => {
    const df = DG.DataFrame.fromCsv('age,sex,height\n25,M,180\n31,F,165');
    df.name = 'paramEvalTableCols';
    grok.shell.addTable(df);
    try {
      const call = tableHost.prepare();
      call.setParamValue('df', df);
      call.setParamValue('cols', [df.columns.byName('age'), df.columns.byName('height')]);
      const v = call.inputs.cols;
      expect(Array.isArray(v), false);
      expect(v instanceof DG.ColumnList, true);
      expectArray(v.names(), ['age', 'height']);
      expect(v.length, 2);
      const list = v.toList();
      expect(list.length, 2);
      expect(list[0] instanceof DG.Column, true);
      expect(list[1] instanceof DG.Column, true);
    } finally {
      grok.shell.closeTable(df);
    }
  });

  test('parentTableParamName answers both annotation link shapes', async () => {
    const get = (window as any).grok_Property_Get;
    expect(typeof get, 'function');
    const s = metaScript();
    const age = s.inputs.find((p) => p.name === 'age')!;
    const cols = s.inputs.find((p) => p.name === 'cols')!;
    expect(get(age.dart, 'parentTableParamName'), 'df');
    expect(age.options['table'] == null, true);
    expect(get(cols.dart, 'parentTableParamName'), 'df');
    expect(cols.options['table'], 'df');
  });

  test('grok_Property_Set round-trips on a registered func param', async () => {
    const get = (window as any).grok_Property_Get;
    const set = (window as any).grok_Property_Set;
    const onCol = regHost.inputs.find((p) => p.name === 'onCol')!;
    expect(get(onCol.dart, 'parentTableParamName') == null, true);
    try {
      set(onCol.dart, 'parentTableParamName', 'data');
      expect(get(onCol.dart, 'parentTableParamName'), 'data');
    } finally {
      // regHost is module-scope: a leftover link would fail the null expectation above on a
      // second same-session run
      set(onCol.dart, 'parentTableParamName', null);
    }
  });

  test('columnTypeFilter derives from the type option', async () => {
    const s = metaScript();
    expect(s.inputs.find((p) => p.name === 'age')!.columnTypeFilter, 'numerical');
    expect(s.inputs.find((p) => p.name === 'cols')!.columnTypeFilter, 'categorical');
  });
}, {owner: 'askalkin@datagrok.ai'});
