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
