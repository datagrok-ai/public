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
}, {owner: 'askalkin@datagrok.ai'});
