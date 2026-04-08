import './date-functions';
import './math-functions';
import './text-functions';
import './logical-functions';
import './conversion-functions';
import './stats-functions';
import {category, test, expect, expectExceptionAsync} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

category('Functions: Multiple outputs', () => {
  test('sync: register and call', async () => {
    const f = grok.functions.register({
      signature: '({int sum, int product}) testMultiOutSync(int a, int b)',
      run: (a: number, b: number) => ({sum: a + b, product: a * b}),
    });
    const call = f.prepare({a: 3, b: 4});
    await call.call();
    expect(call.outputs.get('sum'), 7, 'sum');
    expect(call.outputs.get('product'), 12, 'product');
  });

  test('async: register and call', async () => {
    const f = grok.functions.register({
      signature: '({string greeting, int length}) testMultiOutAsync(string name)',
      run: async (name: string) => ({greeting: `hello ${name}`, length: name.length}),
      isAsync: true,
    });
    const call = f.prepare({name: 'world'});
    await call.call();
    expect(call.outputs.get('greeting'), 'hello world', 'greeting');
    expect(call.outputs.get('length'), 5, 'length');
  });

  test('output params metadata', async () => {
    const f = grok.functions.register({
      signature: '({double x, double y}) testMultiOutMeta(double angle)',
      run: (angle: number) => ({x: Math.cos(angle), y: Math.sin(angle)}),
    });
    expect(f.outputs.length, 2, 'output count');
    expect(f.outputs[0].name, 'x', 'first output name');
    expect(f.outputs[1].name, 'y', 'second output name');
  });

  test('with semantic type', async () => {
    const f = grok.functions.register({
      signature: '({string/Molecule mol, double score}) testMultiOutSemType(string smiles)',
      run: (smiles: string) => ({mol: smiles, score: 0.95}),
    });
    const call = f.prepare({smiles: 'c1ccccc1'});
    await call.call();
    expect(call.outputs.get('mol'), 'c1ccccc1', 'mol');
    expect(call.outputs.get('score'), 0.95, 'score');
  });
});

category('Functions: General', () => {
  test('eval', async () => {
    const dfList: DG.DataFrame[] = await grok.functions
      .eval('OpenServerFile("System:AppData/ApiTests/datasets/demog.csv")');
    expect(dfList[0].columns instanceof DG.ColumnList, true);
  }, {stressTest: true});

  test('call', async () => {
    const dfList: DG.DataFrame[] = await grok.functions
      .call('OpenServerFile', {'fullPath': 'System:AppData/ApiTests/datasets/demog.csv'});
    expect(dfList[0].columns instanceof DG.ColumnList, true);
  }, {stressTest: true});
  
  test('def param', async () => {
    await grok.functions.call('AddNewColumn', {table: grok.data.demo.demog(), expression: 'test', name: 'test'});
  });

  // GROK-13478
  test('call params', async () => {
    expect(await grok.functions.call('sin', {y: 0.5}), null, '{y: 0.5}');
    expect(await grok.functions.call('sin', {x: 0.5}), 0.479425538604203, '{x: 0.5}');
    expect(await grok.functions.call('sin', {X: 0.5}), 0.479425538604203, '{X: 0.5}');
    await expectExceptionAsync(() => grok.functions.call('qqqqqq'));
  });

  test('query params', async () => {
    // expect(await grok.data.query('ApiTests:dummyPackageQuery', {y: 0.5}), null, '{y: 0.5}');
    expect((await grok.data.query('ApiTests:dummyPackageQuery', {x: 0.5})).get('res', 0), 0.5, '{x: 0.5}');
    expect((await grok.data.query('ApiTests:dummyPackageQuery', {X: 0.5})).get('res', 0), 0.5, '{X: 0.5}');
    await expectExceptionAsync(() => grok.data.query('ApiTests:qqqqqq'));
  });

  test(`script's package load`, async () => {
    const sin: DG.Func = await grok.functions.eval('ApiTests:dummyPackageScript');
    expect(sin.package.nqName, 'ApiTests');
  }, {skipReason: 'GROK-15178'});

  test(`package func's package load`, async () => {
    const sin: DG.Func = await grok.functions.eval('ApiTests:dummyDataFrameFunction');
    expect(sin.package.nqName, 'ApiTests');
  }, {skipReason: 'GROK-15178'});

  test(`core package load`, async () => {
    const sin: DG.Func = await grok.functions.eval('Sin');
    expect(sin.package.nqName, 'core');
  }, {skipReason: 'GROK-15178'});
});
