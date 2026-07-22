import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, expectTable, expectExceptionAsync} from '@datagrok-libraries/test/src/test';

import {randomString, escapingTestStrings} from '../utils/test-utils';

category('Celery: node worker', () => {
  test('Scalars round-trip', async () => {
    const scalars: {[func: string]: any} = {jsCvmInt: 42, jsCvmDouble: 0.3, jsCvmBool: true, jsCvmString: 'Datagrok'};
    for (const [func, value] of Object.entries(scalars))
      expect(await grok.functions.call(`CVMTests:${func}`, {x: value}), value);
    const big = '9007199254740993'; // 2^53 + 1 — survives only as a string
    expect(String(await grok.functions.call('CVMTests:jsCvmBigInt', {x: big})), big);
  }, {timeout: 120000 /* long timeout: first call starts the worker container */});

  test('String escaping', async () => {
    for (const s of escapingTestStrings)
      expect(await grok.functions.call('CVMTests:jsCvmString', {x: s}), s);
  }, {timeout: 90000});

  test('Long string', async () => {
    const long = randomString(500000, '0123456789abcdefghijklmnopqrstuvwxyz');
    expect(await grok.functions.call('CVMTests:jsCvmString', {x: long}), long);
  }, {timeout: 90000});

  test('DataFrame round-trip', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'i', [1, 2, 3]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'd', [1.5, 2.5, 3.5]),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 's', ['a', 'x,y', 'Ünï']),
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'b', [true, false, true]),
    ]);
    expectTable(await grok.functions.call('CVMTests:jsCvmDataframe', {df}), df);
  }, {timeout: 90000});

  test('Empty dataframe', async () => {
    const result: DG.DataFrame = await grok.functions.call('CVMTests:jsCvmEmptyDataframe', {});
    expect(result.rowCount, 0);
  }, {timeout: 90000});

  test('Error propagation', async () => {
    await expectExceptionAsync(
      async () => {
        await grok.functions.call('CVMTests:jsCvmError', {});
      },
      (e) => String(e).includes('planned jsCvmError failure'));
  }, {timeout: 90000});

  test('Progress call is a no-op-safe pass-through', async () => {
    expect(await grok.functions.call('CVMTests:jsCvmProgress', {}), 'ok');
  }, {timeout: 90000});

  test('Dapi access with caller token', async () => {
    const login: string = await grok.functions.call('CVMTests:jsCvmCurrentUser', {});
    expect(login, (await grok.dapi.users.current()).login);
  }, {timeout: 90000});

  test('Cancellation', async () => {
    const call = DG.Func.find({package: 'CVMTests', name: 'jsCvmCancel'})[0].prepare({seconds: 60});
    const promise = call.call().catch(() => {});
    await DG.delay(5000);
    await call.cancel();
    await promise;
    expect(call.status, 'Canceled');
  }, {timeout: 120000});

  test('meta.server alias', async () => {
    expect(await grok.functions.call('CVMTests:jsCvmServer', {x: 'srv'}), 'srv');
  }, {timeout: 90000});

  test('Custom container', async () => {
    expect(await grok.functions.call('CVMTests:jsCvmCustomContainer', {}), 'custom');
  }, {timeout: 240000 /* separate container cold start */});
});
