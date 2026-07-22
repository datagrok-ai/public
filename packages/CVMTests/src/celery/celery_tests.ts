import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

import {category, test, expect, expectTable, expectExceptionAsync} from '@datagrok-libraries/test/src/test';

import {randomString, isEqualBytes, escapingTestStrings} from '../utils/test-utils';

category('Celery: datagrok-celery-task', () => {
  test('Scalars round-trip', async () => {
    const scalars: {[func: string]: any} = {cvmInt: 42, cvmDouble: 0.3, cvmBool: true, cvmString: 'Datagrok'};
    for (const [func, value] of Object.entries(scalars))
      expect(await grok.functions.call(`CVMTests:${func}`, {x: value}), value);
    const big = '9007199254740993'; // 2^53 + 1 — survives only as a string
    expect(String(await grok.functions.call('CVMTests:cvmBigInt', {x: big})), big);
  }, {timeout: 120000 /* long timeout: first call starts the worker container */});

  test('String escaping', async () => {
    for (const s of escapingTestStrings)
      expect(await grok.functions.call('CVMTests:cvmString', {x: s}), s);
  }, {timeout: 90000});

  test('Long string', async () => {
    const long = randomString(500000, '0123456789abcdefghijklmnopqrstuvwxyz');
    expect(await grok.functions.call('CVMTests:cvmString', {x: long}), long);
  }, {timeout: 90000});

  test('DataFrame round-trip', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'i', [1, 2, 3]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'd', [1.5, 2.5, 3.5]),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 's', ['a', 'x,y', 'Ünï']),
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'b', [true, false, true]),
    ]);
    expectTable(await grok.functions.call('CVMTests:cvmDataframe', {df}), df);
  }, {timeout: 90000});

  test('DataFrame null fidelity', async () => {
    const result: DG.DataFrame = await grok.functions.call('CVMTests:cvmDataframeNulls', {});
    expect(result.getCol('f').isNone(1), true);
    expect(result.getCol('f').get(0), 1.0);
    expect(result.getCol('f').get(2), 3.0);
  }, {timeout: 90000});

  test('Int column boundary typing', async () => {
    const inBound: DG.DataFrame = await grok.functions.call('CVMTests:cvmIntInBound', {});
    expect(inBound.getCol('col1').type === DG.COLUMN_TYPE.INT, true);
    const outBound: DG.DataFrame = await grok.functions.call('CVMTests:cvmIntOutBound', {});
    expect(outBound.getCol('col1').type === DG.COLUMN_TYPE.BIG_INT, true);
  }, {timeout: 90000});

  test('Empty dataframe', async () => {
    const result: DG.DataFrame = await grok.functions.call('CVMTests:cvmEmptyDataframe', {});
    expect(result.rowCount, 0);
  }, {timeout: 90000});

  test('Blob streaming (multi-chunk)', async () => {
    expect(await roundTripBytes('CVMTests:cvmBlob', 'blob', 5 * 1024 * 1024), true);
  }, {timeout: 120000});

  test('File round-trip', async () => {
    expect(await roundTripBytes('CVMTests:cvmFile', 'file', 100), true);
  }, {timeout: 90000});

  test('Datetime tz round-trip', async () => {
    const now = dayjs();
    const result = await grok.functions.call('CVMTests:cvmDate', {dt: now});
    expect(result.valueOf(), now.add(1, 'day').valueOf());
  }, {timeout: 90000});

  test('USER_API_KEY injection', async () => {
    const key: string = await grok.functions.call('CVMTests:cvmApiKey', {});
    expect(typeof key === 'string' && key.length > 0, true);
  }, {timeout: 90000});

  test('Cancellation', async () => {
    const call = DG.Func.find({package: 'CVMTests', name: 'cvmCancel'})[0].prepare({seconds: 60});
    const promise = call.call().catch(() => {});
    await DG.delay(5000);
    await call.cancel();
    await promise;
    expect(call.status, 'Canceled');
  }, {timeout: 120000});

  test('Single-output guard', async () => {
    await expectExceptionAsync(
      async () => {
        await grok.functions.call('CVMTests:cvmTwoOutputs', {});
      },
      (e) => String(e).includes('Only one return parameter is allowed'));
  }, {timeout: 90000});

  test('Multi-input streaming', async () => {
    const size = 1234;
    const blob = makeBytes(size);
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'a', [1, 2, 3, 4]),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'b', ['p', 'q', 'r', 's']),
    ]);
    const result: string = await grok.functions.call('CVMTests:cvmMultiInput',
      {blob: DG.FileInfo.fromBytes('mi.bin', blob), df});
    expect(result, `${size}:4:2`);
  }, {timeout: 90000});
}, {node: true});

function makeBytes(size: number): Uint8Array {
  const data = new Uint8Array(size);
  for (let i = 0; i < size; i++)
    data[i] = i & 0xff;
  return data;
}

async function roundTripBytes(funcName: string, paramName: string, size: number): Promise<boolean> {
  const data = makeBytes(size);
  const result = await grok.functions.call(funcName,
    {[paramName]: DG.FileInfo.fromBytes(`bytes_${size}.bin`, data)});
  return isEqualBytes(data, (result as DG.FileInfo).data);
}
