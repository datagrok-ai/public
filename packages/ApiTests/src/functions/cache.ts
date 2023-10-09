/* eslint-disable */
// @ts-nocheck
import {before, category, delay, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {DataFrame, Func, Qnum, toDart} from 'datagrok-api/dg';
import * as DG from 'datagrok-api/dg';
import {check} from './utils';
import {expectTable} from '../package';
import dayjs from "dayjs";

const cacheName = 'function_results_cache';
const demogHeavy = grok.data.demo.demog(5000000);// 90mb
const demogLite = grok.data.demo.demog();

const demog = registerFunc('dataframe test_demog(string type)', (type: string) => {
  if (type === 'h')
    return demogHeavy;
  else
    return demogLite;
}, false, '0 0 23 * *');

/** Creates a function that returns the parameter */
function echo(type: string) {
  return registerFunc(`${type} test_return_${type}(${type} x)`, (x: any) => {return x;});
}

const tiny = registerFunc('double tiny(int x)',
    (x: number) => Math.random() * x, false, '0 0 23 * *')

const exc = registerFunc('int test_exc(int x)', (x: number) => {throw 'My exception'});

category('Functions: Client-side cache', () => {

  before(async () => {
    await grok.functions.clientCache.start();
    await grok.functions.clientCache.clear();
  });

  test('Clear, GetRecordCount', async () => {
    await grok.functions.clientCache.clear();
    await tiny.apply({x: 7});
    expect(await grok.functions.clientCache.getRecordCount(), 1);
    await grok.functions.clientCache.clear();
    expect(await grok.functions.clientCache.getRecordCount(), 0);
  });

  test('Clear: Additional', async () => {
    const func = registerFunc('int getEpochDateClear()', () => Date.now(), false, '0 0 23 * *');
    const res = await expectSameResults(func);
    await grok.functions.clientCache.clear();
    const res1 = await func.apply();
    expect(res !== res1, true);
  });

  test('Clear: by func', async () => {
    await grok.functions.clientCache.clear();
    await tiny.apply({'x': 1});
    await tiny.apply({'x': 2});
    await demog.apply({'type': 'l'});
    expect(await grok.functions.clientCache.getRecordCount(), 3);
    await grok.functions.clientCache.clear(tiny.id);
    expect(await grok.functions.clientCache.getRecordCount(), 1);
  });

  test('All scalars, datetime', async () => {
    await grok.functions.clientCache.clear();

    await expectSameResults(echo(DG.TYPE.STRING), {x: 'foo'});
    await expectSameResults(echo(DG.TYPE.INT), {x: 1234});
    await expectSameResults(echo(DG.TYPE.FLOAT), {x: 1234.56});
    await expectSameResults(echo(DG.TYPE.QNUM), {x: Qnum.less(5)});
    // await expectSameResults(echo(DG.TYPE.BIG_INT), {x: 9007199254740991n});
    await expectSameResults(echo(DG.TYPE.DATE_TIME), {x: dayjs.utc()})
  });

  test('Exceptions: shouldn\'t be cached', async () => {
    await grok.functions.clientCache.clear();
    await expectExceptionAsync(async () => await exc.apply({x: 7}));
    const count = await grok.functions.clientCache.getRecordCount();
    expect(count, 0);
  });

  test('Expiration: fastCache', async () => {
    const func = registerFunc('double randDouble()', () => Math.random(), false);
    const res1 = await func.apply();
    await delay(60000);
    const res2 = await func.apply();
    expect(res1 != res2, true)
  }, {timeout: 80000});

  test('Cached DataFrame id diff', async () => {
    const df1 = await demog.apply({'type': 'l'});
    const df2 = await demog.apply({'type': 'l'});
    expect(df1.id !== df2.id, true);
  });

  test('100k calls no cache', async () => {
    return toDart(await runLoop(false, tiny, getTinyGenerator(true, 100000)));
  }, {timeout: 400000});

  test('100k calls cache', async () => {
    await tiny.apply({'x': 1});
    return toDart(await runLoop(true, tiny, getTinyGenerator(true, 100000)));
  }, {timeout: 400000});

  test('5 heavy cached', async () => {
    return toDart(await runLoop(true, demog, getHeavyGenerator()));
  }, {timeout: 10000000000000});

  test('Records limit, tiny', async () => {
    await grok.functions.clientCache.clear();
    await demog.apply({'type': 'h'});

    const emptyCache = await runLoop(true, demog, getHeavyGenerator());

    await runLoop(true, tiny, getTinyGenerator(false, 100000));// populate with 100k records
    const after100k = await runLoop(true, demog, getHeavyGenerator());

    await runLoop(true, tiny, getTinyGenerator(false, 150000));// populate with 100 + 150k records
    const after250k = await runLoop(true, demog, getHeavyGenerator());

    await runLoop(true, tiny, getTinyGenerator(false, 250000));// populate with 250 + 250k records
    const after500k = await runLoop(true, demog, getHeavyGenerator());

    await runLoop(true, tiny, getTinyGenerator(false, 500000));// populate with 500 + 500k records
    const after1kk = await runLoop(true, demog, getHeavyGenerator());

    await runLoop(true, tiny, getTinyGenerator(false, 1000000));// populate with 1kk + 1kk records
    const after2kk = await runLoop(true, demog, getHeavyGenerator());

    return toDart({"empty": emptyCache, "100k": after100k, "250k": after250k, "500k": after500k, "1kk": after1kk, "2kk": after2kk});
  }, {timeout: 10000000000000, skipReason: 'Just for test purposes'});

  test('Records random stress test', async () => {
    await grok.functions.clientCache.clear();
    await demog.apply({'type': 'h'});
    let tiny = registerFunc('double tiny(int x)',
        (x: number) => Math.random() * x, false, '*/5 * * * *'); // every 5 minutes

    let otherDemog = registerFunc('dataframe test_demog(string type)', (type: string) => {
        return demogLite;
    }, false); //every minute invalidate

    const emptyCache = await runLoop(true, demog, getHeavyGenerator());
    await runLoop(true, otherDemog, getHeavyGenerator(5));
    await runLoop(true, tiny, getTinyGenerator(false, 100000));// populate with 100k records
    const after100k = await runLoop(true, demog, getHeavyGenerator());

    tiny = registerFunc('double tiny(int x)',
        (x: number) => Math.random() * x, false, '* * * * *'); // every 1 minute
    await runLoop(true, tiny, getTinyGenerator(false, 150000));// populate with 100 + 150k records
    const after250k = await runLoop(true, demog, getHeavyGenerator());

    await runLoop(true, tiny, getTinyGenerator(false, 250000));// populate with 250 + 250k records
    const after500k = await runLoop(true, otherDemog, getHeavyGenerator());

    await runLoop(true, tiny, getTinyGenerator(false, 500000));// populate with 500 + 500k records
    const after1kk = await runLoop(true, demog, getHeavyGenerator());

    await runLoop(true, tiny, getTinyGenerator(false, 1000000));// populate with 1kk + 1kk records
    const after2kk = await runLoop(true, otherDemog, getHeavyGenerator());

    return toDart({"empty": emptyCache, "100k": after100k, "250k": after250k, "500k": after500k, "1kk": after1kk, "2kk": after2kk});
  }, {timeout: 10000000000000});
});

async function expectSameResults(f: DG.Func, params?: object): Promise<any> {
  const first = await f.apply(params);
  await delay(100);
  const second = await f.apply(params);
  if (first.constructor.name === 'DataFrame')
    expectTable(first, second);
  else if (Array.isArray(first))
    expectArray(first, second);
  else if (first instanceof dayjs)
    expect(first.format('YYYY-MM-DDTHH:mm:ss'), second.format('YYYY-MM-DDTHH:mm:ss'));
  else
    expect(first, second);
  return second;
}

/** Expects an asynchronous {@link action} to throw an exception. Use {@link check} to perform
 * deeper inspection of the exception if necessary. */
async function expectExceptionAsync(action: () => Promise<void>, check?: (exception: any) => boolean): Promise<void> {
  let caught: boolean = false;
  let checked: boolean = false;
  try {
    await action();
  } catch (e) {
    caught = true;
    checked = !check || check(e);
  } finally {
    if (!caught)
      throw 'An exception is expected but not thrown';
    if (!checked)
      throw 'An expected exception is thrown, but it does not satisfy the condition';
  }
}

function registerFunc(signature: string, run: (param: any) => any, isAsync: boolean = false,
                      invalidateOn: string = '* * * * *'): Func {
  return grok.functions.register({
    signature: signature,
    run: run,
    options: {'cache': 'true', 'cache.invalidateOn': invalidateOn},
    isAsync: isAsync,
  });
}

async function getFunctionExecutionTime(f: Func, params: object = {}): Promise<number> {
  const start = Date.now();
  await f.apply(params);
  return Date.now() - start;
}

async function runLoop(cache: boolean, func: Func, argumentsProvider: Generator<void, object, object>): object {
  try {
    if (!cache) {
      grok.functions.clientCache.stop();
    }
    const results = [];
    for (let args of argumentsProvider())
      results.push(await getFunctionExecutionTime(func, args))

    const sum = results.reduce((p, c) => p + c, 0);
    return {'Average time': sum / results.length,
      'Min time': Math.min(...results), 'Max time': Math.max(...results)};
  } catch (e) {
    grok.log.error(e);
  }
}

function getTinyGenerator(same: boolean = true, count: number = 1000000): Generator<void, object, object> {
  return function* () {
    let n = 0;
    while (n < count) {
      yield {'x': same ? 1 : n};
      n++;
    }
  }
}

function getHeavyGenerator(count: number = 5): Generator<void, object, object> {
  return function* () {
    let n = 0;
    while (n < count) {
      yield {'type': 'h'};
      n++;
    }
  }
}

