/* eslint-disable */
import {before, category, delay, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {check} from './utils';
import {expectTable} from '../package';
import dayjs from "dayjs";

const demogHeavy = grok.data.demo.demog(5000000); // 90mb
const demogMiddle = grok.data.demo.demog(1100000); // 20mb
const demogLite = grok.data.demo.demog();

const demog = registerFunc('dataframe test_demog(string type)', (type: string) => {
  if (type === 'h')
    return demogHeavy;
  else if (type === 'm')
    return demogMiddle;
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


category('Benchmarks: Client-side cache', () => {
  before(async () => {
    await grok.functions.clientCache.start();
    await grok.functions.clientCache.clear();
  });

  test('Tiny scalar calls no cache', async () => {
    const iterations = DG.Test.isInBenchmark ? 40000 : 100;
    return DG.toDart(await runLoop(false, tiny, getTinyGenerator(true, iterations)));
  }, {timeout: 400000, benchmark: true});

  test('Tiny scalar calls with cache', async () => {
    await tiny.apply({'x': 1});
    const iterations = DG.Test.isInBenchmark ? 40000 : 100;
    return DG.toDart(await runLoop(true, tiny, getTinyGenerator(true, iterations)));
  }, {timeout: 400000, benchmark: true});

  test('Cached dataframe', async () => {
    const type = DG.Test.isInBenchmark ? 'h' : 'l';
    return DG.toDart(await runLoop(true, demog, getHeavyGenerator(10, type)));
  }, {timeout: 180000, benchmark: true});
}, { owner: 'ppolovyi@datagrok.ai'});

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
    await expectSameResults(echo(DG.TYPE.QNUM), {x: DG.Qnum.less(5)});
    // await expectSameResults(echo(DG.TYPE.BIG_INT), {x: 9007199254740991n});
    await expectSameResults(echo(DG.TYPE.DATE_TIME), {x: dayjs.utc()})
  });

  test('Exceptions: shouldn\'t be cached', async () => {
    await grok.functions.clientCache.clear();
    await expectExceptionAsync(async () => await exc.apply({x: 7}));
    const count = await grok.functions.clientCache.getRecordCount();
    expect(count, 0);
  });

  test('Expiration', async () => {
    const func = registerFunc('dataframe getNowDf()', () => DG.DataFrame.fromCsv(`id,date
id1,${Date.now()}`));
    const res1 = await func.apply();

    let transaction: any;
    let db: any;
    const cacheName = 'datagrok_function_cache';
    try {
      const request = window.indexedDB.open(cacheName);
      request.onsuccess = function() {
        db = request.result;
        transaction = db.transaction('cache_entry', 'readwrite');
        const entryStore = transaction.objectStore('cache_entry');
        entryStore.openCursor().onsuccess = (event: any) => {
          const cursor = event.target.result;
          if (cursor) {
            let updateEntry = cursor.value;
            if (updateEntry.__meta_id === func.id) {
              updateEntry.__expires = new Date("1996-08-26T03:24:00").getTime();
              const request = cursor.update(updateEntry);
              request.onsuccess = () => {};
            }
            cursor.continue();
          }
        };
      };
      if (transaction != null)
        transaction.commit();
    } catch (e) {
      if (transaction != null)
        transaction.abort();
      throw e;
    } finally {
      if (db != null)
        db.close;
    }

    await delay(500);
    const res2 = await func.apply();
    expect(res2.columns.byName('date').get(0) !== res1.columns.byName('date').get(0), true);
  }, {timeout: 80000});

  test('Cached DataFrame id diff', async () => {
    const df1 = await demog.apply({'type': 'l'});
    const df2 = await demog.apply({'type': 'l'});
    expect(df1.id !== df2.id, true);
  });
}, { owner: 'ppolovyi@datagrok.ai'});

async function expectSameResults(f: DG.Func, params?: object): Promise<any> {
  const first = await f.apply(params);
  await delay(100);
  const second = await f.apply(params);
  if (first.constructor.name === 'DataFrame')
    expectTable(first, second);
  else if (Array.isArray(first))
    expectArray(first, second);
  else if (first instanceof dayjs)
    expect((first as dayjs.Dayjs).format('YYYY-MM-DDTHH:mm:ss'), second.format('YYYY-MM-DDTHH:mm:ss'));
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
                      invalidateOn: string = '* * * * *'): DG.Func {
  return grok.functions.register({
    signature: signature,
    run: run,
    options: {'cache': 'true', 'cache.invalidateOn': invalidateOn},
    isAsync: isAsync,
  });
}

async function getFunctionExecutionTime(f: DG.Func, params: object = {}): Promise<number> {
  const start = Date.now();
  await f.apply(params);
  return Date.now() - start;
}

async function runLoop(cache: boolean, func: DG.Func, argumentsProvider: Generator<void, object, object>): Promise<object | undefined>  {
  try {
    if (!cache) {
      grok.functions.clientCache.stop();
    }
    const results = [];
    // @ts-ignore
    for (let args of argumentsProvider())
      results.push(await getFunctionExecutionTime(func, args))

    const sum = results.reduce((p, c) => p + c, 0);
    return {'Iterations' : results.length, 'Average time': sum / results.length,
      'Min time': Math.min(...results), 'Max time': Math.max(...results)};
  } catch (e) {
    grok.log.error(e);
  }
}

function getTinyGenerator(same: boolean = true, count: number = 1000000): Generator<void, object, object> {
  // @ts-ignore
  return function* () {
    let n = 0;
    while (n < count) {
      yield {'x': same ? 1 : n};
      n++;
    }
  }
}

function getHeavyGenerator(count: number = 5, type:string = 'h'): Generator<void, object, object> {
  // @ts-ignore
  return function* () {
    let n = 0;
    while (n < count) {
      yield {'type': type};
      n++;
    }
  }
}

