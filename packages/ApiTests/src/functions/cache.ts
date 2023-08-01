/* eslint-disable */
// @ts-nocheck
import {after, before, category, delay, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {DataFrame, Func, Qnum} from 'datagrok-api/dg';
import * as DG from 'datagrok-api/dg';
import {check} from './utils';
import {expectTable} from '../package';

const cars = grok.data.demo.demog(10000);
const cacheName = 'function_results_cache';

const dup = grok.functions.register({
  signature: 'int test_dup(int x)',
  run: (x: number) => {return 2 * x;},
  options: {cache: 'true'},
});

const demog = grok.functions.register({
  signature: 'dataframe test_demog(int x)',
  run: (x: number) => {return cars;},
  options: {cache: 'true'},
});

const exc = grok.functions.register({
  signature: 'int test_exc(int x)',
  run: (x: number) => {throw 'My exception';},
  options: {cache: 'true'},
});

const allInputTypes = grok.functions.register({
  signature: 'int test_allInputTypes(int i, double d, string s)',
  run: (i: number, d: number, s: string) => {`${i}-${d}-${s}`;},
  options: {cache: 'true'},
});

/** Creates a function that returns the parameter */
function echo(type: string) {
  return grok.functions.register({
    signature: `${type} test_return_${type}(${type} x)`,
    run: (x: any) => {return x;},
    options: {cache: 'true'},
  });
}

category('Functions: Client-side cache', () => {

  // before(async () => {
  //   await grok.functions.clientCache.clear();
  //   await grok.functions.clientCache.start();
  // });

  test('All scalars', async () => {
    await expectSameResults(echo(DG.TYPE.STRING), {x: 'foo'});
    await expectSameResults(echo(DG.TYPE.INT), {x: 1234});
    await expectSameResults(echo(DG.TYPE.FLOAT), {x: 1234.56});
    await expectSameResults(echo(DG.TYPE.QNUM), {x: Qnum.less(5)});
  });

  test('Clear, GetRecordCount', async () => {
    await grok.functions.clientCache.clear();
    await dup.apply({x: 7});
    expect(await grok.functions.clientCache.getRecordCount(), 1);
    await grok.functions.clientCache.clear();
    expect(await grok.functions.clientCache.getRecordCount(), 0);
  });

  test('Exceptions', async () => {
    await expectExceptionAsync(async () => await exc.apply({x: 7}));
    await delay(10);
    await expectExceptionAsync(async () => await exc.apply({x: 7}));
  });

  test('10K calls no cache', async () => {
    grok.shell.settings.clientSideCache = false;
    for (let i = 0; i < 10000; i++)
      await dup.apply({x: i});
    grok.shell.settings.clientSideCache = true;
  });

  test('10K calls with cache', async () => {
    grok.shell.settings.clientSideCache = true;
    for (let i = 0; i < 10000; i++)
      await dup.apply({x: i});
  });

  test('1K demog10K', async () => {
    grok.shell.settings.clientSideCache = true;
    for (let i = 0; i < 1000; i++)
      await demog.apply({x: i});
  });

  test('Expiration', async () => {
    await grok.functions.clientCache.clear();
    await dup.apply({x: 7});
    await delay(10);
    await grok.functions.clientCache.cleanup();
    await dup.apply({x: 7});
  });

  test('Client function cache output', async () => {
    await testOutputCacheFunc('dataframe getNowDf()', () => DataFrame.fromCsv(`id,date
id1,${Date.now()}`));
    await testOutputCacheFunc('int getEpochDate()', () => Date.now());
    await testOutputCacheFunc('double getRandomDouble()', () => Math.random());
    await testOutputCacheFunc('bool getRandomBool()', () => Math.random() < 0.5);
    await testOutputCacheFunc('string getRandomStringDate()', () => new Date().toISOString());
    // await testOutputCacheFunc('List<double> getRandomIntList()', () => {
    //   const result = [];
    //   for (let i = 0; i < 5; i++)
    //     result.push(Math.random());
    //   return result;
    // });
  });

  test('Client function cache primitive argument test', async () => {
    await testFunctionPrimitiveParam('string strParam(string i)', (i: string) => Date.now() + i,
      {i: 'hello'}, {i: 'world'});
    await testFunctionPrimitiveParam('int intParam(int i)', (i: number) => Date.now() * i,
      {i: 1}, {i: 2});
    await testFunctionPrimitiveParam('double doubleParam(double i)', (i: number) => Date.now() * i,
      {i: 3.14}, {i: 22.01});
    await testFunctionPrimitiveParam('string boolParam(bool i)', (i: boolean) => Date.now() +
            (Math.random() < 0.5 && i).toString(),
    {i: true}, {i: false});
  });

  test('Client function cache performance', async () => {
    const func = registerFunc('double wait()', async () => {
      await delay(300);
      return Math.random();
    }, true);
    const firstTime = await getFunctionExecutionTime(func);
    const secondTime = await getFunctionExecutionTime(func);
    expect(firstTime > secondTime * 2, true, `The first execution time ${firstTime} ms
        is no more than twice the second execution time ${secondTime} ms`);
  });

  test('Client function cache clear test', async () => {
    const func = registerFunc('int getEpochDateClear()', () => Date.now());
    const res = await expectSameResults(func);
    await grok.functions.clientCache.clear();
    expect(res !== await func.apply(), true);
  });

  test('Cached DataFrame id diff', async () => {
    const df1 = await demog.apply({x: 1});
    const df2 = await demog.apply({x: 1});
    expect(df1.id !== df2.id, true);
  });
});

async function expectSameResults(f: DG.Func, params?: object): Promise<any> {
  const first = await f.apply(params);
  await delay(100);
  const second = await f.apply(params);
  if (first.constructor.name === 'DataFrame')
    expectTable(first, second);
  else if (Array.isArray(first))
    expectArray(first, second);
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

function registerFunc(signature: string, run: (param: any) => any, isAsync: boolean = false): Func {
  return grok.functions.register({
    signature: signature,
    run: run,
    options: {'cache': 'true', 'cache.invalidateOn': '0 * * * *'},
    isAsync: isAsync,
  });
}

async function testOutputCacheFunc(signature: string, run: () => any): Promise<any> {
  const func = registerFunc(signature, run);
  await expectSameResults(func);
}

async function getFunctionExecutionTime(f: Func, params: object = {}): Promise<number> {
  const start = Date.now();
  await f.apply(params);
  return Date.now() - start;
}

async function testFunctionPrimitiveParam(signature: string, run: (param: any) => any, params: object = {},
  nextParam: object = {}): Promise<any> {
  const func = registerFunc(signature, run);
  const res = await expectSameResults(func, params);
  const resNext = await func.apply(nextParam);
  expect(res !== resNext, true);
}
