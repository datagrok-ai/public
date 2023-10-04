/* eslint-disable */
// @ts-nocheck
import {before, category, delay, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {DataFrame, Func, Qnum, toDart} from 'datagrok-api/dg';
import * as DG from 'datagrok-api/dg';
import {check} from './utils';
import {expectTable} from '../package';

const tiny = registerFunc('double tiny()', () => Math.random(), false, '0 0 23 * *')

category('Functions: Client-side cache', () => {

  before(async () => {
    await grok.functions.clientCache.start();
    await grok.functions.clientCache.clear();
  });

  // test('Clear: Additional', async () => {
  //   const func = registerFunc('int getEpochDateClear()', () => Date.now());
  //   const res = await expectSameResults(func);
  //   await grok.functions.clientCache.clear();
  //   expect(res !== await func.apply(), true);
  // });
  //
  // test('Exceptions: shouldn\'t be cached', async () => {
  //   await grok.functions.clientCache.clear();
  //   await expectExceptionAsync(async () => await exc.apply({x: 7}));
  //   const count = await grok.functions.clientCache.getRecordCount();
  //   expect(count, 0);
  // });
  //


  test('Stopping working', async () => {
    await grok.functions.clientCache.stop();
    expect(grok.functions.clientCache.isRunning, false);
    await grok.functions.clientCache.start();
  });

  // test('Simple function shouldn\'t be cached', async () => {
  //   await grok.functions.clientCache.clear();
  //   await tiny.apply();
  //   expect(await grok.functions.clientCache.getRecordCount(), 0);
  // });
  //
  // test('Expiration', async () => {
  //   const func = registerFunc('double tiny()', () => Math.random(), false);
  //   const res1 = await func.apply();
  //   await delay(300000);
  //   const res2 = await func.apply();
  //   expect(res1 != res2, true)
  // }, {timeout: 80000});
  //
  // test('Cached DataFrame id diff', async () => {
  //   const df1 = await demog.apply({x: 1});
  //   const df2 = await demog.apply({x: 1});
  //   expect(df1.id !== df2.id, true);
  // });
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

function registerFunc(signature: string, run: (param: any) => any, isAsync: boolean = false,
                      invalidateOn: string = '0 * * * *'): Func {
  return grok.functions.register({
    signature: signature,
    run: run,
    options: {'cache': 'true', 'cache.invalidateOn': invalidateOn},
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
