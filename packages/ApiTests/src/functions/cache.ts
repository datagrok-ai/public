import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {DataFrame, DataQuery, Qnum} from 'datagrok-api/dg';
import dayjs from 'dayjs';
import * as DG from 'datagrok-api/dg';
import {check} from "./utils";

const cars = grok.data.demo.demog(10000);
const cacheName = 'function_results_cache';

const dup = grok.functions.register({
  signature: 'int test_dup(int x)',
  run: (x: number) => { return 2 * x; },
  options: { cache: 'true'}
});

const demog = grok.functions.register({
  signature: 'dataframe test_demog(int x)',
  run: (x: number) => { return cars; },
  options: { cache: 'true'}
});

const exc = grok.functions.register({
  signature: 'int test_exc(int x)',
  run: (x: number) => { throw 'My exception'; },
  options: { cache: 'true'}
});

const allInputTypes = grok.functions.register({
  signature: 'int test_allInputTypes(int i, double d, string s)',
  run: (i: number, d: number, s: string) => { `${i}-${d}-${s}`; },
  options: { cache: 'true'}
});

/** Creates a function that returns the parameter */
function echo(type: string) {
  return grok.functions.register({
    signature: `${type} test_return_${type}(${type} x)`,
    run: (x: any) => { return x; },
    options: { cache: 'true'}
  });
}

category('Functions: Client-side cache', () => {

  // before(async () => {
  //   await grok.functions.clientCache.clear();
  //   await grok.functions.clientCache.start();
  // });

  test('String', async () => {
    let f = grok.functions.register({
      signature: 'string getRandomStringDate()',
      run: () => { return new Date().toISOString(); },
      options: {'cache': 'true', 'cache.invalidateOn': '0 * * * *'},
    });

    await expectSameResults(f);
  });

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
  })

  test('10K calls no cache', async () => {
    grok.shell.settings.clientSideCache = false;
    for (let i = 0; i < 10000; i++)
      await dup.apply({x : i});
    grok.shell.settings.clientSideCache = true;
  });

  test('10K calls with cache', async () => {
    grok.shell.settings.clientSideCache = true;
    for (let i = 0; i < 10000; i++)
      await dup.apply({x : i});
  });

  test('1K demog10K', async () => {
    grok.shell.settings.clientSideCache = true;
    for (let i = 0; i < 1000; i++)
      await demog.apply({x : i});
  });

  test('Expiration', async () => {
    await grok.functions.clientCache.clear();
    await dup.apply({x: 7});
    await delay(10);
    await grok.functions.clientCache.cleanup();
    await dup.apply({x: 7});
  });
});


async function expectSameResults(f: DG.Func, params?: object): Promise<any> {
  const first = await f.apply(params);
  await delay(100);
  const second = await f.apply(params);
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
  }
  catch (e) {
    caught = true;
    checked = !check || check(e);
  }
  finally {
    if (!caught)
      throw 'An exception is expected but not thrown';
    if (!checked)
      throw 'An expected exception is thrown, but it does not satisfy the condition';
  }
}

