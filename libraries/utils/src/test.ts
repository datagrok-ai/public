import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable, Subscription} from 'rxjs';
import Timeout = NodeJS.Timeout;
import {DataFrame} from 'datagrok-api/dg';

const STANDART_TIMEOUT = 30000;
const BENCHMARK_TIMEOUT = 10800000;

export const tests: {
  [key: string]: {
    tests?: Test[], before?: () => Promise<void>, after?: () => Promise<void>,
    beforeStatus?: string, afterStatus?: string, clear?: boolean, timeout?: number
  }
} = {};

const autoTestsCatName = 'Auto Tests';
const demoCatName = 'Demo';
const wasRegistered: {[key: string]: boolean} = {};
export let currentCategory: string;

export namespace assure {
  export function notNull(value: any, name?: string) {
    if (value == null)
      throw new Error(`${name == null ? 'Value' : name} not defined`);
  }
}

export interface TestOptions {
  timeout?: number;
  unhandledExceptionTimeout?: number;
  skipReason?: string;
}

export interface CategoryOptions {
  clear?: boolean;
  timeout?: number;
}

export class TestContext {
  catchUnhandled = true;
  report = false;

  constructor(catchUnhandled?: boolean, report?: boolean) {
    if (catchUnhandled !== undefined) this.catchUnhandled = catchUnhandled;
    if (report !== undefined) this.report = report;
  };
}

export class Test {
  test: () => Promise<any>;
  name: string;
  category: string;
  options?: TestOptions;

  constructor(category: string, name: string, test: () => Promise<any>, options?: TestOptions) {
    this.category = category;
    this.name = name;
    options ??= {};
    options.timeout ??= STANDART_TIMEOUT;
    this.options = options;
    this.test = async (): Promise<any> => {
      return new Promise(async (resolve, reject) => {
        let result = '';
        try {
          result = await test();
        } catch (e: any) {
          reject(e);
        }
        resolve(result);
      });
    };
  }
}

export async function testEvent<T>(event: Observable<T>,
  handler: (args: T) => void, trigger: () => void, ms: number = 0): Promise<string> {
  let sub: Subscription;
  return new Promise((resolve, reject) => {
    sub = event.subscribe((args: any) => {
      try {
        handler(args);
      } catch (e) {
        reject(e);
      }
      sub.unsubscribe();
      resolve('OK');
    });
    setTimeout(() => {
      sub.unsubscribe();
      // eslint-disable-next-line prefer-promise-reject-errors
      reject('timeout');
    }, ms);
    trigger();
  });
}

export function test(name: string, test: () => Promise<any>, options?: TestOptions): void {
  if (tests[currentCategory] == undefined)
    tests[currentCategory] = {};
  if (tests[currentCategory].tests == undefined)
    tests[currentCategory].tests = [];
  tests[currentCategory].tests!.push(new Test(currentCategory, name, test, options));
}

/* Tests two objects for equality, throws an exception if they are not equal. */
export function expect(actual: any, expected: any = true, error?: string): void {
  if (error)
    error = `${error}, `;
  else error = '';
  if (actual !== expected)
    throw new Error(`${error}Expected "${expected}", got "${actual}"`);
}

export function expectFloat(actual: number, expected: number, tolerance = 0.001, error?: string): void {
  if ((actual === Number.POSITIVE_INFINITY && expected === Number.POSITIVE_INFINITY) ||
      (actual === Number.NEGATIVE_INFINITY && expected === Number.NEGATIVE_INFINITY) ||
      (actual === Number.NaN && expected === Number.NaN) || (isNaN(actual) && isNaN(expected)))
    return;
  const areEqual = Math.abs(actual - expected) < tolerance;
  expect(areEqual, true, `${error ?? ''} (tolerance = ${tolerance})`);
  if (!areEqual)
    throw new Error(`Expected ${expected}, got ${actual} (tolerance = ${tolerance})`);
}

export function expectTable(actual: DataFrame, expected: DataFrame, error?: string): void {
  const expectedRowCount = expected.rowCount;
  const actualRowCount = actual.rowCount;
  expect(actualRowCount, expectedRowCount, `${error ?? ''}, row count`);

  for (const column of expected.columns) {
    const actualColumn = actual.columns.byName(column.name);
    if (actualColumn == null)
      throw new Error(`Column ${column.name} not found`);
    if (actualColumn.type != column.type)
      throw new Error(`Column ${column.name} type expected ${column.type} got ${actualColumn.type}`);
    for (let i = 0; i < expectedRowCount; i++) {
      const value = column.get(i);
      const actualValue = actualColumn.get(i);
      if (column.type == DG.TYPE.FLOAT)
        expectFloat(actualValue, value, 0.0001, error);
      else if (column.type == DG.TYPE.DATE_TIME)
        expect(actualValue.isSame(value), true, error);
      else
        expect(actualValue, value, error);
    }
  }
}

export function expectObject(actual: { [key: string]: any }, expected: { [key: string]: any }) {
  for (const [expectedKey, expectedValue] of Object.entries(expected)) {
    if (!actual.hasOwnProperty(expectedKey))
      throw new Error(`Expected property "${expectedKey}" not found`);

    const actualValue = actual[expectedKey];
    if (actualValue instanceof Array && expectedValue instanceof Array)
      expectArray(actualValue, expectedValue);
    else if (actualValue instanceof Object && expectedValue instanceof Object)
      expectObject(actualValue, expectedValue);
    else if (Number.isFinite(actualValue) && Number.isFinite(expectedValue))
      expectFloat(actualValue, expectedValue);
    else if (actualValue != expectedValue)
      throw new Error(`Expected (${expectedValue}) for key '${expectedKey}', got (${actualValue})`);
  }
}

export function expectArray(actual: ArrayLike<any>, expected: ArrayLike<any>) {
  const actualLength = actual.length;
  const expectedLength = expected.length;

  if (actualLength != expectedLength) {
    throw new Error(`Arrays are of different length: actual array length is ${actualLength} ` +
      `and expected array length is ${expectedLength}`);
  }

  for (let i = 0; i < actualLength; i++) {
    if (actual[i] instanceof Array && expected[i] instanceof Array)
      expectArray(actual[i], expected[i]);
    else if (actual[i] instanceof Object && expected[i] instanceof Object)
      expectObject(actual[i], expected[i]);
    else if (actual[i] != expected[i])
      throw new Error(`Expected ${expected[i]} at position ${i}, got ${actual[i]}`);
  }
}

/* Defines a test suite. */
export function category(category: string, tests_: () => void, options?: CategoryOptions): void {
  currentCategory = category;
  tests_();
  if (tests[currentCategory]) {
    tests[currentCategory].clear = options?.clear ?? true;
    tests[currentCategory].timeout = options?.timeout;
  }
}

/* Defines a function to be executed before the tests in this category are executed. */
export function before(before: () => Promise<void>): void {
  if (tests[currentCategory] == undefined)
    tests[currentCategory] = {};
  tests[currentCategory].before = before;
}

/* Defines a function to be executed after the tests in this category are executed. */
export function after(after: () => Promise<void>): void {
  if (tests[currentCategory] == undefined)
    tests[currentCategory] = {};
  tests[currentCategory].after = after;
}

function addNamespace(s: string, f: DG.Func): string {
  return s.replace(new RegExp(f.name, 'gi'), f.nqName);
}

export async function initAutoTests(packageId: string, module?: any) {
  if (wasRegistered[packageId]) return;
  const moduleTests = module ? module.tests : tests;
  if (moduleTests[autoTestsCatName] !== undefined ||
      moduleTests[demoCatName] !== undefined ||
      Object.keys(moduleTests).find((c) => c.startsWith(autoTestsCatName))) {
    wasRegistered[packageId] = true;
    return;
  }
  const moduleAutoTests = [];
  const moduleDemo = [];
  const packFunctions = await grok.dapi.functions.filter(`package.id = "${packageId}"`).list();
  const reg = new RegExp(/skip:\s*([^,\s]+)|wait:\s*(\d+)|cat:\s*([^,\s]+)/g);
  for (const f of packFunctions) {
    const tests = f.options['test'];
    const demo = f.options['demoPath'];
    if ((tests && Array.isArray(tests) && tests.length)) {
      for (let i = 0; i < tests.length; i++) {
        const res = (tests[i] as string).matchAll(reg);
        const map: {skip?: string, wait?: number, cat?: string} = {};
        Array.from(res).forEach((arr) => {
          if (arr[0].startsWith('skip')) map['skip'] = arr[1];
          else if (arr[0].startsWith('wait')) map['wait'] = parseInt(arr[2]);
          else if (arr[0].startsWith('cat')) map['cat'] = arr[3];
        });
        const test = new Test(autoTestsCatName, tests.length === 1 ? f.name : `${f.name} ${i + 1}`, async () => {
          const res = await grok.functions.eval(addNamespace(tests[i], f));
          if (map.wait) await delay(map.wait);
          // eslint-disable-next-line no-throw-literal
          if (typeof res === 'boolean' && !res) throw `Failed: ${tests[i]}, expected true, got ${res}`;
        }, {skipReason: map.skip});
        if (map.cat) {
          const cat: string = autoTestsCatName + ': ' + map.cat;
          test.category = cat;
          if (moduleTests[cat] === undefined)
            moduleTests[cat] = {tests: [], clear: true};
          moduleTests[cat].tests.push(test);
        } else {
          moduleAutoTests.push(test);
        }
      }
    }
    if (demo) {
      const wait = f.options['demoWait'] ? parseInt(f.options['demoWait']) : undefined;
      const test = new Test(demoCatName, f.friendlyName, async () => {
        grok.shell.lastError = '';
        await f.apply();
        await delay(wait ? wait : 1000);
        if (grok.shell.lastError)
          throw new Error(grok.shell.lastError);
      });
      moduleDemo.push(test);
    }
  }
  wasRegistered[packageId] = true;
  if (moduleAutoTests.length)
    moduleTests[autoTestsCatName] = {tests: moduleAutoTests, clear: true};
  if (moduleDemo.length)
    moduleTests[demoCatName] = {tests: moduleDemo, clear: true};
}

export async function runTests(options?: {category?: string, test?: string, testContext?: TestContext}) {
  const package_ = grok.functions.getCurrentCall()?.func?.package;
  await initAutoTests(package_.id);
  const results: { category?: string, name?: string, success: boolean,
                   result: string, ms: number, skipped: boolean }[] = [];
  console.log(`Running tests`);
  options ??= {};
  options!.testContext ??= new TestContext();
  grok.shell.lastError = '';
  for (const [key, value] of Object.entries(tests)) {
    if (options?.category != undefined) {
      if (!key.toLowerCase().startsWith(options?.category.toLowerCase()))
        continue;
    }
    console.log(`Started ${key} category`);
    try {
      if (value.before)
        await value.before();
    } catch (x: any) {
      value.beforeStatus = x.toString();
    }
    const t = value.tests ?? [];
    const res = [];
    if (value.clear) {
      for (let i = 0; i < t.length; i++) {
        res.push(await execTest(t[i], options?.test, value.timeout, package_.name));
        grok.shell.closeAll();
        DG.Balloon.closeAll();
      }
    } else {
      for (let i = 0; i < t.length; i++)
        res.push(await execTest(t[i], options?.test, value.timeout, package_.name));
    }
    const data = (await Promise.all(res)).filter((d) => d.result != 'skipped');
    try {
      if (value.after)
        await value.after();
    } catch (x: any) {
      value.afterStatus = x.toString();
    }
    // Clear after category
    // grok.shell.closeAll();
    // DG.Balloon.closeAll();
    if (value.afterStatus)
      data.push({category: key, name: 'init', result: value.afterStatus, success: false, ms: 0, skipped: false});
    if (value.beforeStatus)
      data.push({category: key, name: 'init', result: value.beforeStatus, success: false, ms: 0, skipped: false});
    results.push(...data);
  }
  if (options.testContext.catchUnhandled) {
    await delay(1000);
    if (grok.shell.lastError.length > 0) {
      results.push({
        category: 'Unhandled exceptions',
        name: 'exceptions',
        result: grok.shell.lastError, success: false, ms: 0, skipped: false
      });
    }
  }
  if (options.testContext.report) {
    const logger = new DG.Logger();
    const successful = results.filter((r) => r.success).length;
    const skipped = results.filter((r) => r.skipped).length;
    const failed = results.filter((r) => !r.success);
    const description = 'Package @package tested: @successful successful, @skipped skipped, @failed failed tests';
    const params = {
      successful: successful,
      skipped: skipped,
      failed: failed.length,
      package: package_
    };
    for (const r of failed) Object.assign(params, {[`${r.category} | ${r.name}`]: r.result});
    logger.log(description, params, 'package-tested');
  }
  return results;
}

async function execTest(t: Test, predicate: string | undefined, categoryTimeout?: number, packageName?: string) {
  let r: { category?: string, name?: string, success: boolean, result: any, ms: number, skipped: boolean };
  const filter = predicate != undefined && (t.name.toLowerCase() !== predicate.toLowerCase());
  const skip = t.options?.skipReason || filter;
  const skipReason = filter ? 'skipped' : t.options?.skipReason;
  if (!skip)
    console.log(`Started ${t.category} ${t.name}`);
  const start = new Date();
  try {
    if (skip) {
      r = {success: true, result: skipReason!, ms: 0, skipped: true};
    } else {
      let timeout_ = t.options?.timeout === STANDART_TIMEOUT &&
        categoryTimeout ? categoryTimeout : t.options?.timeout!;
      timeout_ = DG.Test.isInBenchmark && timeout_ === STANDART_TIMEOUT ? BENCHMARK_TIMEOUT : timeout_;
      r = {success: true, result: await timeout(t.test, timeout_) ?? 'OK', ms: 0, skipped: false};
    }
  } catch (x: any) {
    r = {success: false, result: x.toString(), ms: 0, skipped: false};
  }
  const stop = new Date();
  // @ts-ignore
  r.ms = stop - start;
  if (!skip)
    console.log(`Finished ${t.category} ${t.name} for ${r.ms} ms`);
  r.category = t.category;
  r.name = t.name;
  if (!filter) {
    let params = {'success': r.success, 'result': r.result, 'ms': r.ms, 'skipped': r.skipped,
      'type': 'package', packageName, 'category': t.category, 'test': t.name};
    if (r.result.constructor == Object) {
      const res = Object.keys(r.result).reduce((acc, k) => ({...acc, ['result.' + k]: r.result[k]}), {});
      params = {...params, ...res};
    }
    grok.log.usage(`${packageName}: ${t.category}: ${t.name}`,
      params, `test-package ${packageName}: ${t.category}: ${t.name}`);
  }
  return r;
}

/* Waits [ms] milliseconds */
export async function delay(ms: number) {
  await new Promise((r) => setTimeout(r, ms));
}

export async function awaitCheck(checkHandler: () => boolean,
  error: string = 'Timeout exceeded', wait: number = 500, interval: number = 50): Promise<void> {
  return new Promise((resolve, reject) => {
    setTimeout(() => {
      clearInterval(intervalId);
      reject(new Error(error));
    }, wait);
    // @ts-ignore
    const intervalId: Timeout = setInterval(() => {
      if (checkHandler()) {
        clearInterval(intervalId);
        resolve();
      }
    }, interval);
  });
}

async function timeout(func: () => Promise<any>, testTimeout: number): Promise<any> {
  let timeout: Timeout | null = null;
  const timeoutPromise = new Promise<any>((_, reject) => {
    //@ts-ignore
    timeout = setTimeout(() => {
      // eslint-disable-next-line prefer-promise-reject-errors
      reject('EXECUTION TIMEOUT');
    }, testTimeout);
  });
  try {
    return await Promise.race([func(), timeoutPromise]);
  } finally {
    if (timeout)
      clearTimeout(timeout);
  }
}

export function isDialogPresent(dialogTitle: string): boolean {
  const dialogs = DG.Dialog.getOpenDialogs();
  for (let i = 0; i < dialogs.length; i++) {
    if (dialogs[i].title == dialogTitle)
      return true;
  }
  return false;
}

/** Expects an asynchronous {@link action} to throw an exception. Use {@link check} to perform
 * deeper inspection of the exception if necessary.
 * @param  {function(): Promise<void>} action
 * @param  {function(any): boolean} check
 * @return {Promise<void>}
 */
export async function expectExceptionAsync(action: () => Promise<void>,
  check?: (exception: any) => boolean): Promise<void> {
  let caught: boolean = false;
  let checked: boolean = false;
  try {
    await action();
  } catch (e) {
    caught = true;
    checked = !check || check(e);
  } finally {
    if (!caught)
      throw new Error('An exception is expected but not thrown');
    if (!checked)
      throw new Error('An expected exception is thrown, but it does not satisfy the condition');
  }
}

/**
 * Universal test for viewers. It search viewers in DOM by tags: canvas, svg, img, input, h1, a
 * @param  {string} v Viewer name
 * @param  {DG.DataFrame} df Dataframe to use. Should have at least 3 rows
 * @param  {boolean} options.detectSemanticTypes Specify whether to detect semantic types or not
 * @param  {boolean} options.readOnly If set to true, the dataframe will not be modified during the test
 * @param  {object} options List of options (optional)
 * @return {Promise<void>} The test is considered successful if it completes without errors
 */
export async function testViewer(v: string, df: DG.DataFrame,
  options?: {detectSemanticTypes?: boolean, readOnly?: boolean}): Promise<void> {
  if (options?.detectSemanticTypes) await grok.data.detectSemanticTypes(df);
  const tv = grok.shell.addTableView(df);
  const viewerName = `[name=viewer-${v.replace(/\s+/g, '-')} i]`;
  let selector = `${viewerName} canvas,${viewerName} svg,${viewerName} img,
    ${viewerName} input,${viewerName} h1,${viewerName} a`;
  const res = [];
  try {
    let viewer = tv.addViewer(v);
    await awaitCheck(() => document.querySelector(selector) !== null,
      'cannot load viewer', 3000);
    const tag = document.querySelector(selector)?.tagName;
    res.push(Array.from(tv.viewers).length);
    if (!options?.readOnly) {
      Array.from(df.row(0).cells).forEach((c:any) => c.value = null);
      const num = df.rowCount < 20 ? Math.floor(df.rowCount / 2) : 10;
      df.rows.select((row: DG.Row) => row.idx >= 0 && row.idx < num);
      await delay(50);
      for (let i = num; i < num * 2; i++) df.filter.set(i, false);
      await delay(50);
      df.currentRowIdx = 1;
      const df1 = df.clone();
      df.columns.names().slice(0, Math.ceil(df.columns.length / 2)).forEach((c: any) => df.columns.remove(c));
      await delay(100);
      tv.dataFrame = df1;
    }
    const optns = viewer.getOptions(true).look;
    const props = viewer.getProperties();
    const newProps: Record<string, string | boolean> = {};
    Object.keys(optns).filter((k) => typeof optns[k] === 'boolean').forEach((k) => newProps[k] = !optns[k]);
    props.filter((p: DG.Property) => p.choices !== null)
      .forEach((p: DG.Property) => newProps[p.name] = p.choices.find((c: any) => c !== optns[p.name])!);
    viewer.setOptions(newProps);
    await delay(300);
    const layout = tv.saveLayout();
    const oldProps = viewer.getOptions().look;
    tv.resetLayout();
    res.push(Array.from(tv.viewers).length);
    tv.loadLayout(layout);
    selector = `${viewerName} ${tag}`;
    await awaitCheck(() => document.querySelector(selector) !== null,
      'cannot load viewer from layout', 3000);
    res.push(Array.from(tv.viewers).length);
    viewer = Array.from(tv.viewers).find((v: any) => v.type !== 'Grid')!;
    expectArray(res, [2, 1, 2]);
    expect(JSON.stringify(viewer.getOptions().look), JSON.stringify(oldProps));
  } finally {
    tv.close();
    grok.shell.closeTable(df);
    DG.Balloon.closeAll();
  }
}
