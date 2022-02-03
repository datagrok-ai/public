import * as grok from "datagrok-api/grok";

export let tests: {[key: string]: {tests?: Test[], before?: () => Promise<void>, after?: () => Promise<void>, beforeStatus?: string, afterStatus?: string}} = {};

export let currentCategory: string;

export namespace assure {

  export function notNull(value: any, name?: string) {
    if (value == null)
      throw `${name == null ? 'Value' : name} not defined`;
  }
}

export class Test {
  test: () => Promise<any>;
  name: string;
  category: string;

  constructor(category: string, name: string, test: () => Promise<any>) {
    this.category = category;
    this.name = name;
    this.test = test;
  }
}

export function test(name: string, test: () => Promise<any>): void {
  if (tests[currentCategory] == undefined)
    tests[currentCategory] = {};
  if (tests[currentCategory].tests == undefined)
    tests[currentCategory].tests = [];
  tests[currentCategory].tests!.push(new Test(currentCategory, name , test));
}

/** Awaits for a while checking the error status of the platform */
export async function testAssureNoError(ms: number): Promise<boolean> {
  const timer = (ms: number) => new Promise(res => setTimeout(res, ms));
  for (let i = 0; i < ms / 500; ++i) {
    await timer(500);
    if (grok.shell.lastError.length !== 0)
      return false;
  }
  return true;
}

/** Does best effort to catch all exceptional situations. Currently doesn't catch "Uncaught in promise" */
export function testExpectFinish(name: string, foo: () => Promise<any>, ms: number = 5000): void {
  test(name, async (): Promise<any> => {
    try {
      grok.shell.lastError = '';
      await foo();
    } catch (e: any) {
      throw "Exception is not expected";
    }
    let noError = await testAssureNoError(ms);
    if (!noError)
      throw "Exceptions while waiting for results";
  });
}

/** Tests two objects for equality, throws an exception if they are not equal. */
export function expect(actual: any, expected: any): void {
  if (actual !== expected)
    throw `Expected "${expected}", got "${actual}"`;
}

export function expectFloat(actual: number, expected: number, tolerance = 0.001): void {
  const areEqual = Math.abs(actual - expected) < 0.001;
  if (!areEqual)
    throw `Expected ${expected}, got ${actual} (tolerance = ${tolerance})"`;
}

/** Defines a test suite. */
export function category(category: string, tests: () => void): void {
  currentCategory = category;
  tests();
}

/** Defines a function to be executed before the tests in this category are executed. */
export function before(before: () => Promise<void>): void {
  if (tests[currentCategory] == undefined)
    tests[currentCategory] = {};
  tests[currentCategory].before = before;
}

/** Defines a function to be executed after the tests in this category are executed. */
export function after(after: () => Promise<void>): void {
  if (tests[currentCategory] == undefined)
    tests[currentCategory] = {};
  tests[currentCategory].after = after;
}


export async function runTests(options?: {category?: string, test?: string}) {
  let results: { category?: string, name?: string, success: boolean, result: string}[] = [];

  for (const [key, value] of Object.entries(tests)) {
    if (options?.category != undefined) {
      if (!key.toLowerCase().startsWith(options?.category.toLowerCase()))
        continue;
    }
    try {
      if (value.before)
        await value.before();
    } catch (x: any) {
      value.beforeStatus = x.toString();
    }
    let t = value.tests ?? [];
    let res = [];
    for(let i = 0; i < t.length; i++) {
      res.push(await execTest(t[i], options?.test));
    }

    let data = (await Promise.all(res)).filter((d) => d.result != 'skipped');
    try {
      if (value.after)
        await value.after();
    } catch (x: any) {
      value.afterStatus = x.toString();
    }
    if (value.afterStatus)
      data.push({category: key, name: 'init', result: value.afterStatus, success: false});
    if (value.beforeStatus)
      data.push({category: key, name: 'init', result: value.beforeStatus, success: false});
    results.push(...data);
  }

  return results;
}

async function execTest(t: Test, predicate: string | undefined) {
  let r: { category?: string, name?: string, success: boolean, result: string };
  try {
    if (predicate != undefined && (!t.name.toLowerCase().startsWith(predicate.toLowerCase())))
      r = {success: true, result: 'skipped'};
    else
      r = {success: true, result: await t.test() ?? 'OK'};
  } catch (x: any) {
    r = {success: false, result: x.toString()};
  }
  r.category = t.category;
  r.name = t.name;
  return r;
}

export async function delay(ms: number) {
  await  new Promise(r => setTimeout(r, ms));
}
