
export let tests: {[key: string]: {tests?: Test[], before?: () => Promise<void>, after?: () => Promise<void>, beforeStatus?: string, afterStatus?: string}} = {};

export let currentCategory: string;

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

export function expect(a: any, b: any) {
  if (a != b)
    throw `Value "${a}" is not equal. Expected "${b}"`;
}

export function category(category: string, tests: () => void): void {
  currentCategory = category;
  tests();
}

export function before(before: () => Promise<void>): void {
  if (tests[currentCategory] == undefined)
    tests[currentCategory] = {};
  tests[currentCategory].before = before;
}

export function after(after: () => Promise<void>): void {
  if (tests[currentCategory] == undefined)
    tests[currentCategory] = {};
  tests[currentCategory].after = after;
}


export async function runTests() {
  let results: { category?: string, name?: string, success: boolean, result: string}[] = [];

  for (const [key, value] of Object.entries(tests)) {
    try {
      if (value.before)
        await value.before();
    } catch (x: any) {
       value.beforeStatus = x.toString();
    }
    let t = value.tests ?? [];
    let res = t.map(async (t) => {
      let r: { category?: string, name?: string, success: boolean, result: string };
      try {
        r = {success: true, result: await t.test() ?? 'OK'};
      } catch (x: any) {
        r = {success: false, result: x.toString()};
      }
      r.category = t.category;
      r.name = t.name;
      return r;
    });

    let data = await Promise.all(res);
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

export async function delay(ms: number) {
  await  new Promise(r => setTimeout(r, ms));
}
