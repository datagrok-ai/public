
export let tests: Test[] = [];
export let preconditions: {[key: string]: {before?: () => Promise<void>, after?: () => Promise<void>}} = {};

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
  tests.push(new Test(currentCategory, name , test));
}

export function category(category: string, tests: () => void): void {
  currentCategory = category;
  tests();
}

export function before(before: () => Promise<void>): void {
  if (preconditions[currentCategory] == undefined)
    preconditions[currentCategory] = {};
  preconditions[currentCategory].before = before;
}

export function after(after: () => Promise<void>): void {
  if (preconditions[currentCategory] == undefined)
    preconditions[currentCategory] = {};
  preconditions[currentCategory].after = after;
}


export async function runTests() {
  for (const [key, value] of Object.entries(preconditions)) {
    if (value.before)
      await value.before();
  }
  let results = tests.map(async (t) => {
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

  let data = await Promise.all(results);

  for (const [key, value] of Object.entries(preconditions)) {
    if (value.after)
      await value.after();
  }
  return data;
}