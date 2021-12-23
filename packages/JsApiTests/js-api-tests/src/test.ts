
export let tests: Test[] = [];
export let currentCategory: String;

export class Test {
  test: Function;
  name: String;
  category: String;

  constructor(category: String, name: String, test: Function) {
    this.category = category;
    this.name = name;
    this.test = test;
  }
}

export function test(name: String, test: Function): void {
  tests.push(new Test(currentCategory, name , test));
}

export function category(category: String, tests: Function): void {
  currentCategory = category;
  tests();
}