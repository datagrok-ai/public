import * as DG from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import {runTests, tests, TestContext, initAutoTests as initTests, expect, expectTable as _expectTable } from '@datagrok-libraries/test/src/test';

import './tests/funccall';

export let _package = new DG.Package();
export { tests };
import dayjs from "dayjs";

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({ category, test, testContext });
  return DG.DataFrame.fromObjects(data)!;
}

//name: initAutoTests
export async function initAutoTests() {
  await initTests(_package, _package.getModule('package-test.js'));
}

//name: expectTable
//input: dataframe actual
//input: dataframe expected
//output: bool result
export function expectTable(actual: DG.DataFrame, expected: DG.DataFrame): boolean {
  _expectTable(actual, expected);
  return true;
}

//name: expectDate
//input: datetime actual
//input: datetime expected
export function expectDate(actual: dayjs.Dayjs, expected: dayjs.Dayjs): void {
  expect(actual.valueOf(), expected.valueOf());
}
