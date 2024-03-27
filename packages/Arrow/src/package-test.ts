import * as DG from 'datagrok-api/dg';
import {runTests, tests, TestContext} from '@datagrok-libraries/utils/src/test';
import {default as init} from "parquet-wasm/esm/arrow1";

import './tests/arrow-tests';

export const _package = new DG.Package();
export {tests}

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  await init(_package.webRoot + 'dist/arrow1_bg.wasm');
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}