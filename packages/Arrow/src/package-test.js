import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import {runTests, tests} from "@datagrok-libraries/utils/src/test";
import './tests/arrow-tests.ts'

export let _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//output: dataframe result
export async function test(category, test) {
  const data = await runTests({category, test});
  return DG.DataFrame.fromObjects(data);
}
