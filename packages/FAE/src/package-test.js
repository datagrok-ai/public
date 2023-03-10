import * as DG from "datagrok-api/dg";
import {runTests} from "@datagrok-libraries/utils/src/test";

export let _package = new DG.Package();

//name: test
//output: dataframe result
export async function test() {
  let data = await runTests();
  return DG.DataFrame.fromObjects(data);
}
