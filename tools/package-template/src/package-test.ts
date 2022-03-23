import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import {runTests} from "@datagrok-libraries/utils/src/test";
import "./tests/chem-tests";

export let _package = new DG.Package();

//name: test
//output: dataframe result
export async function test(): Promise<DG.DataFrame> {
  let data = await runTests();
  return DG.DataFrame.fromObjects(data)!;
}