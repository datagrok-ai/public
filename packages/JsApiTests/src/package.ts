/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './dataframe/df';
import './shell/shell';
import './viewer/viewer';
import {runTests} from "./test";

export let _package = new DG.Package();


//name: test
export async function test() {
  let data = await runTests();
  grok.shell.addTableView(DG.DataFrame.fromObjects(data)!);
}
