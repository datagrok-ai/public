/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './dataframe/data_frame';
import './dataframe/calculated-columns';
import './shell/shell';
import './shell/windows';
import './viewer/viewer';
import './views/layouts';
import './dapi/files';
import './dapi/fetch';
import './dapi/admin';
import './dapi/groups';
import './ui/inputs';
import './dapi/dapi';
import './dapi/entities';
import './dapi/layouts';
import './dapi/projects';
import './dapi/tables';
import './dapi/user-data-storage';
import './dapi/users';
import './shell/ml';
import './ui/divs';

import {runTests} from "./test";
export let _package = new DG.Package();


//name: testJsApi
//top-menu: Tools | Dev | JS API Tests
export async function testJsApi() {
  let data = await runTests();
  grok.shell.addTableView(DG.DataFrame.fromObjects(data)!);
}
