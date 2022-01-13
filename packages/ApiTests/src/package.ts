/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './dataframe/data_frame';
import './dataframe/calculated-columns';
import './shell/shell';
import './shell/windows';
import './viewer/viewer';
import './views/events';
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
import './ui/buttons';
import './widgets/legend';

import {runTests} from "@datagrok-libraries/utils/src/test";
export let _package = new DG.Package();


//name: testJsApi
//output: dataframe result
//top-menu: Tools | Dev | JS API Tests
export async function testJsApi(): Promise<DG.DataFrame> {
  let data = await runTests();
  return DG.DataFrame.fromObjects(data)!;
}
