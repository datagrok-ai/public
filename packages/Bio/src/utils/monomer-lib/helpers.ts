import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {LIB_PATH} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

import {_package} from '../../package';

export async function getLibFileNameList(): Promise<string[]> {
  // list files recursively because permissions are available for folders only
  const res: string[] = await Promise.all((await grok.dapi.files.list(LIB_PATH, true, ''))
    .map(async (it) => {
      // Get relative path (to LIB_PATH)
      return it.fullPath.substring(LIB_PATH.length);
    }));
  return res;
}

export async function manageFiles() {
  const a = ui.dialog({title: 'Manage files'})
    .add(ui.fileBrowser({path: 'System:AppData/Bio/libraries'}).root)
    .addButton('OK', () => a.close())
    .show();
}

