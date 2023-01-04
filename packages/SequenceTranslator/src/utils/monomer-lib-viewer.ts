/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {LIB_PATH, DEFAULT_LIB_FILENAME} from './const';

export async function viewMonomerLib(): Promise<void> {
  const table = await parseMonomerLib(LIB_PATH, DEFAULT_LIB_FILENAME);
  grok.shell.addTableView(table);
}

async function parseMonomerLib(path: string, fileName: string): Promise<DG.DataFrame> {
  const fileSource = new DG.FileSource(path);
  const file = await fileSource.readAsText(fileName);
  const obj = JSON.parse(file);
  const df = DG.DataFrame.fromObjects(obj)!;
  return df;
}
