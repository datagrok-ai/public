import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {LIB_PATH} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

import {MonomerLibFileManager} from './lib-file-manager';

export async function getLibFileNameList(): Promise<string[]> {
  const fileManager = await MonomerLibFileManager.getInstance();
  return fileManager.getValidFiles();
}
