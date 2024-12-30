import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {delay, testEvent} from '@datagrok-libraries/utils/src/test';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {_package} from '../package-test';


export async function initHelmMainPackage(): Promise<void> {
  await getHelmHelper();
}

export async function awaitGrid(grid: DG.Grid, timeout: number = 5000): Promise<void> {
  await delay(100);
  await testEvent(grid.onAfterDrawContent, () => {},
    () => { grid.invalidate(); }, timeout);
}

export async function readDataframe(tableName: string): Promise<DG.DataFrame> {
  const file = await loadFileAsText(tableName);
  const df = DG.DataFrame.fromCsv(file);
  df.name = tableName.replace('.csv', '');
  return df;
}

export async function loadFileAsText(name: string): Promise<string> {
  return await _package.files.readAsText(name);
}
