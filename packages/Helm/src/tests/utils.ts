import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {delay, testEvent} from '@datagrok-libraries/utils/src/test';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';


export async function initHelmMainPackage(): Promise<void> {
  await getHelmHelper();
}

export async function awaitGrid(grid: DG.Grid, timeout: number = 5000): Promise<void> {
  await delay(100);
  await testEvent(grid.onAfterDrawContent, () => {},
    () => { grid.invalidate(); }, timeout);
}
