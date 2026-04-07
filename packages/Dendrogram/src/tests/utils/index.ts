import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {delay, testEvent} from '@datagrok-libraries/test/src/test';

export async function awaitGrid(grid: DG.Grid, timeout: number = 5000): Promise<void> {
  await delay(0);
  await testEvent(grid.onAfterDrawContent, () => {},
    () => { grid.invalidate(); }, timeout);
}
