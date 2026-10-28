import * as DG from 'datagrok-api/dg';
import {ViewHost} from '@datagrok-libraries/u2/dg';

export function campaignRibbon(tv: DG.TableView, session: {refresh(): Promise<void>}): void {
  ViewHost.of(tv).ribbon({
    name: 'campaign-actions',
    actions: [
      {icon: 'plus', tooltip: 'Add new row', run: () => tv.dataFrame.rows.addNew(null, true)},
      {icon: 'sync', tooltip: 'Refresh Campaign', run: () => session.refresh()},
    ],
  });
}
