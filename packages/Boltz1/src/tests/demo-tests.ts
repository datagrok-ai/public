import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {take} from 'rxjs/operators';

import {openBoltzDemo} from '../demo/boltz-app';

category('Demo', () => {
  test('opensComplexCellWhenAnotherViewIsCurrent', async () => {
    grok.shell.closeAll();
    grok.shell.addTableView(DG.DataFrame.fromColumns([DG.Column.fromStrings('Complex', ['decoy'])]));
    const sub = grok.events.onProjectOpened.pipe(take(1))
      .subscribe(() => grok.shell.v = grok.shell.addView(DG.View.create()));
    try {
      await openBoltzDemo();
      await delay(100);
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.name === 'boltz_demo_data');
      expect(tv != null, true);
      expect(grok.shell.v === tv, true);
      expect(tv!.dataFrame.currentCell?.column?.name, 'Complex');
    }
    finally {
      sub.unsubscribe();
      grok.shell.closeAll();
    }
  }, {timeout: 60000});
});
