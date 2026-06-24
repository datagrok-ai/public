import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, until, wait, withAttachedViewer} from '../helpers';

// BarChart auto legend is suppressed when stack column equals split column (would duplicate the X axis).
category('AI: Viewers: BarChart legend visibility', () => {
  const legendPresent = (v: DG.Viewer): boolean => v.root.querySelector('.d4-legend') != null;

  test('legend hidden when stack column equals split column', async () => {
    await withAttachedViewer<DG.BarChartViewer>(demog(), DG.VIEWER.BAR_CHART,
      {value: 'age', split: 'race', stack: 'race'}, async (v) => {
        expectLook(v, {splitColumnName: 'race', stackColumnName: 'race'});
        await wait();
        expect(legendPresent(v), false);
      });
  });

  test('legend visible when stack and split columns differ', async () => {
    await withAttachedViewer<DG.BarChartViewer>(demog(), DG.VIEWER.BAR_CHART,
      {value: 'age', split: 'race', stack: 'sex'}, async (v) => {
        expectLook(v, {splitColumnName: 'race', stackColumnName: 'sex'});
        await until(() => legendPresent(v));
        expect(legendPresent(v), true);
      });
  });

  test('legend toggles when stack switches to/from split column', async () => {
    await withAttachedViewer<DG.BarChartViewer>(demog(), DG.VIEWER.BAR_CHART,
      {value: 'age', split: 'race', stack: 'sex'}, async (v) => {
        await until(() => legendPresent(v));
        expect(legendPresent(v), true);
        v.setOptions({stackColumnName: 'race'});
        await wait();
        expect(legendPresent(v), false);
        v.setOptions({stackColumnName: 'sex'});
        await until(() => legendPresent(v));
        expect(legendPresent(v), true);
      });
  });
}, {owner: 'agolovko@datagrok.ai'});
