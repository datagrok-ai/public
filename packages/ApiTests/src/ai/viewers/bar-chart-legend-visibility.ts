import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, wait, withAttachedViewer} from '../helpers';

// Source: core/client/d4/lib/src/viewers/bar_chart/bar_chart_core.dart:298
// (canShowLegend) and :69 (legendCol => stackCol when stackCol is
// categorical/datetime). When stackColumnName == splitColumnName, the
// legend would duplicate the categorical X axis, so the auto legend is
// suppressed. Logic: canShowLegend() returns false iff
//   legendVisibility == Auto && sameLegendColumnDisplayed(splitCol)
// where sameLegendColumnDisplayed(col) checks col == legendCol == stackCol.
// Test relies on real attach + repaint: detached viewers skip legend layout
// (see _refreshLegend early-return on dataFrame == null), so we use a real
// TableView and a 300ms delay to let the deferred legend refresh land.
// Same-column case: .d4-legend should be removed from the DOM
// (see legend.remove() at legend_mixin.dart:530 and legend.dart:168
// 'root?.remove()'). Different-column case: .d4-legend should be present.
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
        await wait();
        expect(legendPresent(v), true);
      });
  });

  test('legend toggles when stack switches to/from split column', async () => {
    await withAttachedViewer<DG.BarChartViewer>(demog(), DG.VIEWER.BAR_CHART,
      {value: 'age', split: 'race', stack: 'sex'}, async (v) => {
        await wait();
        expect(legendPresent(v), true);
        v.setOptions({stackColumnName: 'race'});
        await wait();
        expect(legendPresent(v), false);
        v.setOptions({stackColumnName: 'sex'});
        await wait();
        expect(legendPresent(v), true);
      });
  });
}, {owner: 'agolovko@datagrok.ai'});
