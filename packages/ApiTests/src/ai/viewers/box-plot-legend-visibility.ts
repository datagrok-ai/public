import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, look, until, withAttachedViewer} from '../helpers';

// BoxPlot auto legend is suppressed when the single category column equals markerColorColumnName.
// BoxPlot auto-syncs markerColorColumnName to the category column; tests flip allowColorSynchronization
// off before assigning a different markerColor so the deferred sync respects the override.
category('AI: Viewers: BoxPlot legend visibility', () => {
  // legendVisibility must be Auto for the "single category == markerColor" suppression to apply;
  // 'Always' would force the legend on and contradict the hidden-legend assertions below.
  const opts = {valueColumnName: 'age', category1ColumnName: 'race',
    markerColorColumnName: 'race', legendVisibility: 'Auto'};
  const legendPresent = (v: DG.Viewer): boolean => v.root.querySelector('.d4-legend') != null;

  test('legend hidden when single category equals markerColor column', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT, opts, async (v) => {
      await until(() => look(v)['markerColorColumnName'] === 'race');
      expectLook(v, {category1ColumnName: 'race', markerColorColumnName: 'race'});
      expect(legendPresent(v), false);
    });
  });

  test('legend visible when category and markerColor columns differ', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT, opts, async (v) => {
      await until(() => look(v)['category1ColumnName'] === 'race');
      v.setOptions({allowColorSynchronization: false, markerColorColumnName: 'sex'});
      await until(() => look(v)['markerColorColumnName'] === 'sex' && legendPresent(v));
      expectLook(v, {category1ColumnName: 'race', markerColorColumnName: 'sex'});
      expect(legendPresent(v), true);
    });
  });

  test('legend toggles when markerColor switches to/from category column', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT, opts, async (v) => {
      await until(() => look(v)['category1ColumnName'] === 'race');
      v.setOptions({allowColorSynchronization: false, markerColorColumnName: 'sex'});
      await until(() => legendPresent(v) === true);
      expect(legendPresent(v), true);
      v.setOptions({markerColorColumnName: 'race'});
      await until(() => legendPresent(v) === false);
      expect(legendPresent(v), false);
      v.setOptions({markerColorColumnName: 'sex'});
      await until(() => legendPresent(v) === true);
      expect(legendPresent(v), true);
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
