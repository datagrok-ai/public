import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, withAttachedViewer} from '../helpers';

// BoxPlot annotation-regions JS reach plus the ResetView cast bugfix.
category('AI: Viewers: BoxPlot Annotations', () => {
  test('resetView() fires onResetView after cast bugfix', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT,
      {value: 'age', category1: 'race'}, async (v) => {
        await expectFiresWithin(v.onResetView, () => v.resetView(), 2000);
      });
  });

  test('meta.annotationRegions on BoxPlot surfaces platform gap', async () => {
    // Pins a known gap: passes today, will FAIL once BoxPlotLook declares annotationRegions.
    const v = DG.Viewer.boxPlot(demog(), {value: 'age', category1: 'race'});
    let pinned = false;
    try {
      const ar = v.meta.annotationRegions;
      ar.add({header: 'will not stick', description: 'gap pin', opacity: 50});
      pinned = ar.items.length === 0;
    } catch (_e) {
      pinned = true;
    }
    expect(pinned, true);
  });

  test('enableAnnotationRegionDrawing surfaces wrap gap', async () => {
    // Bug-pin: passes today (TypeError), will FAIL once the Dart wrap lands.
    const v = DG.Viewer.boxPlot(demog(), {value: 'age', category1: 'race'});
    let threw = false;
    try {
      (v as any).enableAnnotationRegionDrawing((_r: any) => {});
    } catch (_e) {
      threw = true;
    }
    expect(threw, true);
  });

  test('disableAnnotationRegionDrawing surfaces wrap gap', async () => {
    // Bug-pin: passes today (TypeError), will FAIL once the Dart wrap lands.
    const v = DG.Viewer.boxPlot(demog(), {value: 'age', category1: 'race'});
    let threw = false;
    try {
      (v as any).disableAnnotationRegionDrawing();
    } catch (_e) {
      threw = true;
    }
    expect(threw, true);
  });

  test('viewport getter/setter survives resetView()', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT,
      {value: 'age', category1: 'race'}, async (v) => {
        const vp = v.viewport;
        expect(vp != null, true);
        v.viewport = new DG.Rect(vp.x, vp.y, vp.width / 2, vp.height / 2);
        v.resetView();
        const reset = v.viewport;
        expect(reset != null, true);
        expect(reset.width > 0, true);
        expect(reset.height > 0, true);
      });
  });
}, {owner: 'agolovko@datagrok.ai'});
