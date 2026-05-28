import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, withAttachedViewer} from '../helpers';

// DG.BoxPlot — public/js-api/src/viewer.ts (BoxPlot class); reset cast bugfix:
// core/client/xamgle/lib/src/interop/grok_api.dart (BoxPlotViewer_ResetView)
// (scenario: box-plot-annotations).
// Covers BoxPlot annotation-regions JS reach + the ResetView cast bugfix.
// Two tests document missing Dart wraps for enable/disableAnnotationRegionDrawing
// — they pass today (TypeError → threw=true) and will FAIL once the wraps land
// (intended regression signal).
category('AI: Viewers: BoxPlot Annotations', () => {
  test('resetView() fires onResetView after cast bugfix', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT,
      {value: 'age', category1: 'race'}, async (v) => {
        await expectFiresWithin(v.onResetView, () => v.resetView(), 2000);
      });
  });

  test('meta.annotationRegions on BoxPlot surfaces platform gap', async () => {
    // Confirmed gap: BoxPlotLook does not declare `annotationRegions` as a @Prop,
    // and BoxPlotCore does not compose AnnotationRegionsMixin. Accessing
    // v.meta.annotationRegions, or calling .add() on it, throws on BoxPlot
    // (works fine on ScatterPlot where the prop IS declared). This test
    // pins the gap: passes today (gap exists), FAILS once BoxPlotLook is wired.
    const v = DG.Viewer.boxPlot(demog(), {value: 'age', category1: 'race'});
    let pinned = false;
    try {
      const ar = v.meta.annotationRegions;
      ar.add({header: 'will not stick', description: 'gap pin', opacity: 50});
      // If the platform gap is fixed, ar.items.length will be 1 and this stays false.
      pinned = ar.items.length === 0;
    } catch (_e) {
      pinned = true;
    }
    expect(pinned, true);
  });

  test('enableAnnotationRegionDrawing surfaces wrap gap', async () => {
    // Dangling JS reference: viewer.ts calls grok_BoxPlotViewer_EnableAnnotationRegionDrawing
    // but no reg(...) exists in grok_api.dart -> TypeError on call. Bug-pin:
    // passes today (threw=true). Will FAIL (threw=false) once the wrap lands.
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
    // Dangling JS reference: viewer.ts calls grok_BoxPlotViewer_DisableAnnotationRegionDrawing
    // but no reg(...) exists in grok_api.dart -> TypeError on call. Bug-pin
    // symmetric to enableAnnotationRegionDrawing above.
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
        // Shrink the viewport, then assert resetView snaps back to a non-trivial box.
        v.viewport = new DG.Rect(vp.x, vp.y, vp.width / 2, vp.height / 2);
        v.resetView();
        const reset = v.viewport;
        expect(reset != null, true);
        expect(reset.width > 0, true);
        expect(reset.height > 0, true);
      });
  });
}, {owner: 'agolovko@datagrok.ai'});
