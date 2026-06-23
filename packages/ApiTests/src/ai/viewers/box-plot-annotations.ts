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

  test('meta.annotationRegions round-trips on BoxPlot', async () => {
    const v = DG.Viewer.boxPlot(demog(), {value: 'age', category1: 'race'});
    const ar = v.meta.annotationRegions;
    const before = ar.items.length;
    ar.add({header: 'AI band', description: 'box annotation', opacity: 50});
    expect(ar.items.length, before + 1);
    expect(ar.items[ar.items.length - 1].header, 'AI band');
  });

  test('enableAnnotationRegionDrawing is callable without throwing', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT,
      {value: 'age', category1: 'race'}, async (v) => {
        expect(typeof (v as any).enableAnnotationRegionDrawing, 'function');
        (v as any).enableAnnotationRegionDrawing(false, (_r: any) => {});
        (v as any).disableAnnotationRegionDrawing();
      });
  });

  test('disableAnnotationRegionDrawing is callable without throwing', async () => {
    await withAttachedViewer<DG.BoxPlot>(demog(), DG.VIEWER.BOX_PLOT,
      {value: 'age', category1: 'race'}, async (v) => {
        expect(typeof (v as any).disableAnnotationRegionDrawing, 'function');
        (v as any).disableAnnotationRegionDrawing();
      });
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
