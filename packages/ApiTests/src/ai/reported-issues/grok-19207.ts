import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, withTableView} from '../helpers';

// Regression coverage for GROK-19207: drawing a line-chart annotation
// region, then closing/deleting the chart and re-opening another, threw
// "Index out of range". JS cannot synthesize the mouse drag, so we pin
// the lifecycle invariants the fix guarantees: the annotation-drawing
// helpers can be invoked on an attached Line chart, the chart can be
// closed via the attached viewer, and a fresh Line chart on the same
// TableView (or a second df.plot.line round) survives without throwing.
category('AI: GROK-19207: Annotation regions after delete', () => {
  test('enable region drawing, close attached viewer, add new line chart', async () => {
    await withTableView(demog(), async (tv) => {
      const v1 = tv.addViewer(DG.VIEWER.LINE_CHART,
        {xColumnName: 'age', yColumnNames: ['height']}) as DG.LineChartViewer;
      v1.enableAnnotationRegionDrawing();
      v1.disableAnnotationRegionDrawing();
      v1.close();
      let v2: DG.LineChartViewer | undefined;
      expectNoThrow(() => {
        v2 = tv.addViewer(DG.VIEWER.LINE_CHART, {xColumnName: 'age', yColumnNames: ['weight']}) as DG.LineChartViewer;
      });
      expect(v2 instanceof DG.LineChartViewer, true);
      expect(v2!.type, DG.VIEWER.LINE_CHART);
    });
  });

  test('meta.annotationRegions.add survives close and re-open on the same tv', async () => {
    await withTableView(demog(), async (tv) => {
      const v1 = tv.addViewer(DG.VIEWER.LINE_CHART,
        {xColumnName: 'age', yColumnNames: ['height']}) as DG.LineChartViewer;
      const ar = v1.meta.annotationRegions;
      expect(ar.items.length, 0);
      ar.add({header: 'GROK-19207 region', description: 'lifecycle pin', opacity: 50});
      expect(ar.items.length, 1);
      const parsed = JSON.parse(v1.props['annotationRegions']);
      expect(parsed.length, 1);
      v1.close();
      let v2: DG.LineChartViewer | undefined;
      expectNoThrow(() => {
        v2 = tv.addViewer(DG.VIEWER.LINE_CHART, {xColumnName: 'age', yColumnNames: ['weight']}) as DG.LineChartViewer;
      });
      expect(v2!.meta.annotationRegions.items.length, 0);
    });
  });

  test('two df.plot.line lifecycles in sequence (attached close path)', async () => {
    await withTableView(demog(), async (tv) => {
      const df = tv.dataFrame;
      const v1 = df.plot.line({xColumnName: 'age', yColumnNames: ['height']}) as DG.LineChartViewer;
      tv.dockManager.dock(v1, DG.DOCK_TYPE.RIGHT);
      v1.enableAnnotationRegionDrawing();
      v1.disableAnnotationRegionDrawing();
      v1.close();
      let v2: DG.LineChartViewer | undefined;
      expectNoThrow(() => {
        v2 = df.plot.line({xColumnName: 'age', yColumnNames: ['weight']}) as DG.LineChartViewer;
        tv.dockManager.dock(v2!, DG.DOCK_TYPE.RIGHT);
      });
      expect(v2 instanceof DG.LineChartViewer, true);
      v2!.close();
    });
  });
});
