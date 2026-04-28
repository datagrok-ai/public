import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-19207: drawing a line-chart annotation
// region, then closing/deleting the chart and re-opening another, threw
// "Index out of range". JS cannot synthesize the mouse drag, so we pin
// the lifecycle invariants the fix guarantees: the annotation-drawing
// helpers can be invoked on an attached Line chart, the chart can be
// closed via the attached viewer, and a fresh Line chart on the same
// TableView (or a second df.plot.line round) survives without throwing.
category('AI: GROK-19207: Annotation regions after delete', () => {
  test('enable region drawing, close attached viewer, add new line chart', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v1 = tv.addViewer(DG.VIEWER.LINE_CHART,
        {xColumnName: 'age', yColumnNames: ['height']}) as DG.LineChartViewer;
      expect(v1 instanceof DG.LineChartViewer, true);
      expect(v1.type, DG.VIEWER.LINE_CHART);
      // enableAnnotationRegionDrawing is the JS surface the fix touches.
      v1.enableAnnotationRegionDrawing();
      v1.disableAnnotationRegionDrawing();
      // Close the attached viewer (round-5: detached close may throw).
      v1.close();
      // The fix is "platform stays usable after close" — adding a new
      // Line chart to the same tv must not throw.
      var threw = false;
      var v2: DG.LineChartViewer | undefined;
      try {
        v2 = tv.addViewer(DG.VIEWER.LINE_CHART,
          {xColumnName: 'age', yColumnNames: ['weight']}) as DG.LineChartViewer;
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      expect(v2 != null, true);
      expect(v2 instanceof DG.LineChartViewer, true);
      expect(v2!.type, DG.VIEWER.LINE_CHART);
    }
    finally {
      tv.close();
    }
  });

  test('meta.annotationRegions.add survives close and re-open on the same tv', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v1 = tv.addViewer(DG.VIEWER.LINE_CHART,
        {xColumnName: 'age', yColumnNames: ['height']}) as DG.LineChartViewer;
      expect(v1 instanceof DG.LineChartViewer, true);
      const ar = v1.meta.annotationRegions;
      expect(ar != null, true);
      expect(Array.isArray(ar.items), true);
      expect(ar.items.length, 0);
      // Helper add() is the JS-reachable substitute for UI mouse drawing.
      ar.add({header: 'GROK-19207 region', description: 'lifecycle pin', opacity: 50});
      expect(ar.items.length, 1);
      expect(ar.items[0].header, 'GROK-19207 region');
      const raw = v1.props['annotationRegions'];
      expect(typeof raw, 'string');
      const parsed = JSON.parse(raw);
      expect(Array.isArray(parsed), true);
      expect(parsed.length, 1);
      // Close and re-open: the previously-buggy path.
      v1.close();
      var threw = false;
      var v2: DG.LineChartViewer | undefined;
      try {
        v2 = tv.addViewer(DG.VIEWER.LINE_CHART,
          {xColumnName: 'age', yColumnNames: ['weight']}) as DG.LineChartViewer;
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      expect(v2 instanceof DG.LineChartViewer, true);
      // Fresh chart starts with an empty annotation-regions helper.
      const ar2 = v2!.meta.annotationRegions;
      expect(ar2 != null, true);
      expect(Array.isArray(ar2.items), true);
      expect(ar2.items.length, 0);
    }
    finally {
      tv.close();
    }
  });

  test('two df.plot.line lifecycles in sequence (attached close path)', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      // First Line chart via df.plot.line; route through tv so close goes
      // through the attached path (detached close may throw per round-5).
      const v1 = df.plot.line({xColumnName: 'age', yColumnNames: ['height']}) as DG.LineChartViewer;
      expect(v1 instanceof DG.LineChartViewer, true);
      tv.dockManager.dock(v1, DG.DOCK_TYPE.RIGHT);
      v1.enableAnnotationRegionDrawing();
      v1.disableAnnotationRegionDrawing();
      v1.close();
      // Second Line chart on the same df: must not throw on creation.
      var threw = false;
      var v2: DG.LineChartViewer | undefined;
      try {
        v2 = df.plot.line({xColumnName: 'age', yColumnNames: ['weight']}) as DG.LineChartViewer;
        tv.dockManager.dock(v2, DG.DOCK_TYPE.RIGHT);
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      expect(v2 instanceof DG.LineChartViewer, true);
      expect(v2!.type, DG.VIEWER.LINE_CHART);
      if (v2 != null)
        v2.close();
    }
    finally {
      tv.close();
    }
  });
});
