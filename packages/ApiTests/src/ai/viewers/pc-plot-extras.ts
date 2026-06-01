import * as DG from 'datagrok-api/dg';
import {category, expect, expectFloat, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, expectNoThrow, subscribeAll, until, withAttachedViewer} from '../helpers';

// PcPlot layout getters, canvas/overlay handles, render side effects, and render events.
category('AI: Viewers: PcPlot Extras', () => {
  const columns = ['age', 'height', 'weight'];

  test('resetView() round-trips and fires onResetView', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: columns}, async (v) => {
      await until(() => v.chartBox != null);
      await expectFiresWithin(v.onResetView, () => v.resetView(), 2000);
    });
  });

  test('invalidateCanvas() is callable without throwing', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: columns}, async (v) => {
      await until(() => v.chartBox != null);
      expectNoThrow(() => v.invalidateCanvas());
    });
  });

  test('chartBox is a Rect with sane geometry on an attached viewer', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: columns}, async (v) => {
      await until(() => v.chartBox != null);
      const box = v.chartBox;
      expect(box != null, true);
      expect(box.width > 0, true);
      expect(box.height > 0, true);
    });
  });

  test('chartW equals chartBox.width / (columnNames.length - 1)', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: columns}, async (v) => {
      await until(() => v.chartBox != null);
      const expected = v.chartBox.width / (columns.length - 1);
      expectFloat(v.chartW, expected);
    });
  });

  test('getColX(i) is monotonically increasing and matches chartBox endpoints', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: columns}, async (v) => {
      await until(() => v.chartBox != null);
      const box = v.chartBox;
      const xs: number[] = [];
      for (let i = 0; i < columns.length; i++)
        xs.push(v.getColX(i));
      expectFloat(xs[0], box.left);
      expectFloat(xs[columns.length - 1], box.right);
      for (let j = 1; j < xs.length; j++)
        expect(xs[j] > xs[j - 1], true);
    });
  });

  test('activeFrame returns the source DataFrame when no transformation applied', async () => {
    const df = demog(50);
    await withAttachedViewer<DG.PcPlot>(df, DG.VIEWER.PC_PLOT, {columnNames: columns}, async (v) => {
      await until(() => v.chartBox != null);
      expect(v.activeFrame != null, true);
      expect(v.activeFrame.rowCount, df.rowCount);
    });
  });

  test('canvas and overlay are HTMLCanvasElement on attached viewer', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: columns}, async (v) => {
      await until(() => v.chartBox != null);
      expect(v.canvas instanceof HTMLCanvasElement, true);
      expect(v.overlay instanceof HTMLCanvasElement, true);
    });
  });

  test('onBeforeDrawOverlay and onAfterDrawOverlay fire on render', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: columns}, async (v) => {
      await until(() => v.chartBox != null);
      await expectFiresWithin(v.onAfterDrawOverlay, () => v.invalidateCanvas(), 2000);
      await expectFiresWithin(v.onBeforeDrawOverlay, () => v.invalidateCanvas(), 2000);
    });
  });

  test('onViewerRendered fires on first render', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: columns}, async (v) => {
      await expectFiresWithin(v.onViewerRendered, () => v.invalidateCanvas(), 2000);
    });
  });

  test('Full event sweep via subscribeAll', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: columns}, async (v) => {
      await until(() => v.chartBox != null);
      subscribeAll([v.onResetView, v.onBeforeDrawOverlay, v.onAfterDrawOverlay, v.onViewerRendered])();
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
