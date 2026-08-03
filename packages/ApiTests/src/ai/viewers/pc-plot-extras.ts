import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, expectFloat, test} from '@datagrok-libraries/test/src/test';
import {demog, expectFiresWithin, expectNoThrow, subscribeAll, until, withAttachedViewer} from '../helpers';

// PcPlot layout getters, canvas/overlay handles, render side effects, and render events.
// Read-only geometry/handle/event checks share one attached+rendered plot (clear:false) to
// avoid re-attaching per test; the first-render event keeps its own fresh viewer.
category('AI: Viewers: PcPlot Extras', () => {
  const columns = ['age', 'height', 'weight'];
  let src: DG.DataFrame;
  let tv: DG.TableView;
  let v: DG.PcPlot;

  before(async () => {
    src = demog(50);
    tv = grok.shell.addTableView(src);
    v = tv.addViewer(DG.VIEWER.PC_PLOT, {columnNames: columns}) as DG.PcPlot;
    await until(() => v.chartBox != null);
  });

  after(async () => {
    tv.close();
  });

  test('resetView() round-trips and fires onResetView', async () => {
    await expectFiresWithin(v.onResetView, () => v.resetView(), 2000);
  });

  test('invalidateCanvas() is callable without throwing', async () => {
    expectNoThrow(() => v.invalidateCanvas());
  });

  test('chartBox is a Rect with sane geometry on an attached viewer', async () => {
    const box = v.chartBox;
    expect(box != null, true);
    expect(box.width > 0, true);
    expect(box.height > 0, true);
  });

  test('getColX(i) is monotonically increasing and matches chartBox endpoints', async () => {
    const box = v.chartBox;
    const xs: number[] = [];
    for (let i = 0; i < columns.length; i++)
      xs.push(v.getColX(i));
    expectFloat(xs[0], box.left);
    expectFloat(xs[columns.length - 1], box.right);
    for (let j = 1; j < xs.length; j++)
      expect(xs[j] > xs[j - 1], true);
  });

  test('activeFrame returns the source DataFrame when no transformation applied', async () => {
    expect(v.activeFrame != null, true);
    expect(v.activeFrame.rowCount, src.rowCount);
  });

  test('canvas and overlay are HTMLCanvasElement on attached viewer', async () => {
    expect(v.canvas instanceof HTMLCanvasElement, true);
    expect(v.overlay instanceof HTMLCanvasElement, true);
  });

  test('onBeforeDrawOverlay and onAfterDrawOverlay fire on render', async () => {
    await expectFiresWithin(v.onAfterDrawOverlay, () => v.invalidateCanvas(), 2000);
    await expectFiresWithin(v.onBeforeDrawOverlay, () => v.invalidateCanvas(), 2000);
  });

  // First-render event: needs a fresh viewer, so it can't share the pre-rendered fixture.
  test('onViewerRendered fires on first render', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: columns}, async (fresh) => {
      await expectFiresWithin(fresh.onViewerRendered, () => fresh.invalidateCanvas(), 2000);
    });
  });

  test('Full event sweep via subscribeAll', async () => {
    subscribeAll([v.onResetView, v.onBeforeDrawOverlay, v.onAfterDrawOverlay, v.onViewerRendered])();
  });
}, {owner: 'agolovko@datagrok.ai', clear: false});
