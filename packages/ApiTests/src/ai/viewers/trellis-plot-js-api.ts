import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {df, expectRoundTripPropAndLook, expectNoThrow, withAttachedViewer, until} from '../helpers';

// DG.TrellisPlotViewer: x/yCategoriesCount are the product of resolved X/Y column category counts,
// null until columns resolve on attach. oneColumnOnly reads look only (no attach needed).

// g2 has 2 categories (a,b); g3 has 3 (x,y,z); val feeds the inner viewer.
function trellisDf(): DG.DataFrame {
  return df([
    ['g2', DG.COLUMN_TYPE.STRING, ['a', 'a', 'b', 'b', 'a', 'b']],
    ['g3', DG.COLUMN_TYPE.STRING, ['x', 'y', 'z', 'x', 'y', 'z']],
    ['val', DG.COLUMN_TYPE.FLOAT, [1, 2, 3, 4, 5, 6]],
  ]);
}

const ATTACH = {xColumnNames: ['g2'], yColumnNames: ['g3'], viewerType: DG.VIEWER.HISTOGRAM};

category('AI: Viewers: Trellis Plot', () => {
  test('trellisPlot factory returns typed DG.TrellisPlotViewer of the right type', async () => {
    const v = DG.Viewer.trellisPlot(trellisDf(), ATTACH);
    expect(v instanceof DG.TrellisPlotViewer, true);
    expect(v.type, DG.VIEWER.TRELLIS_PLOT);
  });

  test('viewerType prop round-trips', async () => {
    const v = DG.Viewer.trellisPlot(trellisDf(), ATTACH);
    expectRoundTripPropAndLook(v, {viewerType: DG.VIEWER.BAR_CHART});
  });

  test('oneColumnOnly is true with a single category axis, false with both', async () => {
    const single = DG.Viewer.trellisPlot(trellisDf(), {xColumnNames: ['g2'], yColumnNames: []}) as DG.TrellisPlotViewer;
    expect(single.oneColumnOnly, true);
    const both = DG.Viewer.trellisPlot(trellisDf(), ATTACH) as DG.TrellisPlotViewer;
    expect(both.oneColumnOnly, false);
  });

  test('x/y category counts equal the product of column category counts', async () => {
    await withAttachedViewer<DG.TrellisPlotViewer>(trellisDf(), DG.VIEWER.TRELLIS_PLOT, ATTACH, async (v) => {
      await until(() => v.xCategoriesCount != null && v.yCategoriesCount != null);
      expect(v.xCategoriesCount, 2); // g2 => {a,b}
      expect(v.yCategoriesCount, 3); // g3 => {x,y,z}
    });
  });

  test('two X columns multiply the X category count', async () => {
    const opts = {xColumnNames: ['g2', 'g3'], yColumnNames: ['g3'], viewerType: DG.VIEWER.HISTOGRAM};
    await withAttachedViewer<DG.TrellisPlotViewer>(trellisDf(), DG.VIEWER.TRELLIS_PLOT, opts, async (v) => {
      await until(() => v.xCategoriesCount != null);
      expect(v.xCategoriesCount, 6); // |g2| * |g3| = 2 * 3
    });
  });

  test('boundary: counts read null-safe before attach, no throw', async () => {
    const v = DG.Viewer.trellisPlot(trellisDf(), ATTACH) as DG.TrellisPlotViewer;
    // Columns may be unresolved pre-attach: the getter must degrade to null, not throw.
    expectNoThrow(() => v.xCategoriesCount);
    expectNoThrow(() => v.yCategoriesCount);
  });
}, {owner: 'agolovko@datagrok.ai'});
