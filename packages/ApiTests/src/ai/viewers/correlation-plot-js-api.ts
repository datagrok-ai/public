import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, expectFloat, before, after} from '@datagrok-libraries/test/src/test';
import {df, expectNoThrow, withAttachedViewer, until} from '../helpers';

// DG.CorrelationPlot: getCorrelation honors the viewer's combined filter, so attach before reading.

// x=[1..5]; yPos = 2*x (perfect +1 Pearson); yNeg = 6-x (perfect -1 Pearson).
function corrDf(): DG.DataFrame {
  return df([
    ['x', DG.COLUMN_TYPE.FLOAT, [1, 2, 3, 4, 5]],
    ['yPos', DG.COLUMN_TYPE.FLOAT, [2, 4, 6, 8, 10]],
    ['yNeg', DG.COLUMN_TYPE.FLOAT, [5, 4, 3, 2, 1]],
  ]);
}

const ATTACH = {xColumnNames: ['x', 'yPos', 'yNeg'], yColumnNames: ['x', 'yPos', 'yNeg']};

category('AI: Viewers: Correlation Plot', () => {
  let tv: DG.TableView; let v: DG.CorrelationPlot; let d: DG.DataFrame;
  before(async () => {
    d = corrDf();
    tv = grok.shell.addTableView(d);
    v = tv.addViewer(DG.VIEWER.CORR_PLOT, ATTACH) as DG.CorrelationPlot;
    await until(() => v != null);
  });
  after(async () => {
    tv.close();
  });

  test('correlationPlot factory returns typed DG.CorrelationPlot of the right type', async () => {
    const v = DG.Viewer.correlationPlot(corrDf(), ATTACH);
    expect(v instanceof DG.CorrelationPlot, true);
    expect(v.type, DG.VIEWER.CORR_PLOT);
  });

  test('perfect positive linear relation -> Pearson r == 1.0', async () => {
    expectFloat(v.getCorrelation(d.col('x')!, d.col('yPos')!), 1.0);
  });

  test('perfect negative linear relation -> Pearson r == -1.0', async () => {
    expectFloat(v.getCorrelation(d.col('x')!, d.col('yNeg')!), -1.0);
  });

  test('correlation of a column with itself == 1.0', async () => {
    expectFloat(v.getCorrelation(d.col('x')!, d.col('x')!), 1.0);
  });

  test('correlation is symmetric in its arguments', async () => {
    const ab = v.getCorrelation(d.col('x')!, d.col('yNeg')!);
    const ba = v.getCorrelation(d.col('yNeg')!, d.col('x')!);
    expectFloat(ab, ba);
  });

  test('boundary: single-row frame attaches and reading a coefficient does not throw', async () => {
    const one = df([['x', DG.COLUMN_TYPE.FLOAT, [1.0]], ['y', DG.COLUMN_TYPE.FLOAT, [2.0]]]);
    await withAttachedViewer<DG.CorrelationPlot>(one, DG.VIEWER.CORR_PLOT,
      {xColumnNames: ['x', 'y'], yColumnNames: ['x', 'y']}, async (v) => {
        // Single row => undefined variance; coefficient may be NaN. Scope to no-throw.
        expectNoThrow(() => v.getCorrelation(one.col('x')!, one.col('y')!));
      });
  });
}, {owner: 'agolovko@datagrok.ai', clear: false});
