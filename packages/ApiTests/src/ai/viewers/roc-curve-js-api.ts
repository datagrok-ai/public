import * as DG from 'datagrok-api/dg';
import {category, expect, test, expectFloat} from '@datagrok-libraries/test/src/test';
import {df, expectRoundTripPropAndLook, expectNoThrow, withAttachedViewer} from '../helpers';

// DG.RocCurve
// AUC is reachable only via the static RocCurveCore.calculateImpl -> RocCurveCalculationResult.auc; v.auc() wraps it.

// Build a binary target + numeric score frame.
function rocDf(target: string[], score: number[]): DG.DataFrame {
  return df([['actual', DG.COLUMN_TYPE.STRING, target], ['score', DG.COLUMN_TYPE.FLOAT, score]]);
}

const ATTACH = {targetColumn: 'actual', positiveClass: 'yes', predictionColumnNames: ['score']};

category('AI: Viewers: Roc Curve', () => {
  test('rocCurve factory returns typed DG.RocCurve of the right type', async () => {
    const v = DG.Viewer.rocCurve(rocDf(['yes', 'yes', 'no', 'no'], [0.9, 0.8, 0.2, 0.1]), ATTACH);
    expect(v instanceof DG.RocCurve, true);
    expect(v.type, DG.VIEWER.ROC_CURVE);
  });

  test('settings round-trip (targetColumn/positiveClass/showThreshold)', async () => {
    const v = DG.Viewer.rocCurve(rocDf(['yes', 'yes', 'no', 'no'], [0.9, 0.8, 0.2, 0.1]), ATTACH);
    // predictionColumnNames is list-valued and can't be equality-checked; it's exercised via the
    // factory ATTACH options and AUC cases below, so round-trip only the scalar props here.
    expectRoundTripPropAndLook(v, {targetColumn: 'actual', positiveClass: 'yes', showThreshold: true});
  });

  test('perfect separation -> AUC == 1.0', async () => {
    const d = rocDf(['yes', 'yes', 'no', 'no'], [0.9, 0.8, 0.2, 0.1]);
    const v = DG.Viewer.rocCurve(d, ATTACH);
    expectFloat(v.auc(d.col('actual')!, d.col('score')!, 'yes'), 1.0);
  });

  test('inverted separation -> AUC == 0.0', async () => {
    const d = rocDf(['yes', 'yes', 'no', 'no'], [0.1, 0.2, 0.8, 0.9]);
    const v = DG.Viewer.rocCurve(d, ATTACH);
    expectFloat(v.auc(d.col('actual')!, d.col('score')!, 'yes'), 0.0);
  });

  test('partial separation -> deterministic AUC == 0.75', async () => {
    // Desc-sorted scores: 0.9(yes),0.6(no),0.3(yes),0.1(no) => one positive outranks the negatives,
    // one does not. tpr=[0,.5,.5,1,1], fpr=[0,0,.5,.5,1] => trapezoidal AUC = 0.75.
    const d = rocDf(['yes', 'yes', 'no', 'no'], [0.9, 0.3, 0.6, 0.1]);
    const v = DG.Viewer.rocCurve(d, ATTACH);
    expectFloat(v.auc(d.col('actual')!, d.col('score')!, 'yes'), 0.75);
  });

  test('tied scores -> AUC is a valid bounded value', async () => {
    // All-equal scores: AUC depends on input row order (no separating power), so only assert bounds.
    const d = rocDf(['yes', 'no', 'yes', 'no'], [0.5, 0.5, 0.5, 0.5]);
    const v = DG.Viewer.rocCurve(d, ATTACH);
    const auc = v.auc(d.col('actual')!, d.col('score')!, 'yes');
    expect(auc >= 0 && auc <= 1, true);
  });

  test('boundary: attaches and degenerate single-class AUC does not throw', async () => {
    const d = rocDf(['yes', 'yes', 'yes', 'yes'], [0.9, 0.8, 0.2, 0.1]);
    await withAttachedViewer<DG.RocCurve>(d, DG.VIEWER.ROC_CURVE, ATTACH, async (v) => {
      // Single-class frame: totalNegatives == 0 => fpr divides by zero, AUC may be NaN. Scope to no-throw.
      expectNoThrow(() => v.auc(d.col('actual')!, d.col('score')!, 'yes'));
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
