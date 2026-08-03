import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test, expectFloat} from '@datagrok-libraries/test/src/test';
import {df, expectRoundTripPropAndLook} from '../helpers';

// DG.RocCurve
// auc() takes the target/prediction columns and positive class directly and computes the AUC from
// that data alone (api.grok_RocCurve_CalculateAuc) — it does not read the viewer's own state. So a
// single viewer instance can be reused to evaluate every AUC scenario below.

// Build a binary target + numeric score frame.
function rocDf(target: string[], score: number[]): DG.DataFrame {
  return df([['actual', DG.COLUMN_TYPE.STRING, target], ['score', DG.COLUMN_TYPE.FLOAT, score]]);
}

const ATTACH = {targetColumn: 'actual', positiveClass: 'yes', predictionColumnNames: ['score']};

function auc(v: DG.RocCurve, target: string[], score: number[]): number {
  const d = rocDf(target, score);
  return v.auc(d.col('actual')!, d.col('score')!, 'yes');
}

category('AI: Viewers: Roc Curve', () => {
  let v: DG.RocCurve;

  before(async () => {
    v = DG.Viewer.rocCurve(rocDf(['yes', 'yes', 'no', 'no'], [0.9, 0.8, 0.2, 0.1]), ATTACH);
  });

  after(async () => {
    v?.detach();
  });

  test('rocCurve factory returns typed DG.RocCurve of the right type', async () => {
    expect(v instanceof DG.RocCurve, true);
    expect(v.type, DG.VIEWER.ROC_CURVE);
  });

  test('settings round-trip (targetColumn/positiveClass/showThreshold)', async () => {
    // predictionColumnNames is list-valued and can't be equality-checked; it's exercised via the
    // factory ATTACH options and AUC cases below, so round-trip only the scalar props here.
    expectRoundTripPropAndLook(v, {targetColumn: 'actual', positiveClass: 'yes', showThreshold: true});
  });

  test('perfect separation -> AUC == 1.0', async () => {
    expectFloat(auc(v, ['yes', 'yes', 'no', 'no'], [0.9, 0.8, 0.2, 0.1]), 1.0);
  });

  test('inverted separation -> AUC == 0.0', async () => {
    expectFloat(auc(v, ['yes', 'yes', 'no', 'no'], [0.1, 0.2, 0.8, 0.9]), 0.0);
  });

  test('partial separation -> deterministic AUC == 0.75', async () => {
    // Desc-sorted scores: 0.9(yes),0.6(no),0.3(yes),0.1(no) => one positive outranks the negatives,
    // one does not. tpr=[0,.5,.5,1,1], fpr=[0,0,.5,.5,1] => trapezoidal AUC = 0.75.
    expectFloat(auc(v, ['yes', 'yes', 'no', 'no'], [0.9, 0.3, 0.6, 0.1]), 0.75);
  });

  test('tied scores -> AUC is a valid bounded value', async () => {
    // All-equal scores have no separating power; the AUC depends on input row order, so assert bounds.
    const a = auc(v, ['yes', 'no', 'yes', 'no'], [0.5, 0.5, 0.5, 0.5]);
    expect(a >= 0 && a <= 1, true);
  });
}, {owner: 'agolovko@datagrok.ai'});
