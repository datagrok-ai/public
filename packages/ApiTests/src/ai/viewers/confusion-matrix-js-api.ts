import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, expectFloat, expectArray, before, after} from '@datagrok-libraries/test/src/test';
import {df, expectRoundTripPropAndLook, expectNoThrow, withAttachedViewer, until} from '../helpers';

// DG.ConfusionMatrix
// Stats compute on frame attach (_refreshValues). Always attach + await categories before asserting.

// Build an actual/predicted DataFrame from a count map keyed by `${actual};${predicted}`.
function matrixDf(counts: {[pair: string]: number}): DG.DataFrame {
  const actual: string[] = [];
  const predicted: string[] = [];
  for (const pair of Object.keys(counts)) {
    const [a, p] = pair.split(';');
    for (let i = 0; i < counts[pair]; i++) {
      actual.push(a);
      predicted.push(p);
    }
  }
  return df([['actual', DG.COLUMN_TYPE.STRING, actual], ['prediction', DG.COLUMN_TYPE.STRING, predicted]]);
}

// Known 2x2 matrix: actual 'yes' is positive, 'no' is negative.
// (yes;yes)=TP=3 (yes;no)=FN=2 (no;yes)=FP=1 (no;no)=TN=4 ; total=10
const BIN = {'yes;yes': 3, 'yes;no': 2, 'no;yes': 1, 'no;no': 4};
const ATTACH = {xColumnName: 'actual', yColumnName: 'prediction'};

// Row share for class c where actual==c — diagonal(c) / (sum over predicted of count[c;*]).
function rowShare(c: string): number {
  let num = 0; let den = 0;
  for (const pair of Object.keys(BIN)) {
    const [a, p] = pair.split(';');
    if (a === c) {
      den += BIN[pair as keyof typeof BIN];
      if (p === c)
        num += BIN[pair as keyof typeof BIN];
    }
  }
  return num / den;
}

// Column share for class c where predicted==c — diagonal(c) / (sum over actual of count[*;c]).
function colShare(c: string): number {
  let num = 0; let den = 0;
  for (const pair of Object.keys(BIN)) {
    const [a, p] = pair.split(';');
    if (p === c) {
      den += BIN[pair as keyof typeof BIN];
      if (a === c)
        num += BIN[pair as keyof typeof BIN];
    }
  }
  return num / den;
}

category('AI: Viewers: Confusion Matrix', () => {
  let tv: DG.TableView; let v: DG.ConfusionMatrix; let d: DG.DataFrame;
  before(async () => {
    d = matrixDf(BIN);
    tv = grok.shell.addTableView(d);
    v = tv.addViewer(DG.VIEWER.CONFUSION_MATRIX, ATTACH) as DG.ConfusionMatrix;
    await until(() => v != null && v.categories.length > 0);
  });
  after(async () => {
    tv.close();
  });

  test('confusionMatrix factory returns typed DG.ConfusionMatrix of the right type', async () => {
    const v = DG.Viewer.confusionMatrix(matrixDf(BIN), ATTACH);
    expect(v instanceof DG.ConfusionMatrix, true);
    expect(v.type, DG.VIEWER.CONFUSION_MATRIX);
  });

  test('x/y column props round-trip', async () => {
    const v = DG.Viewer.confusionMatrix(matrixDf(BIN), ATTACH);
    expectRoundTripPropAndLook(v, {xColumnName: 'actual', yColumnName: 'prediction'});
  });

  test('binary frame: two categories and isBinary true', async () => {
    expect(v.categories.length, 2);
    expect(v.isBinary, true);
    expectArray([...v.categories].sort(), ['no', 'yes']);
  });

  test('accuracy is exact on a known matrix', async () => {
    expectFloat(v.accuracy, 0.7); // (TP+TN)/total = (3+4)/10
  });

  test('sensitivity/specificity/precision/NPV exact on a known matrix', async () => {
    // Index 0 is the "positive" class; derive expectations from the live category order.
    const pos = v.categories[0];
    const neg = v.categories[1];
    expectFloat(v.sensitivity, rowShare(pos));
    expectFloat(v.specificity, rowShare(neg));
    expectFloat(v.precision, colShare(pos));
    expectFloat(v.negativePredictedValue, colShare(neg));
  });

  test('getRowShare and getColumnShare per index match hand-computed shares', async () => {
    expectFloat(v.getRowShare(0), rowShare(v.categories[0]));
    expectFloat(v.getRowShare(1), rowShare(v.categories[1]));
    expectFloat(v.getColumnShare(0), colShare(v.categories[0]));
    expectFloat(v.getColumnShare(1), colShare(v.categories[1]));
  });

  test('non-binary 3-class frame: isBinary false, three categories, binary stats degrade to 0', async () => {
    const three = matrixDf({'a;a': 2, 'b;b': 2, 'c;c': 2, 'a;b': 1, 'b;c': 1});
    await withAttachedViewer<DG.ConfusionMatrix>(three, DG.VIEWER.CONFUSION_MATRIX, ATTACH, async (v) => {
      await until(() => v.categories.length > 0);
      expect(v.categories.length, 3);
      expect(v.isBinary, false);
      expectFloat(v.sensitivity, 0.0); // core returns 0.0 when !isBinary
    });
  });

  test('boundary: single-row frame attaches without throwing', async () => {
    const one = df([['actual', DG.COLUMN_TYPE.STRING, ['yes']], ['prediction', DG.COLUMN_TYPE.STRING, ['yes']]]);
    await withAttachedViewer<DG.ConfusionMatrix>(one, DG.VIEWER.CONFUSION_MATRIX, ATTACH, async (v) => {
      await until(() => v.categories.length > 0);
      expect(v.categories.length >= 1, true);
      // Single class => not binary; reading stats must not throw (shares may be null on empty off-diagonal).
      expectNoThrow(() => v.accuracy);
      expectNoThrow(() => v.getColumnShare(0));
    });
  });
}, {owner: 'agolovko@datagrok.ai', clear: false});
