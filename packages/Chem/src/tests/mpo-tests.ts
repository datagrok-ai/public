import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {mpo} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {NumericalDesirability, QualifierPolicy, MpoScale, MissingValueConfig, DesirabilityLine}
  from '@datagrok-libraries/statistics/src/mpo/mpo-types';

/* Computation tests for MPO desirability scoring — what the BDD features (Chem/bdd/features/mpo) must not
   assert, since it is the numbers an analysis computes (the lead's ruling; MPO MISSING.md item 6). The
   scoring reads each column's `desirabilityTemplate` tag and folds the per-column desirabilities into a
   score; qualified numbers (<x, >x) score by the Score as policy (worst / limit / best). The ADME-like
   matrix numbers are the same ones Chem/bdd/features/mpo/qnum-scoring.feature claims end to end. */
category('MPO: qnum scoring', () => {
  const numerical = (p: Partial<NumericalDesirability> & {line: DesirabilityLine}): NumericalDesirability => ({
    functionType: 'numerical', weight: 1, missingValues: {strategy: 'exclude'}, ...p});

  // fromList rejects the qnum type, so qnum columns go through DG.Column.qnum with exact=false (the values are
  // already packed by DG.Qnum.less/greater/exact; null is read as an empty cell).
  const qnum = (name: string, values: (number | null)[], t: NumericalDesirability): DG.Column => {
    const col = DG.Column.qnum(name, values.length, values as number[], false);
    col.setTag('desirabilityTemplate', JSON.stringify(t));
    return col;
  };
  const num = (name: string, values: number[], t: NumericalDesirability): DG.Column => {
    const col = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values);
    col.setTag('desirabilityTemplate', JSON.stringify(t));
    return col;
  };

  /** Scores the Caco2/Lipophilicity/Solubility table of qnum-scoring.feature with the given Caco2 policy and
   *  missing-value strategy. Lipophilicity 2 and Solubility 0 are weight 1 and desirability 1 (flat curve),
   *  Caco2 weighs 0.8 on the curve -8→0, -6→1, -4→0, so a row scores (2 + 0.8 d) / 2.8, as the feature derives. */
  function scoreCaco2(policy: QualifierPolicy, missing: MissingValueConfig = {strategy: 'exclude'}): DG.Column {
    const caco2 = qnum('Caco2',
      [DG.Qnum.less(-5), DG.Qnum.greater(-5), DG.Qnum.exact(-6), DG.Qnum.less(-9), null],
      numerical({weight: 0.8, line: [[-8, 0], [-6, 1], [-4, 0]], qualifiers: policy, missingValues: missing}));
    const flat1: DesirabilityLine = [[-1000, 1], [1000, 1]];
    const lipo = num('Lipophilicity', [2, 2, 2, 2, 2], numerical({line: flat1}));
    const sol = num('Solubility', [0, 0, 0, 0, 0], numerical({line: flat1}));
    const df = DG.DataFrame.fromColumns([caco2, lipo, sol]);
    return mpo(df, [caco2, lipo, sol], 'ADME', 'Average').scoreColumn;
  }

  test('Least favorable value (the default) scores a bound by the worst desirability it allows', async () => {
    const s = scoreCaco2(QualifierPolicy.Worst);
    expectFloat(s.get(0), 0.714286, 1e-4, '<-5');   // worst over [-8,-5] is 0
    expectFloat(s.get(1), 0.714286, 1e-4, '>-5');   // worst over [-5,-4] is 0
    expectFloat(s.get(2), 1.0, 1e-4, '-6 exact');   // on the peak
    expectFloat(s.get(3), 0.714286, 1e-4, '<-9');   // wholly left of the curve, 0
    expect(s.isNone(4), true, 'empty Caco2 excludes the row');
  });

  test('Limit value scores a bound as if it were exact', async () => {
    const s = scoreCaco2(QualifierPolicy.Limit);
    expectFloat(s.get(0), 0.857143, 1e-4, '<-5 at -5 is 0.5');
    expectFloat(s.get(1), 0.857143, 1e-4, '>-5 at -5 is 0.5');
    expectFloat(s.get(2), 1.0, 1e-4, '-6 exact');
    expectFloat(s.get(3), 0.714286, 1e-4, '<-9 still 0 (limit off the curve clamps to 0)');
  });

  test('Most favorable value scores a bound by the best desirability it allows', async () => {
    const s = scoreCaco2(QualifierPolicy.Best);
    expectFloat(s.get(0), 1.0, 1e-4, '<-5 allows up to 1');
    expectFloat(s.get(1), 0.857143, 1e-4, '>-5 allows up to 0.5');
    expectFloat(s.get(2), 1.0, 1e-4, '-6 exact');
    expectFloat(s.get(3), 0.714286, 1e-4, '<-9 still 0');
  });

  test('Missing value: exclude bails the row, default substitutes, skip drops the property', async () => {
    expect(scoreCaco2(QualifierPolicy.Worst, {strategy: 'exclude'}).isNone(4), true, 'exclude → empty');
    // default 0.5 for the missing Caco2: (2 + 0.8*0.5) / 2.8
    expectFloat(scoreCaco2(QualifierPolicy.Worst, {strategy: 'default', score: 0.5}).get(4), 0.857143, 1e-4, 'default');
    // skip drops Caco2 from row 5 only: (1 + 1) / (1 + 1)
    expectFloat(scoreCaco2(QualifierPolicy.Worst, {strategy: 'skip'}).get(4), 1.0, 1e-4, 'skip');
  });

  test('An inverted curve scores 1 - d(x)', async () => {
    // line rises 0→1 over 0→10; inverted, so x=2 (d 0.2) scores 0.8, x=8 (d 0.8) scores 0.2
    const col = num('x', [2, 8, 5],
      numerical({line: [[0, 0], [10, 1]], inverted: true}));
    const s = mpo(DG.DataFrame.fromColumns([col]), [col], 'inv', 'Average').scoreColumn;
    expectFloat(s.get(0), 0.8, 1e-4, 'x=2 inverted');
    expectFloat(s.get(1), 0.2, 1e-4, 'x=8 inverted');
    expectFloat(s.get(2), 0.5, 1e-4, 'x=5 inverted');
  });

  test('A log-scale curve: positive exacts interpolate in log space, <=0 and bounds at <=0 are empty', async () => {
    // raw line 1→0, 100→1, read in log10 space (xs 0..2); 10 → log10 10 = 1 → d 0.5
    const col = qnum('IC50',
      [DG.Qnum.exact(10), DG.Qnum.exact(0), DG.Qnum.less(0), DG.Qnum.exact(-2), DG.Qnum.greater(0)],
      numerical({line: [[1, 0], [100, 1]], scale: MpoScale.Log, qualifiers: QualifierPolicy.Worst}));
    const s = mpo(DG.DataFrame.fromColumns([col]), [col], 'log', 'Average').scoreColumn;
    expectFloat(s.get(0), 0.5, 1e-4, '10 → log10 is 1, halfway');
    expect(s.isNone(1), true, 'exact 0 has no log');
    expect(s.isNone(2), true, '<0 has no log');
    expect(s.isNone(3), true, 'exact -2 has no log');
    expectFloat(s.get(4), 0.0, 1e-4, '>0 spans the whole curve, worst is 0');
  });

  test('A bound and an exact value wholly past the end of the curve both score 0', async () => {
    // Lipophilicity curve ends at 5 with 0.5; 6 and >6 are past the end
    const col = qnum('Lipophilicity',
      [DG.Qnum.exact(6), DG.Qnum.greater(6), DG.Qnum.exact(3)],
      numerical({line: [[-2, 0], [2, 1], [5, 0.5]], qualifiers: QualifierPolicy.Limit}));
    const s = mpo(DG.DataFrame.fromColumns([col]), [col], 'past', 'Average').scoreColumn;
    expectFloat(s.get(0), 0.0, 1e-4, 'exact 6 is past the end');
    expectFloat(s.get(1), 0.0, 1e-4, '>6 is wholly past the end');
    expect(s.get(2) > 0, true, '3 is on the curve');
  });

  test('The smallest exact value on a flat curve scores the flat value, like the rest', async () => {
    // a flat curve at 0.5 over the column's range; the minimum exact value is on the edge, not off it
    const col = qnum('IC50',
      [DG.Qnum.exact(0.5), DG.Qnum.exact(2), DG.Qnum.greater(10), DG.Qnum.exact(5)],
      numerical({line: [[0.5, 0.5], [10, 0.5]], qualifiers: QualifierPolicy.Worst}));
    const s = mpo(DG.DataFrame.fromColumns([col]), [col], 'flat', 'Average').scoreColumn;
    expectFloat(s.get(0), 0.5, 1e-4, 'smallest exact on the edge');
    expectFloat(s.get(1), 0.5, 1e-4);
    expectFloat(s.get(3), 0.5, 1e-4);
  });
});

category('MPO: profile handler', () => {
  // other packages call Chem:init as well (HitTriage does), and each call used to register another handler
  test('Calling Chem:init again does not register another profile handler', async () => {
    const handlers = () => DG.ObjectHandler.list().filter((h) => h.type === 'mpo-profile').length;
    await grok.functions.call('Chem:init', {});
    const before = handlers();
    await grok.functions.call('Chem:init', {});
    await grok.functions.call('Chem:init', {});
    expect(handlers(), before);
    expect(before > 0, true, 'no profile handler registered');
  });
});
