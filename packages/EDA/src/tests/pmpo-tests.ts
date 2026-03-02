// Tests for Probabilistic MPO (pMPO)
// Reference scores are pre-computed and stored in the 'drugs-props-train-scores.csv' file.
// This scores are computed using the library https://github.com/Merck/pmpo

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {Pmpo} from '../probabilistic-scoring/prob-scoring';
import {P_VAL_TRES_DEFAULT, Q_CUTOFF_DEFAULT, R2_DEFAULT, SCORES_PATH,
  SOURCE_PATH, EQUALITY_SIGN} from '../probabilistic-scoring/pmpo-defs';
import {getSynteticPmpoData} from '../probabilistic-scoring/data-generator';

const TIMEOUT = 10000;
const MAD_THRESH = 1E-6;

const DESIRABILITY_COL_NAME = 'CNS';
const DESCRIPTOR_NAMES = ['TPSA', 'TPSA_S', 'HBA', 'HBD', 'MW', 'nAtoms',
  'cLogD_ACD_v15', 'mapKa', 'cLogP_Biobyte', 'mbpKa', 'cLogP_ACD_v15', 'ALogP98'];
const SCORES_NAME = 'Score';
const DRUG = 'Drug';

const SIGMOIDAL = 'Sigmoidal';
const GAUSSIAN = 'Gaussian';
const PMPO_MODES = [SIGMOIDAL, GAUSSIAN];

const SAMPLES_K = 100;
const SAMPLES_COUNT = 1000 * SAMPLES_K;

/** Computes the maximum absolute deviation between pMPO scores in two data frames */
function getScoreMaxDeviation(sourceDrugCol: DG.Column, sourceScores: DG.Column,
  referenceDrugCol: DG.Column, referenceScores: DG.Column): number {
  let mad = 0;

  const sourceDrugList = sourceDrugCol.toList();
  const referenceDrugList = referenceDrugCol.toList();

  const sourceScoresRaw = sourceScores.getRawData();
  const referenceScoresRaw = referenceScores.getRawData();

  sourceDrugList.forEach((name, idx) => {
    const refIdx = referenceDrugList.indexOf(name);

    if (refIdx < 0)
      throw new Error(`Failed to compare pMPO scores: the "${name}" drug is missing in the reference data.`);

    mad = Math.max(mad, Math.abs(sourceScoresRaw[idx] - referenceScoresRaw[refIdx]));
  });

  return mad;
} // getScoreMaxDeviation

category('Probabilistic MPO: Computation', () => {
  // Correctness tests: compare pMPO scores with reference scores
  PMPO_MODES.forEach((refScoreName) => {
    const useSigmoid = (refScoreName == SIGMOIDAL);

    test('Correctness: ' + refScoreName, async () => {
      let sourceDf: DG.DataFrame | null = null;
      let referenceDf: DG.DataFrame | null = null;
      let desirability: DG.Column | null = null;
      let descriptors: DG.Column[] = [];
      let sourceDrugCol: DG.Column | null = null;
      let referenceDrugCol: DG.Column | null = null;
      let referencePrediction: DG.Column | null = null;
      let mad: number | null = null;

      try {
        // Load data
        sourceDf = await grok.dapi.files.readCsv(SOURCE_PATH);
        referenceDf = await grok.dapi.files.readCsv(SCORES_PATH);

        // Extract training items
        desirability = sourceDf.col(DESIRABILITY_COL_NAME);
        descriptors = sourceDf.columns.byNames(DESCRIPTOR_NAMES);

        if (desirability == null)
          throw new Error();

        // Train pMPO model
        const trainRes = Pmpo.fit(
          sourceDf,
          DG.DataFrame.fromColumns(descriptors).columns,
          desirability,
          P_VAL_TRES_DEFAULT,
          R2_DEFAULT,
          Q_CUTOFF_DEFAULT,
        );

        // Apply pMPO
        const prediction = Pmpo.predict(sourceDf, trainRes.params, useSigmoid, SCORES_NAME);

        // Compare with reference scores
        sourceDrugCol = sourceDf.col(DRUG);
        referenceDrugCol = referenceDf.col(DRUG);
        referencePrediction = referenceDf.col(refScoreName);

        mad = getScoreMaxDeviation(sourceDrugCol!, prediction, referenceDrugCol!, referencePrediction!);

        //console.log(refScoreName, ': max absolute deviation of pMPO scores:', mad);
      } catch (error) {
        grok.shell.error((error as Error).message);
      }

      expect(sourceDf !== null, true, 'Failed to load the source data: ' + SOURCE_PATH);
      expect(referenceDf !== null, true, 'Failed to load the scores data: ' + SCORES_PATH);
      expect(desirability !== null, true, 'Inconsistent source data: no column ' + DESIRABILITY_COL_NAME);
      expect(descriptors.length, DESCRIPTOR_NAMES.length, 'Inconsistent source data: no enough of columns');
      expect(sourceDrugCol !== null, true, 'Inconsistent source data: no column ' + DRUG);
      expect(referenceDrugCol !== null, true, 'Inconsistent reference data: no column ' + DRUG);
      expect(referencePrediction !== null, true, 'Inconsistent reference data: no column ' + SCORES_NAME);
      expect(mad !== null, true, 'Failed to compare pMPO scores with the reference data');
      expect(mad! < MAD_THRESH, true, `Max absolute deviation of pMPO scores exceeds the threshold (${MAD_THRESH})`);
    }, {timeout: TIMEOUT});
  });

  // Performance tests: measure time of pMPO training
  test('Performance: ' + SAMPLES_K + 'K drugs, ' + DESCRIPTOR_NAMES.length + ' descriptors', async () => {
    let sourceDf: DG.DataFrame | null = null;
    let desirability: DG.Column | null = null;
    let descriptors: DG.Column[] = [];

    try {
      // Generate synthetic data
      sourceDf = await getSynteticPmpoData(SAMPLES_COUNT);

      // Extract training items
      desirability = sourceDf.col(DESIRABILITY_COL_NAME);
      descriptors = sourceDf.columns.byNames(DESCRIPTOR_NAMES);

      if (desirability == null)
        throw new Error();

      // Train pMPO model
      const trainRes = Pmpo.fit(
        sourceDf,
        DG.DataFrame.fromColumns(descriptors).columns,
        desirability,
        P_VAL_TRES_DEFAULT,
        R2_DEFAULT,
        Q_CUTOFF_DEFAULT,
      );

      // Apply pMPO
      Pmpo.predict(sourceDf, trainRes.params, true, SCORES_NAME);
    } catch (error) {
      grok.shell.error((error as Error).message);
    }

    expect(sourceDf !== null, true, 'Failed to load the source data: ' + SOURCE_PATH);
    expect(desirability !== null, true, 'Inconsistent source data: no column ' + DESIRABILITY_COL_NAME);
    expect(descriptors.length, DESCRIPTOR_NAMES.length, 'Inconsistent source data: no enough of columns');
  }, {timeout: TIMEOUT});
});

/** Creates a test DataFrame with clearly separated desired/non-desired groups */
function createValidTestDf(rowCount: number = 20): DG.DataFrame {
  const half = Math.floor(rowCount / 2);
  const desList: boolean[] = [];
  const d1 = new Float64Array(rowCount);
  const d2 = new Float64Array(rowCount);
  const d3 = new Float64Array(rowCount);

  for (let i = 0; i < rowCount; i++) {
    desList.push(i < half);
    const j = i < half ? i : i - half;
    const t = j / Math.max(half - 1, 1);
    // d1, d2: clearly separated groups → low p-value
    d1[i] = i < half ? 9 + 2 * t : 1 + 2 * t;
    d2[i] = i < half ? 18 + 4 * t : 3 + 4 * t;
    // d3: same distribution in both groups → high p-value
    d3[i] = 5 + 0.1 * j;
  }

  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'des', desList),
    DG.Column.fromFloat64Array('d1', d1),
    DG.Column.fromFloat64Array('d2', d2),
    DG.Column.fromFloat64Array('d3', d3),
  ]);
}

/** Extracts a ColumnList of named descriptors from the DataFrame */
function getDescrCols(df: DG.DataFrame, names: string[]): DG.ColumnList {
  return DG.DataFrame.fromColumns(df.columns.byNames(names)).columns;
}

category('Probabilistic MPO: API', () => {
  // --- isApplicable: validates input thresholds, sample count, desirability, and descriptor quality ---

  // pValThresh = 0.0001 < P_VAL_TRES_MIN (0.001) → rejected
  test('isApplicable: rejects p-value below minimum', async () => {
    const df = createValidTestDf();
    const des = df.col('des')!;
    const descr = getDescrCols(df, ['d1', 'd2']);
    expect(Pmpo.isApplicable(descr, des, 0.0001, R2_DEFAULT, Q_CUTOFF_DEFAULT), false);
  });

  // r2Tresh = 0.001 < R2_MIN (0.01) → rejected
  test('isApplicable: rejects R² below minimum', async () => {
    const df = createValidTestDf();
    const des = df.col('des')!;
    const descr = getDescrCols(df, ['d1', 'd2']);
    expect(Pmpo.isApplicable(descr, des, P_VAL_TRES_DEFAULT, 0.001, Q_CUTOFF_DEFAULT), false);
  });

  // qCutoff = 0.001 < Q_CUTOFF_MIN (0.01) → rejected
  test('isApplicable: rejects q-cutoff below minimum', async () => {
    const df = createValidTestDf();
    const des = df.col('des')!;
    const descr = getDescrCols(df, ['d1', 'd2']);
    expect(Pmpo.isApplicable(descr, des, P_VAL_TRES_DEFAULT, R2_DEFAULT, 0.001), false);
  });

  // 8 rows < MIN_SAMPLES_COUNT (10) → rejected
  test('isApplicable: rejects too few samples', async () => {
    const df = createValidTestDf(8);
    const des = df.col('des')!;
    const descr = getDescrCols(df, ['d1', 'd2']);
    expect(Pmpo.isApplicable(descr, des, P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT), false);
  });

  // All-true desirability → stdev = 0 → rejected
  test('isApplicable: rejects single-category desirability', async () => {
    const n = 20;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'des', new Array(n).fill(true)),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => i * 2)),
    ]);
    const des = df.col('des')!;
    const descr = getDescrCols(df, ['d1', 'd2']);
    expect(Pmpo.isApplicable(descr, des, P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT), false);
  });

  // String column among descriptors → not numerical → rejected
  test('isApplicable: rejects non-numerical descriptor', async () => {
    const n = 20;
    const half = n / 2;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'des',
        Array.from({length: n}, (_, i) => i < half)),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromStrings('strCol', Array.from({length: n}, (_, i) => 'a' + i)),
    ]);
    const des = df.col('des')!;
    const descr = getDescrCols(df, ['d1', 'strCol']);
    expect(Pmpo.isApplicable(descr, des, P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT), false);
  });

  // Both descriptors are constant (stdev = 0) → no non-constant columns → rejected
  test('isApplicable: rejects all-constant descriptors', async () => {
    const n = 20;
    const half = n / 2;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'des',
        Array.from({length: n}, (_, i) => i < half)),
      DG.Column.fromFloat64Array('c1', new Float64Array(n).fill(5)),
      DG.Column.fromFloat64Array('c2', new Float64Array(n).fill(3)),
    ]);
    const des = df.col('des')!;
    const descr = getDescrCols(df, ['c1', 'c2']);
    expect(Pmpo.isApplicable(descr, des, P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT), false);
  });

  // Exactly MIN_SAMPLES_COUNT (10) rows with valid descriptors → accepted
  test('isApplicable: accepts valid data at minimum sample count', async () => {
    const df = createValidTestDf(10);
    const des = df.col('des')!;
    const descr = getDescrCols(df, ['d1', 'd2']);
    expect(Pmpo.isApplicable(descr, des, P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT), true);
  });

  // --- isTableValid: validates table structure (row count and numeric column variance) ---

  // 1 row < minimum of 2 → rejected
  test('isTableValid: rejects table with 1 row', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('a', new Float64Array([1])),
      DG.Column.fromFloat64Array('b', new Float64Array([2])),
    ]);
    expect(Pmpo.isTableValid(df, false), false);
  });

  // All columns filled with a single value → 0 non-constant columns < 2 → rejected
  test('isTableValid: rejects all-constant numeric columns', async () => {
    const n = 10;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('a', new Float64Array(n).fill(5)),
      DG.Column.fromFloat64Array('b', new Float64Array(n).fill(3)),
    ]);
    expect(Pmpo.isTableValid(df, false), false);
  });

  // Only 1 column with variance > 0, need at least 2 → rejected
  test('isTableValid: rejects single non-constant numeric column', async () => {
    const n = 10;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('a', Float64Array.from({length: n}, (_, i) => i)),
      DG.Column.fromFloat64Array('b', new Float64Array(n).fill(3)),
    ]);
    expect(Pmpo.isTableValid(df, false), false);
  });

  // Exactly 2 columns with variance > 0 → minimum met → accepted
  test('isTableValid: accepts two non-constant numeric columns', async () => {
    const n = 10;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('a', Float64Array.from({length: n}, (_, i) => i)),
      DG.Column.fromFloat64Array('b', Float64Array.from({length: n}, (_, i) => i * 2)),
    ]);
    expect(Pmpo.isTableValid(df, false), true);
  });

  // --- fit: trains pMPO model, computes statistics, filters descriptors by p-value and correlation ---

  // Valid data with well-separated groups → at least one descriptor selected
  test('fit: returns non-empty params', async () => {
    const df = createValidTestDf();
    const trainRes = Pmpo.fit(df, getDescrCols(df, ['d1', 'd2']), df.col('des')!,
      P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    expect(trainRes.params.size > 0, true, 'Expected non-empty params');
  });

  // Weights are z-scores normalized by their sum → must equal 1
  test('fit: weights sum to 1', async () => {
    const df = createValidTestDf();
    const trainRes = Pmpo.fit(df, getDescrCols(df, ['d1', 'd2']), df.col('des')!,
      P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    let sum = 0;
    trainRes.params.forEach((p) => sum += p.weight);
    expect(Math.abs(sum - 1.0) < 1e-10, true, `Weights sum ${sum} should equal 1.0`);
  });

  // Correlation filtering can only remove from p-value-selected set, never add
  test('fit: selectedByCorr is subset of selectedByPvalue', async () => {
    const df = createValidTestDf();
    const trainRes = Pmpo.fit(df, getDescrCols(df, ['d1', 'd2', 'd3']), df.col('des')!,
      P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    const allInPvalue = trainRes.selectedByCorr.every((d) => trainRes.selectedByPvalue.includes(d));
    expect(allInPvalue, true, 'selectedByCorr must be a subset of selectedByPvalue');
  });

  // Statistics table should contain one row per input descriptor (3)
  test('fit: statistics table row count matches descriptor count', async () => {
    const descrNames = ['d1', 'd2', 'd3'];
    const df = createValidTestDf();
    const trainRes = Pmpo.fit(df, getDescrCols(df, descrNames), df.col('des')!,
      P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    expect(trainRes.descrStatsTable.rowCount, descrNames.length);
  });

  // 8 rows < MIN_SAMPLES_COUNT → isApplicable fails → fit throws
  test('fit: throws on non-applicable data', async () => {
    const df = createValidTestDf(8); // too few samples
    let threw = false;
    try {
      Pmpo.fit(df, getDescrCols(df, ['d1', 'd2']), df.col('des')!,
        P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    } catch (_) {
      threw = true;
    }
    expect(threw, true, 'Expected fit to throw on non-applicable data');
  });

  // Both groups have identical distributions → t ≈ 0, p ≈ 1 → all filtered → throws
  test('fit: throws when no descriptors pass p-value filter', async () => {
    const n = 20;
    const half = n / 2;
    // Same distribution in both groups → p-value ≈ 1 → all filtered
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'des', Array.from({length: n}, (_, i) => i < half)),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => (i % half) + 1)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => ((i % half) + 1) * 2)),
    ]);
    let threw = false;
    try {
      Pmpo.fit(df, getDescrCols(df, ['d1', 'd2']), df.col('des')!,
        P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    } catch (_) {
      threw = true;
    }
    expect(threw, true, 'Expected fit to throw when no descriptors pass p-value filter');
  });

  // --- predict: applies trained pMPO parameters to compute scores ---

  // Output column length must match input DataFrame row count
  test('predict: returns column with correct length', async () => {
    const df = createValidTestDf();
    const trainRes = Pmpo.fit(df, getDescrCols(df, ['d1', 'd2']), df.col('des')!,
      P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    const prediction = Pmpo.predict(df, trainRes.params, true, SCORES_NAME);
    expect(prediction.length, df.rowCount);
  });

  // Scores = sum of weight * gaussian * sigmoid, all components >= 0
  test('predict: scores are non-negative', async () => {
    const df = createValidTestDf();
    const trainRes = Pmpo.fit(df, getDescrCols(df, ['d1', 'd2']), df.col('des')!,
      P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    const prediction = Pmpo.predict(df, trainRes.params, true, SCORES_NAME);
    const raw = prediction.getRawData();
    let allNonNeg = true;
    for (let i = 0; i < raw.length; i++)
      if (raw[i] < 0) {allNonNeg = false; break;}
    expect(allNonNeg, true, 'All scores should be non-negative');
  });

  // Weights sum to 1, gaussian in [0,1], sigmoid in [0,1] → score <= 1
  test('predict: scores do not exceed 1', async () => {
    const df = createValidTestDf();
    const trainRes = Pmpo.fit(df, getDescrCols(df, ['d1', 'd2']), df.col('des')!,
      P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    const prediction = Pmpo.predict(df, trainRes.params, true, SCORES_NAME);
    const raw = prediction.getRawData();
    let maxScore = 0;
    for (let i = 0; i < raw.length; i++) maxScore = Math.max(maxScore, raw[i]);
    expect(maxScore <= 1.0 + 1e-10, true, `Max score ${maxScore} should not exceed 1.0`);
  });

  // Sigmoid correction divides by (1 + b*c^(-dx)) → different from pure Gaussian
  test('predict: sigmoid and Gaussian modes produce different scores', async () => {
    const df = createValidTestDf();
    const trainRes = Pmpo.fit(df, getDescrCols(df, ['d1', 'd2']), df.col('des')!,
      P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    const sigScores = Pmpo.predict(df, trainRes.params, true, 'sig').getRawData();
    const gauScores = Pmpo.predict(df, trainRes.params, false, 'gau').getRawData();
    let differ = false;
    for (let i = 0; i < df.rowCount; i++)
      if (Math.abs(sigScores[i] - gauScores[i]) > 1e-12) {differ = true; break;}
    expect(differ, true, 'Sigmoid and Gaussian modes should produce different scores');
  });

  // Params reference descriptor columns not present in the target DataFrame → throws
  test('predict: throws for missing column', async () => {
    const df = createValidTestDf();
    const trainRes = Pmpo.fit(df, getDescrCols(df, ['d1', 'd2']), df.col('des')!,
      P_VAL_TRES_DEFAULT, R2_DEFAULT, Q_CUTOFF_DEFAULT);
    const incompleteDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('other', Float64Array.from({length: 5}, (_, i) => i)),
    ]);
    let threw = false;
    try {
      Pmpo.predict(incompleteDf, trainRes.params, true, SCORES_NAME);
    } catch (_) {
      threw = true;
    }
    expect(threw, true, 'Expected predict to throw for missing column');
  });
});

/** Returns default valid params for validateInputs */
function getValidInputParams(): {
  descriptors: DG.Column[] | null,
  desirability: DG.Column | null,
  threshold: number | null,
  sign: EQUALITY_SIGN,
  desirableCategories: string[] | null,
  pValue: number | null,
  r2: number | null,
  qCutoff: number | null,
} {
  const df = createValidTestDf();
  return {
    descriptors: df.columns.byNames(['d1', 'd2']),
    desirability: df.col('des')!,
    threshold: null,
    sign: EQUALITY_SIGN.DEFAULT,
    desirableCategories: null,
    pValue: P_VAL_TRES_DEFAULT,
    r2: R2_DEFAULT,
    qCutoff: Q_CUTOFF_DEFAULT,
  };
}

category('Probabilistic MPO: validateInputs', () => {
  // --- Settings validation ---

  test('validateInputs: rejects null p-value', async () => {
    const params = getValidInputParams();
    params.pValue = null;
    const result = Pmpo.validateInputs(params);
    expect(result.valid, false);
    expect(result.errors.size, 0, 'No input-specific errors for null settings');
  });

  test('validateInputs: rejects null R²', async () => {
    const params = getValidInputParams();
    params.r2 = null;
    expect(Pmpo.validateInputs(params).valid, false);
  });

  test('validateInputs: rejects null q-cutoff', async () => {
    const params = getValidInputParams();
    params.qCutoff = null;
    expect(Pmpo.validateInputs(params).valid, false);
  });

  test('validateInputs: rejects p-value out of range', async () => {
    const params = getValidInputParams();
    params.pValue = 0;
    expect(Pmpo.validateInputs(params).valid, false);
    params.pValue = 1.5;
    expect(Pmpo.validateInputs(params).valid, false);
  });

  test('validateInputs: rejects R² out of range', async () => {
    const params = getValidInputParams();
    params.r2 = -0.1;
    expect(Pmpo.validateInputs(params).valid, false);
    params.r2 = 1.5;
    expect(Pmpo.validateInputs(params).valid, false);
  });

  test('validateInputs: rejects q-cutoff out of range', async () => {
    const params = getValidInputParams();
    params.qCutoff = 0;
    expect(Pmpo.validateInputs(params).valid, false);
    params.qCutoff = 1.5;
    expect(Pmpo.validateInputs(params).valid, false);
  });

  // --- Column input validation ---

  test('validateInputs: rejects null descriptors', async () => {
    const params = getValidInputParams();
    params.descriptors = null;
    expect(Pmpo.validateInputs(params).valid, false);
  });

  test('validateInputs: rejects null desirability', async () => {
    const params = getValidInputParams();
    params.desirability = null;
    expect(Pmpo.validateInputs(params).valid, false);
  });

  test('validateInputs: rejects empty descriptors', async () => {
    const params = getValidInputParams();
    params.descriptors = [];
    const result = Pmpo.validateInputs(params);
    expect(result.valid, false);
    expect(result.errors.has('descriptors'), true);
  });

  // --- Descriptor quality validation ---

  test('validateInputs: rejects desirability among descriptors', async () => {
    const df = createValidTestDf();
    const des = df.col('des')!;
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, des],
      desirability: des,
      threshold: null,
      sign: EQUALITY_SIGN.DEFAULT,
      desirableCategories: null,
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, false);
    expect(result.errors.has('descriptors'), true);
    expect(result.errors.has('desirability'), true);
  });

  test('validateInputs: rejects zero-variance descriptors', async () => {
    const n = 20;
    const half = n / 2;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'des', Array.from({length: n}, (_, i) => i < half)),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromFloat64Array('constCol', new Float64Array(n).fill(5)),
    ]);
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, df.col('constCol')!],
      desirability: df.col('des')!,
      threshold: null,
      sign: EQUALITY_SIGN.DEFAULT,
      desirableCategories: null,
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, false);
    expect(result.errors.has('descriptors'), true);
  });

  // --- Boolean desirability validation ---

  test('validateInputs: accepts valid boolean desirability', async () => {
    const params = getValidInputParams();
    const result = Pmpo.validateInputs(params);
    expect(result.valid, true);
    expect(result.errors.size, 0);
  });

  test('validateInputs: rejects all-true boolean desirability', async () => {
    const n = 20;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'des', new Array(n).fill(true)),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => i * 2)),
    ]);
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, df.col('d2')!],
      desirability: df.col('des')!,
      threshold: null,
      sign: EQUALITY_SIGN.DEFAULT,
      desirableCategories: null,
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, false);
    expect(result.errors.has('desirability'), true);
  });

  // --- String desirability validation ---

  test('validateInputs: rejects string desirability with single category', async () => {
    const n = 20;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('des', new Array(n).fill('active')),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => i * 2)),
    ]);
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, df.col('d2')!],
      desirability: df.col('des')!,
      threshold: null,
      sign: EQUALITY_SIGN.DEFAULT,
      desirableCategories: ['active'],
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, false);
    expect(result.errors.has('desirability'), true);
  });

  test('validateInputs: rejects no selected categories', async () => {
    const n = 20;
    const half = n / 2;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('des', Array.from({length: n}, (_, i) => i < half ? 'active' : 'inactive')),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => i * 2)),
    ]);
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, df.col('d2')!],
      desirability: df.col('des')!,
      threshold: null,
      sign: EQUALITY_SIGN.DEFAULT,
      desirableCategories: [],
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, false);
    expect(result.errors.has('desirability'), true);
  });

  test('validateInputs: rejects all categories selected', async () => {
    const n = 20;
    const half = n / 2;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('des', Array.from({length: n}, (_, i) => i < half ? 'active' : 'inactive')),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => i * 2)),
    ]);
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, df.col('d2')!],
      desirability: df.col('des')!,
      threshold: null,
      sign: EQUALITY_SIGN.DEFAULT,
      desirableCategories: ['active', 'inactive'],
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, false);
    expect(result.errors.has('desirability'), true);
  });

  test('validateInputs: accepts valid string desirability', async () => {
    const n = 20;
    const half = n / 2;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('des', Array.from({length: n}, (_, i) => i < half ? 'active' : 'inactive')),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i < half ? i + 10 : i)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => i < half ? i * 3 : i)),
    ]);
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, df.col('d2')!],
      desirability: df.col('des')!,
      threshold: null,
      sign: EQUALITY_SIGN.DEFAULT,
      desirableCategories: ['active'],
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, true);
    expect(result.errors.size, 0);
  });

  // --- Numeric desirability validation ---

  test('validateInputs: rejects constant numeric desirability', async () => {
    const n = 20;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('des', new Float64Array(n).fill(5)),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => i * 2)),
    ]);
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, df.col('d2')!],
      desirability: df.col('des')!,
      threshold: 5,
      sign: EQUALITY_SIGN.DEFAULT,
      desirableCategories: null,
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, false);
    expect(result.errors.has('desirability'), true);
  });

  test('validateInputs: rejects null threshold for numeric desirability', async () => {
    const n = 20;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('des', Float64Array.from({length: n}, (_, i) => i)),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => i * 2)),
    ]);
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, df.col('d2')!],
      desirability: df.col('des')!,
      threshold: null,
      sign: EQUALITY_SIGN.DEFAULT,
      desirableCategories: null,
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, false);
    expect(result.errors.has('desirability'), true);
  });

  test('validateInputs: rejects threshold producing single group', async () => {
    const n = 20;
    // All values in [0, 19], threshold 100 with <= → all desired, none non-desired
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('des', Float64Array.from({length: n}, (_, i) => i)),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => i * 2)),
    ]);
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, df.col('d2')!],
      desirability: df.col('des')!,
      threshold: 100,
      sign: EQUALITY_SIGN.LESS_OR_EQUAL,
      desirableCategories: null,
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, false);
    expect(result.errors.has('desirability'), true);
    expect(result.errors.has('threshold'), true);
  });

  test('validateInputs: accepts valid numeric desirability with threshold', async () => {
    const n = 20;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('des', Float64Array.from({length: n}, (_, i) => i)),
      DG.Column.fromFloat64Array('d1', Float64Array.from({length: n}, (_, i) => i + 1)),
      DG.Column.fromFloat64Array('d2', Float64Array.from({length: n}, (_, i) => i * 2)),
    ]);
    const result = Pmpo.validateInputs({
      descriptors: [df.col('d1')!, df.col('d2')!],
      desirability: df.col('des')!,
      threshold: 10,
      sign: EQUALITY_SIGN.LESS_OR_EQUAL,
      desirableCategories: null,
      pValue: P_VAL_TRES_DEFAULT,
      r2: R2_DEFAULT,
      qCutoff: Q_CUTOFF_DEFAULT,
    });
    expect(result.valid, true);
    expect(result.errors.size, 0);
  });
});
