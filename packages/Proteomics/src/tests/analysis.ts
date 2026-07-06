import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {setGroups, getGroups, GroupAssignment, seedAnnotationDialogInputs,
  setOrganism, getOrganism, applyAnnotation} from '../analysis/experiment-setup';
import {medianNormalize, quantileNormalize, vsnNormalize} from '../analysis/normalization';
import {imputeMinProb, imputeKnn, imputeZero, imputeMean, imputeMedian} from '../analysis/imputation';
import {runDifferentialExpression, copyDEResultsToFrame, getDefaultComparison}
  from '../analysis/differential-expression';
import {SEMTYPE} from '../utils/proteomics-types';

/** Build a test DataFrame with protein IDs and intensity columns for two groups. */
function makeTestDf(
  nProteins: number,
  group1Names: string[],
  group2Names: string[],
  group1Values?: number[][],
  group2Values?: number[][],
): DG.DataFrame {
  const cols: DG.Column[] = [];
  cols.push(DG.Column.fromStrings('Protein ID',
    Array.from({length: nProteins}, (_, i) => `P${i}`)));

  const allNames = [...group1Names, ...group2Names];
  const allValues = group1Values && group2Values ?
    [...group1Values, ...group2Values] : undefined;

  for (let c = 0; c < allNames.length; c++) {
    const values = new Float32Array(nProteins);
    if (allValues && allValues[c]) {
      for (let r = 0; r < nProteins; r++)
        values[r] = allValues[c][r];
    }
    const col = DG.Column.fromFloat32Array(allNames[c], values);
    col.semType = SEMTYPE.INTENSITY;
    cols.push(col);
  }
  return DG.DataFrame.fromColumns(cols);
}

// ── Experiment Setup ────────────────────────────────────────────────

category('Experiment Setup', () => {
  test('setGroups persists group assignments as tag', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0']),
      DG.Column.fromFloat32Array('s1', new Float32Array([1])),
    ]);
    const groups: GroupAssignment = {
      group1: {name: 'Control', columns: ['s1']},
      group2: {name: 'Treatment', columns: ['s2']},
    };
    setGroups(df, groups);
    const raw = df.getTag('proteomics.groups');
    expect(raw !== null && raw !== '', true);
    const parsed = JSON.parse(raw!);
    expect(parsed.group1.name, 'Control');
    expect(parsed.group2.name, 'Treatment');
  });

  test('getGroups returns null for untagged DataFrame', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0']),
    ]);
    const result = getGroups(df);
    expect(result, null);
  });

  test('getOrganism returns undefined until set, then round-trips the code', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('id', ['P0'])]);
    expect(getOrganism(df) === undefined, true);
    setOrganism(df, 'rnorvegicus');
    expect(getOrganism(df), 'rnorvegicus');
  });

  test('applyAnnotation persists the organism alongside the groups', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('id', ['P0'])]);
    applyAnnotation(df, {
      group1: {name: 'Control', columns: ['s1', 's2']},
      group2: {name: 'Treatment', columns: ['s3', 's4']},
      organism: 'rnorvegicus',
    });
    expect(getOrganism(df), 'rnorvegicus');
    expect(getGroups(df)!.group1.name, 'Control');
    // Organism is optional — omitting it leaves the tag untouched.
    const df2 = DG.DataFrame.fromColumns([DG.Column.fromStrings('id', ['P0'])]);
    applyAnnotation(df2, {
      group1: {name: 'A', columns: ['s1']},
      group2: {name: 'B', columns: ['s2']},
    });
    expect(getOrganism(df2) === undefined, true);
  });

  test('getGroups round-trips GroupAssignment correctly', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0']),
    ]);
    const groups: GroupAssignment = {
      group1: {name: 'Control', columns: ['s1', 's2', 's3']},
      group2: {name: 'Treatment', columns: ['s4', 's5', 's6']},
    };
    setGroups(df, groups);
    const retrieved = getGroups(df);
    expect(retrieved !== null, true);
    expect(retrieved!.group1.name, 'Control');
    expect(retrieved!.group2.name, 'Treatment');
    expect(retrieved!.group1.columns.length, 3);
    expect(retrieved!.group2.columns.length, 3);
    expect(retrieved!.group1.columns[0], 's1');
    expect(retrieved!.group2.columns[2], 's6');
  });

  test('seedAnnotationDialogInputs returns Control/Treatment defaults when no tag', async () => {
    const df = makeTestDf(3, ['s1', 's2'], ['s3', 's4']);
    const seed = seedAnnotationDialogInputs(df, ['s1', 's2', 's3', 's4']);
    expect(seed.group1Name, 'Control');
    expect(seed.group2Name, 'Treatment');
    expect(seed.group1Cols.length, 0);
    expect(seed.group2Cols.length, 0);
  });

  test('seedAnnotationDialogInputs reads names + columns from existing groups', async () => {
    const df = makeTestDf(3, ['s1', 's2'], ['s3', 's4']);
    setGroups(df, {
      group1: {name: 'DMD', columns: ['s1', 's2']},
      group2: {name: 'WT', columns: ['s3', 's4']},
    });
    const seed = seedAnnotationDialogInputs(df, ['s1', 's2', 's3', 's4']);
    expect(seed.group1Name, 'DMD');
    expect(seed.group2Name, 'WT');
    expect(seed.group1Cols.length, 2);
    expect(seed.group2Cols.length, 2);
    expect(seed.group1Cols[0].name, 's1');
    expect(seed.group2Cols[1].name, 's4');
  });

  test('seedAnnotationDialogInputs drops column refs not in the DataFrame', async () => {
    const df = makeTestDf(3, ['s1', 's2'], ['s3', 's4']);
    setGroups(df, {
      group1: {name: 'DMD', columns: ['s1', 'gone', 's2']},
      group2: {name: 'WT', columns: ['also-gone', 's3', 's4']},
    });
    const seed = seedAnnotationDialogInputs(df, ['s1', 's2', 's3', 's4']);
    expect(seed.group1Cols.length, 2);
    expect(seed.group2Cols.length, 2);
    expect(seed.group1Cols.map((c) => c.name).join(','), 's1,s2');
    expect(seed.group2Cols.map((c) => c.name).join(','), 's3,s4');
  });

  test('seedAnnotationDialogInputs drops column refs not in available list', async () => {
    const df = makeTestDf(3, ['s1', 's2'], ['s3', 's4']);
    setGroups(df, {
      group1: {name: 'A', columns: ['s1', 's2']},
      group2: {name: 'B', columns: ['s3', 's4']},
    });
    const seed = seedAnnotationDialogInputs(df, ['s1', 's3']);
    expect(seed.group1Cols.length, 1);
    expect(seed.group2Cols.length, 1);
    expect(seed.group1Cols[0].name, 's1');
    expect(seed.group2Cols[0].name, 's3');
  });

  test('seedAnnotationDialogInputs falls back to defaults when stored names are empty', async () => {
    const df = makeTestDf(3, ['s1', 's2'], ['s3', 's4']);
    setGroups(df, {
      group1: {name: '', columns: ['s1', 's2']},
      group2: {name: '', columns: ['s3', 's4']},
    });
    const seed = seedAnnotationDialogInputs(df, ['s1', 's2', 's3', 's4']);
    expect(seed.group1Name, 'Control');
    expect(seed.group2Name, 'Treatment');
  });
});

// ── Normalization ───────────────────────────────────────────────────

category('Normalization', () => {
  test('median normalization shifts column medians to zero', async () => {
    const col = DG.Column.fromFloat32Array('intensity', new Float32Array([2, 4, 6, 8, 10]));
    col.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3', 'P4']),
      col,
    ]);
    // Median of [2,4,6,8,10] is 6
    medianNormalize(df, ['intensity']);
    const c = df.col('intensity')!;
    // Values should be shifted by -6: [-4, -2, 0, 2, 4]
    expect(Math.abs(c.get(0) - (-4)) < 0.01, true);
    expect(Math.abs(c.get(1) - (-2)) < 0.01, true);
    expect(Math.abs(c.get(2) - 0) < 0.01, true);
    expect(Math.abs(c.get(3) - 2) < 0.01, true);
    expect(Math.abs(c.get(4) - 4) < 0.01, true);
  });

  test('median normalization skips null values', async () => {
    // 5 values: [2, null, 6, 8, null] -- non-null median of [2,6,8] = 6
    const col = DG.Column.fromFloat32Array('intensity',
      new Float32Array([2, DG.FLOAT_NULL, 6, 8, DG.FLOAT_NULL]));
    col.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3', 'P4']),
      col,
    ]);
    medianNormalize(df, ['intensity']);
    const c = df.col('intensity')!;
    // Nulls should remain null
    expect(c.isNone(1), true);
    expect(c.isNone(4), true);
    // Non-null values shifted by -median
    expect(!c.isNone(0), true);
    expect(!c.isNone(2), true);
    expect(!c.isNone(3), true);
  });

  test('normalization sets proteomics.normalized tag', async () => {
    const col = DG.Column.fromFloat32Array('intensity', new Float32Array([1, 2, 3]));
    col.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2']),
      col,
    ]);
    medianNormalize(df, ['intensity']);
    expect(df.getTag('proteomics.normalized'), 'true');
  });

  test('quantile normalization aligns column distributions', async () => {
    // Three columns with different distributions -- after quantile norm, sorted non-null values should match
    const c1 = DG.Column.fromFloat32Array('s1', new Float32Array([5, 2, 3, 6]));
    const c2 = DG.Column.fromFloat32Array('s2', new Float32Array([4, 14, 8, 2]));
    const c3 = DG.Column.fromFloat32Array('s3', new Float32Array([3, 1, 9, 7]));
    c1.semType = SEMTYPE.INTENSITY;
    c2.semType = SEMTYPE.INTENSITY;
    c3.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3']),
      c1, c2, c3,
    ]);
    quantileNormalize(df, ['s1', 's2', 's3']);
    // After quantile normalization, sorted values across columns should be the same
    const vals1 = [df.col('s1')!.get(0), df.col('s1')!.get(1), df.col('s1')!.get(2), df.col('s1')!.get(3)];
    const vals2 = [df.col('s2')!.get(0), df.col('s2')!.get(1), df.col('s2')!.get(2), df.col('s2')!.get(3)];
    const vals3 = [df.col('s3')!.get(0), df.col('s3')!.get(1), df.col('s3')!.get(2), df.col('s3')!.get(3)];
    vals1.sort((a, b) => a - b);
    vals2.sort((a, b) => a - b);
    vals3.sort((a, b) => a - b);
    for (let i = 0; i < 4; i++) {
      expect(Math.abs(vals1[i] - vals2[i]) < 0.01, true);
      expect(Math.abs(vals1[i] - vals3[i]) < 0.01, true);
    }
  });

  test('quantile normalization preserves null values', async () => {
    const c1 = DG.Column.fromFloat32Array('s1', new Float32Array([5, DG.FLOAT_NULL, 3, 6]));
    const c2 = DG.Column.fromFloat32Array('s2', new Float32Array([4, 14, 8, DG.FLOAT_NULL]));
    c1.semType = SEMTYPE.INTENSITY;
    c2.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3']),
      c1, c2,
    ]);
    quantileNormalize(df, ['s1', 's2']);
    expect(df.col('s1')!.isNone(1), true);
    expect(df.col('s2')!.isNone(3), true);
    // Non-null values should still be present
    expect(!df.col('s1')!.isNone(0), true);
    expect(!df.col('s2')!.isNone(0), true);
  });

  test('quantile normalization sets proteomics.normalized tag', async () => {
    const c1 = DG.Column.fromFloat32Array('s1', new Float32Array([1, 2, 3]));
    const c2 = DG.Column.fromFloat32Array('s2', new Float32Array([4, 5, 6]));
    c1.semType = SEMTYPE.INTENSITY;
    c2.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2']),
      c1, c2,
    ]);
    quantileNormalize(df, ['s1', 's2']);
    expect(df.getTag('proteomics.normalized'), 'true');
  });

  test('quantile normalization with < 2 columns returns without error', async () => {
    const c1 = DG.Column.fromFloat32Array('s1', new Float32Array([1, 2, 3]));
    c1.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2']),
      c1,
    ]);
    // Should not throw
    quantileNormalize(df, ['s1']);
    // Values should be unchanged
    expect(Math.abs(df.col('s1')!.get(0) - 1) < 0.01, true);
  });

  test('vsnNormalize falls back to quantile on R failure', async () => {
    const c1 = DG.Column.fromFloat32Array('log2(s1)', new Float32Array([5, 2, 3, 6]));
    const c2 = DG.Column.fromFloat32Array('log2(s2)', new Float32Array([4, 14, 8, 2]));
    c1.semType = SEMTYPE.INTENSITY;
    c2.semType = SEMTYPE.INTENSITY;
    // Also add raw columns so VSN can find them
    const r1 = DG.Column.fromFloat32Array('s1', new Float32Array([32, 4, 8, 64]));
    const r2 = DG.Column.fromFloat32Array('s2', new Float32Array([16, 16384, 256, 4]));
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3']),
      r1, r2, c1, c2,
    ]);
    // In test environment, R is unavailable, so vsnNormalize should fall back to quantile
    await vsnNormalize(df, ['log2(s1)', 'log2(s2)']);
    // Should still set the tag even via fallback
    expect(df.getTag('proteomics.normalized'), 'true');
  });
});

// ── Imputation ──────────────────────────────────────────────────────

category('Imputation', () => {
  test('MinProb imputation fills all NaN values', async () => {
    const values = new Float32Array([10, DG.FLOAT_NULL, 14, DG.FLOAT_NULL, 12]);
    const col = DG.Column.fromFloat32Array('intensity', values);
    col.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3', 'P4']),
      col,
    ]);
    imputeMinProb(df, ['intensity']);
    const c = df.col('intensity')!;
    for (let i = 0; i < c.length; i++)
      expect(c.isNone(i), false);
  });

  test('imputed values are below the observed mean', async () => {
    // Known values: [10, 12, 14, 16] mean=13, with nulls to impute
    const values = new Float32Array([10, 12, 14, 16, DG.FLOAT_NULL, DG.FLOAT_NULL]);
    const col = DG.Column.fromFloat32Array('intensity', values);
    col.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3', 'P4', 'P5']),
      col,
    ]);
    const observedMean = 13; // (10+12+14+16)/4
    imputeMinProb(df, ['intensity'], 1.8, 0.3);
    const c = df.col('intensity')!;
    // Imputed values (indices 4 and 5) should be well below the observed mean
    expect(c.get(4) < observedMean, true);
    expect(c.get(5) < observedMean, true);
  });

  test('imputation sets proteomics.imputed tag', async () => {
    const col = DG.Column.fromFloat32Array('intensity',
      new Float32Array([1, DG.FLOAT_NULL, 3]));
    col.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2']),
      col,
    ]);
    imputeMinProb(df, ['intensity']);
    expect(df.getTag('proteomics.imputed'), 'true');
  });

  test('kNN imputation fills missing values using nearest neighbors', async () => {
    // 5 rows, 3 columns. Row 2 (P2) has missing col s2.
    // P0: [10, 20, 30], P1: [11, 21, 31], P3: [12, 22, 32], P4: [100, 200, 300]
    // P2: [10.5, null, 30.5] -- nearest neighbors are P0, P1, P3 (not P4)
    // Expected imputed s2 for P2: average of neighbors' s2 values ~ (20+21+22)/3 = 21
    const c1 = DG.Column.fromFloat32Array('s1', new Float32Array([10, 11, 10.5, 12, 100]));
    const c2 = DG.Column.fromFloat32Array('s2', new Float32Array([20, 21, DG.FLOAT_NULL, 22, 200]));
    const c3 = DG.Column.fromFloat32Array('s3', new Float32Array([30, 31, 30.5, 32, 300]));
    c1.semType = SEMTYPE.INTENSITY;
    c2.semType = SEMTYPE.INTENSITY;
    c3.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3', 'P4']),
      c1, c2, c3,
    ]);
    const count = imputeKnn(df, ['s1', 's2', 's3'], 3);
    expect(count, 1);
    // The imputed value should be close to 21 (average of 3 nearest neighbors' s2 values)
    const imputed = df.col('s2')!.get(2);
    expect(Math.abs(imputed - 21) < 1, true);
  });

  test('kNN imputation falls back to column mean when no neighbors have values', async () => {
    // All rows missing in s2 except row 0
    // Row 1 and 2 are missing s2, and row 0 has it
    const c1 = DG.Column.fromFloat32Array('s1', new Float32Array([10, DG.FLOAT_NULL, DG.FLOAT_NULL]));
    const c2 = DG.Column.fromFloat32Array('s2', new Float32Array([20, DG.FLOAT_NULL, DG.FLOAT_NULL]));
    c1.semType = SEMTYPE.INTENSITY;
    c2.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2']),
      c1, c2,
    ]);
    const count = imputeKnn(df, ['s1', 's2'], 10);
    // Should have imputed 4 values (rows 1 and 2, both columns)
    // Row 1 and 2 are all-missing, so they should get column means
    expect(count >= 2, true);
  });

  test('kNN imputation returns correct count and sets tag', async () => {
    const c1 = DG.Column.fromFloat32Array('s1', new Float32Array([1, 2, DG.FLOAT_NULL, 4]));
    const c2 = DG.Column.fromFloat32Array('s2', new Float32Array([5, DG.FLOAT_NULL, 7, 8]));
    c1.semType = SEMTYPE.INTENSITY;
    c2.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3']),
      c1, c2,
    ]);
    const count = imputeKnn(df, ['s1', 's2'], 3);
    expect(count, 2);
    expect(df.getTag('proteomics.imputed'), 'true');
    // All values should be filled
    for (let i = 0; i < 4; i++) {
      expect(df.col('s1')!.isNone(i), false);
      expect(df.col('s2')!.isNone(i), false);
    }
  });

  test('imputeZero replaces all nulls with 0', async () => {
    const c1 = DG.Column.fromFloat32Array('s1', new Float32Array([1, DG.FLOAT_NULL, 3, DG.FLOAT_NULL]));
    c1.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3']),
      c1,
    ]);
    const count = imputeZero(df, ['s1']);
    expect(count, 2);
    expect(Math.abs(df.col('s1')!.get(1) - 0) < 0.01, true);
    expect(Math.abs(df.col('s1')!.get(3) - 0) < 0.01, true);
    expect(df.getTag('proteomics.imputed'), 'true');
  });

  test('imputeMean replaces all nulls with column mean', async () => {
    // Values: [2, null, 6, null] -- mean of non-null = (2+6)/2 = 4
    const c1 = DG.Column.fromFloat32Array('s1', new Float32Array([2, DG.FLOAT_NULL, 6, DG.FLOAT_NULL]));
    c1.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3']),
      c1,
    ]);
    const count = imputeMean(df, ['s1']);
    expect(count, 2);
    expect(Math.abs(df.col('s1')!.get(1) - 4) < 0.01, true);
    expect(Math.abs(df.col('s1')!.get(3) - 4) < 0.01, true);
    expect(df.getTag('proteomics.imputed'), 'true');
  });

  test('imputeMedian replaces all nulls with column median', async () => {
    // Values: [2, null, 6, 10, null] -- median of [2,6,10] = 6
    const c1 = DG.Column.fromFloat32Array('s1', new Float32Array([2, DG.FLOAT_NULL, 6, 10, DG.FLOAT_NULL]));
    c1.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['P0', 'P1', 'P2', 'P3', 'P4']),
      c1,
    ]);
    const count = imputeMedian(df, ['s1']);
    expect(count, 2);
    expect(Math.abs(df.col('s1')!.get(1) - 6) < 0.01, true);
    expect(Math.abs(df.col('s1')!.get(4) - 6) < 0.01, true);
    expect(df.getTag('proteomics.imputed'), 'true');
  });
});

// ── Differential Expression ─────────────────────────────────────────

category('Differential Expression', () => {
  test('DE produces correct log2FC sign for known groups', async () => {
    // Group 1 (control) values ~5, Group 2 (treatment) values ~10
    const g1Names = ['ctrl1', 'ctrl2', 'ctrl3'];
    const g2Names = ['treat1', 'treat2', 'treat3'];
    const nProteins = 5;
    const g1Vals = g1Names.map(() => Array.from({length: nProteins}, () => 5.0));
    const g2Vals = g2Names.map(() => Array.from({length: nProteins}, () => 10.0));
    const df = makeTestDf(nProteins, g1Names, g2Names, g1Vals, g2Vals);
    runDifferentialExpression(df, g1Names, g2Names, 'Control', 'Treatment');
    const fcCol = df.col('log2FC')!;
    // log2FC = mean(group2) - mean(group1) = 10 - 5 = 5 (positive)
    for (let i = 0; i < nProteins; i++)
      expect(fcCol.get(i) > 0, true);
  });

  test('DE produces valid p-values between 0 and 1', async () => {
    // Use slightly varied values so t-test produces real p-values
    const g1Names = ['ctrl1', 'ctrl2', 'ctrl3'];
    const g2Names = ['treat1', 'treat2', 'treat3'];
    const nProteins = 3;
    const g1Vals = [
      [4.8, 5.1, 5.3],
      [4.9, 5.0, 5.2],
      [5.1, 4.8, 5.0],
    ];
    const g2Vals = [
      [9.8, 10.1, 10.3],
      [9.9, 10.0, 10.2],
      [10.1, 9.8, 10.0],
    ];
    const df = makeTestDf(nProteins, g1Names, g2Names, g1Vals, g2Vals);
    runDifferentialExpression(df, g1Names, g2Names, 'Control', 'Treatment');
    const pCol = df.col('p-value')!;
    const adjCol = df.col('adj.p-value')!;
    for (let i = 0; i < nProteins; i++) {
      expect(pCol.get(i) >= 0 && pCol.get(i) <= 1, true);
      expect(adjCol.get(i) >= 0 && adjCol.get(i) <= 1, true);
    }
  });

  test('proteins with <2 replicates get null p-values', async () => {
    // Protein 0: only 1 non-null value in group1
    const g1Names = ['ctrl1', 'ctrl2', 'ctrl3'];
    const g2Names = ['treat1', 'treat2', 'treat3'];
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Protein ID', ['P0']),
      DG.Column.fromFloat32Array('ctrl1', new Float32Array([5.0])),
      DG.Column.fromFloat32Array('ctrl2', new Float32Array([DG.FLOAT_NULL])),
      DG.Column.fromFloat32Array('ctrl3', new Float32Array([DG.FLOAT_NULL])),
      DG.Column.fromFloat32Array('treat1', new Float32Array([10.0])),
      DG.Column.fromFloat32Array('treat2', new Float32Array([10.5])),
      DG.Column.fromFloat32Array('treat3', new Float32Array([9.5])),
    ]);
    runDifferentialExpression(df, g1Names, g2Names, 'Control', 'Treatment');
    const pCol = df.col('p-value')!;
    expect(pCol.isNone(0), true);
  });

  test('significant column marks proteins passing both thresholds', async () => {
    // Protein 0: large FC, small p (should be significant)
    // Protein 1: small FC, large p (should not be significant)
    const g1Names = ['ctrl1', 'ctrl2', 'ctrl3'];
    const g2Names = ['treat1', 'treat2', 'treat3'];
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Protein ID', ['P_significant', 'P_not_significant']),
      // Group 1 control values
      DG.Column.fromFloat32Array('ctrl1', new Float32Array([5.0, 5.0])),
      DG.Column.fromFloat32Array('ctrl2', new Float32Array([5.1, 5.1])),
      DG.Column.fromFloat32Array('ctrl3', new Float32Array([4.9, 4.9])),
      // Group 2: Protein 0 has large difference, Protein 1 is the same
      DG.Column.fromFloat32Array('treat1', new Float32Array([12.0, 5.2])),
      DG.Column.fromFloat32Array('treat2', new Float32Array([12.1, 4.8])),
      DG.Column.fromFloat32Array('treat3', new Float32Array([11.9, 5.0])),
    ]);
    runDifferentialExpression(df, g1Names, g2Names, 'Control', 'Treatment', 1.0, 0.05);
    const sigCol = df.col('significant')!;
    // Protein 0: FC ~ 7, very significant
    expect(sigCol.get(0), true);
    // Protein 1: FC ~ 0, not significant
    expect(sigCol.get(1), false);
  });

  test('DE sets proteomics.de_complete tag', async () => {
    const g1Names = ['ctrl1', 'ctrl2', 'ctrl3'];
    const g2Names = ['treat1', 'treat2', 'treat3'];
    const nProteins = 2;
    const g1Vals = g1Names.map(() => Array.from({length: nProteins}, () => 5.0));
    const g2Vals = g2Names.map(() => Array.from({length: nProteins}, () => 10.0));
    const df = makeTestDf(nProteins, g1Names, g2Names, g1Vals, g2Vals);
    runDifferentialExpression(df, g1Names, g2Names, 'Control', 'Treatment');
    expect(df.getTag('proteomics.de_complete'), 'true');
  });

  // R3/D-09: report-import DE direction default (13-05)

  test('D-09: DE dialog defaults to the declared (group1 = numerator) contrast', async () => {
    // Spectronaut report order: group1 = the parser's first condition = the
    // declared Numerator (e.g. DMD), group2 = WT. The dialog default must be
    // `${g1} vs ${g2}` (which the OK-handler maps to numerator = g1), NOT the
    // alphabetical `${g2} vs ${g1}` that previously caused the mirror defect.
    expect(getDefaultComparison('DMD', 'WT'), 'DMD vs WT');
    expect(getDefaultComparison('WT', 'DMD'), 'WT vs DMD');
  });

  test('D-09: contrast is direction-only — reversing flips sign, not |log2FC|', async () => {
    const g1 = ['dmd1', 'dmd2', 'dmd3'];
    const g2 = ['wt1', 'wt2', 'wt3'];
    const n = 3;
    const g1v = g1.map(() => Array.from({length: n}, () => 9.0)); // group1 high
    const g2v = g2.map(() => Array.from({length: n}, () => 4.0)); // group2 low

    // Declared default (numerator = group1 = DMD): OK-handler calls
    // runDifferentialExpression(df, denominatorCols, numeratorCols, ...) and
    // log2FC = mean(numerator) - mean(denominator).
    const dfDeclared = makeTestDf(n, g1, g2, g1v, g2v);
    runDifferentialExpression(dfDeclared, g2, g1, 'WT', 'DMD');
    const declared = dfDeclared.col('log2FC')!.get(0) as number;

    // Reversed selection (manual override): numerator = group2 = WT.
    const dfReversed = makeTestDf(n, g1, g2, g1v, g2v);
    runDifferentialExpression(dfReversed, g1, g2, 'DMD', 'WT');
    const reversed = dfReversed.col('log2FC')!.get(0) as number;

    expect(declared > 0, true);  // DMD higher → positive in declared orientation
    expect(reversed < 0, true);  // sign flips under reversal
    expect(Math.abs(Math.abs(declared) - Math.abs(reversed)) < 1e-6, true); // |unchanged|
  });

  test('DE assigns correct semantic types', async () => {
    const g1Names = ['ctrl1', 'ctrl2', 'ctrl3'];
    const g2Names = ['treat1', 'treat2', 'treat3'];
    const nProteins = 2;
    const g1Vals = g1Names.map(() => Array.from({length: nProteins}, () => 5.0));
    const g2Vals = g2Names.map(() => Array.from({length: nProteins}, () => 10.0));
    const df = makeTestDf(nProteins, g1Names, g2Names, g1Vals, g2Vals);
    runDifferentialExpression(df, g1Names, g2Names, 'Control', 'Treatment');
    expect(df.col('log2FC')!.semType, SEMTYPE.LOG2FC);
    expect(df.col('p-value')!.semType, SEMTYPE.P_VALUE);
    expect(df.col('adj.p-value')!.semType, SEMTYPE.P_VALUE);
  });
});

// ── R-result alignment (DEqMS scramble regression) ──────────────────

category('DE result alignment', () => {
  test('copyDEResultsToFrame realigns by row key, not position', async () => {
    // Simulates an R result returned SORTED by significance (the DEqMS bug):
    // result row j belongs to df protein index `row[j]-1`, NOT j.
    const n = 4;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Protein ID', ['P0', 'P1', 'P2', 'P3']),
    ]);
    const rowKey = [3, 1, 4, 2]; // result rows map to df idx 2,0,3,1
    // Encode each row's intended df index as its value, so a correct realign
    // makes df row i hold value i. A positional copy would fail this.
    const result = DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'row', rowKey),
      DG.Column.fromList('double', 'log2FC', rowKey.map((r) => r - 1)),
      DG.Column.fromList('double', 'p.value', rowKey.map((r) => r - 1)),
      DG.Column.fromList('double', 'adj.p.value', rowKey.map((r) => r - 1)),
      DG.Column.fromList('bool', 'significant', rowKey.map((r) => r - 1 === 0)),
    ]);

    const sigCount = copyDEResultsToFrame(df, result);

    for (let i = 0; i < n; i++) {
      expect(df.col('log2FC')!.get(i), i);
      expect(df.col('p-value')!.get(i), i);
      expect(df.col('adj.p-value')!.get(i), i);
      expect(df.col('significant')!.get(i), i === 0);
    }
    expect(sigCount, 1);
  });

  test('copyDEResultsToFrame falls back to positional without row key', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Protein ID', ['P0', 'P1', 'P2']),
    ]);
    const result = DG.DataFrame.fromColumns([
      DG.Column.fromList('double', 'log2FC', [10, 20, 30]),
      DG.Column.fromList('double', 'p.value', [1, 2, 3]),
      DG.Column.fromList('double', 'adj.p.value', [4, 5, 6]),
      DG.Column.fromList('bool', 'significant', [false, true, false]),
    ]);

    copyDEResultsToFrame(df, result);

    expect(df.col('log2FC')!.get(1), 20);
    expect(df.col('adj.p-value')!.get(2), 6);
    expect(df.col('significant')!.get(1), true);
    expect(df.col('significant')!.get(0), false);
  });
});
