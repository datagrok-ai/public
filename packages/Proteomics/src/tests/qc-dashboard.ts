import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {
  computeMA,
  computeCV,
  computeLoessTrend,
  createMissingnessMatrix,
  unpivotIntensities,
  computeMissingBarData,
  getIntensityColumns,
} from '../viewers/qc-computations';
import {GroupAssignment} from '../analysis/experiment-setup';
import {SEMTYPE} from '../utils/proteomics-types';

category('QC Dashboard', () => {
  test('getIntensityColumns', async () => {
    const col1 = DG.Column.fromFloat32Array('log2(sample1)', new Float32Array([1, 2]));
    col1.semType = SEMTYPE.INTENSITY;
    const col2 = DG.Column.fromFloat32Array('other', new Float32Array([3, 4]));
    const col3 = DG.Column.fromFloat32Array('log2(sample2)', new Float32Array([5, 6]));
    col3.semType = SEMTYPE.INTENSITY;
    const df = DG.DataFrame.fromColumns([col1, col2, col3]);

    const result = getIntensityColumns(df);
    expect(result.length, 2);
    expect(result[0], 'log2(sample1)');
    expect(result[1], 'log2(sample2)');
  });

  test('MA computation', async () => {
    // 3 rows, 4 intensity columns: 2 per group
    const c1 = DG.Column.fromFloat32Array('log2(ctrl1)', new Float32Array([10, 5, 8]));
    c1.semType = SEMTYPE.INTENSITY;
    const c2 = DG.Column.fromFloat32Array('log2(ctrl2)', new Float32Array([12, 7, 6]));
    c2.semType = SEMTYPE.INTENSITY;
    const c3 = DG.Column.fromFloat32Array('log2(treat1)', new Float32Array([8, 3, DG.FLOAT_NULL]));
    c3.semType = SEMTYPE.INTENSITY;
    const c4 = DG.Column.fromFloat32Array('log2(treat2)', new Float32Array([6, 5, DG.FLOAT_NULL]));
    c4.semType = SEMTYPE.INTENSITY;

    const df = DG.DataFrame.fromColumns([c1, c2, c3, c4]);
    const groups: GroupAssignment = {
      group1: {name: 'Control', columns: ['log2(ctrl1)', 'log2(ctrl2)']},
      group2: {name: 'Treatment', columns: ['log2(treat1)', 'log2(treat2)']},
    };

    computeMA(df, groups);

    const mCol = df.col('M')!;
    const aCol = df.col('A')!;

    // Row 0: g1Mean = (10+12)/2 = 11, g2Mean = (8+6)/2 = 7
    // M = g2 - g1 = 7 - 11 = -4, A = (11+7)/2 = 9
    expect(Math.abs(mCol.get(0) - (-4)) < 0.01, true);
    expect(Math.abs(aCol.get(0) - 9) < 0.01, true);

    // Row 2: all null in group2 -> M and A should be null
    expect(mCol.isNone(2), true);
    expect(aCol.isNone(2), true);
  });

  test('CV computation', async () => {
    // 2 rows, 3 columns for one group
    // Row 0: all identical values (2, 2, 2) -> CV = 0
    // Row 1: varied values -> CV > 0
    const c1 = DG.Column.fromFloat32Array('log2(s1)', new Float32Array([2, 3]));
    c1.semType = SEMTYPE.INTENSITY;
    const c2 = DG.Column.fromFloat32Array('log2(s2)', new Float32Array([2, 5]));
    c2.semType = SEMTYPE.INTENSITY;
    const c3 = DG.Column.fromFloat32Array('log2(s3)', new Float32Array([2, 7]));
    c3.semType = SEMTYPE.INTENSITY;

    const df = DG.DataFrame.fromColumns([c1, c2, c3]);
    computeCV(df, ['log2(s1)', 'log2(s2)', 'log2(s3)'], 'CV', 'Mean');

    const cvCol = df.col('CV')!;
    const meanCol = df.col('Mean')!;

    // Row 0: all log2 values = 2, raw = 2^2 = 4 for all. sd = 0, CV = 0. Mean = 4.
    expect(Math.abs(cvCol.get(0)) < 0.001, true);
    expect(Math.abs(meanCol.get(0) - 4) < 0.01, true);

    // Row 1: varied values -> CV > 0
    expect(cvCol.get(1) > 0, true);
  });

  test('Loess/moving-average trend', async () => {
    // Create a DataFrame with M and A columns, 20+ rows
    const n = 25;
    // Use a seeded LCG instead of Math.random so the smoothness assertion
    // below isn't flaky on adversarial RNG draws.
    let seed = 12345;
    const rand = () => {
      seed = (seed * 1664525 + 1013904223) >>> 0;
      return seed / 0xFFFFFFFF;
    };
    const aValues = new Float32Array(n);
    const mValues = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      aValues[i] = i + 1; // A from 1 to 25
      mValues[i] = 0.5 * Math.sin(i / 5) + (rand() - 0.5) * 0.3; // pattern + seeded noise
    }

    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array('A', aValues),
      DG.Column.fromFloat32Array('M', mValues),
    ]);

    computeLoessTrend(df, 0.5);

    // Assert MA_trend column exists
    const trendCol = df.col('MA_trend');
    expect(trendCol !== null, true);

    // Assert not all null
    let nonNullCount = 0;
    for (let i = 0; i < n; i++) {
      if (!trendCol!.isNone(i))
        nonNullCount++;
    }
    expect(nonNullCount > 0, true);

    // Assert trend is smoother than raw M: variance(trend) < variance(M)
    let sumM = 0; let sumM2 = 0; let countM = 0;
    let sumT = 0; let sumT2 = 0; let countT = 0;
    for (let i = 0; i < n; i++) {
      const m = mValues[i];
      sumM += m; sumM2 += m * m; countM++;
      if (!trendCol!.isNone(i)) {
        const t = trendCol!.get(i) as number;
        sumT += t; sumT2 += t * t; countT++;
      }
    }
    const varM = (sumM2 / countM) - (sumM / countM) * (sumM / countM);
    const varT = (sumT2 / countT) - (sumT / countT) * (sumT / countT);
    expect(varT < varM, true);
  });

  test('Missingness matrix', async () => {
    const c1 = DG.Column.fromFloat32Array('log2(s1)', new Float32Array([5.0, DG.FLOAT_NULL]));
    c1.semType = SEMTYPE.INTENSITY;
    const c2 = DG.Column.fromFloat32Array('log2(s2)', new Float32Array([DG.FLOAT_NULL, 3.0]));
    c2.semType = SEMTYPE.INTENSITY;

    const df = DG.DataFrame.fromColumns([c1, c2]);
    const result = createMissingnessMatrix(df, ['log2(s1)', 'log2(s2)']);

    expect(result.rowCount, 2);
    expect(result.name, 'Missing Values');

    // Row 0: s1=present(1), s2=missing(0)
    const rc1 = result.col('log2(s1)')!;
    const rc2 = result.col('log2(s2)')!;
    expect(rc1.get(0), 1);
    expect(rc2.get(0), 0);
    // Row 1: s1=missing(0), s2=present(1)
    expect(rc1.get(1), 0);
    expect(rc2.get(1), 1);
  });

  test('Unpivot intensities', async () => {
    const idCol = DG.Column.fromStrings('Protein ID', ['P0', 'P1']);
    idCol.semType = SEMTYPE.PROTEIN_ID;
    const c1 = DG.Column.fromFloat32Array('log2(s1)', new Float32Array([10, DG.FLOAT_NULL]));
    c1.semType = SEMTYPE.INTENSITY;
    const c2 = DG.Column.fromFloat32Array('log2(s2)', new Float32Array([20, 30]));
    c2.semType = SEMTYPE.INTENSITY;

    const df = DG.DataFrame.fromColumns([idCol, c1, c2]);
    const result = unpivotIntensities(df, ['log2(s1)', 'log2(s2)']);

    expect(result.name, 'Intensity Distributions');
    // 3 non-null values: P0/s1(10), P0/s2(20), P1/s2(30)
    expect(result.rowCount, 3);
    expect(result.col('ProteinId') !== null, true);
    expect(result.col('Sample') !== null, true);
    expect(result.col('Intensity') !== null, true);
  });

  test('Missing bar data', async () => {
    // 3 rows, 2 intensity columns (1 per group)
    // Col 0: 1 null out of 3 (33.3%)
    // Col 1: 0 null out of 3 (0%)
    const c1 = DG.Column.fromFloat32Array('log2(ctrl1)', new Float32Array([5, DG.FLOAT_NULL, 7]));
    c1.semType = SEMTYPE.INTENSITY;
    const c2 = DG.Column.fromFloat32Array('log2(treat1)', new Float32Array([10, 11, 12]));
    c2.semType = SEMTYPE.INTENSITY;

    const df = DG.DataFrame.fromColumns([c1, c2]);
    const groups: GroupAssignment = {
      group1: {name: 'Control', columns: ['log2(ctrl1)']},
      group2: {name: 'Treatment', columns: ['log2(treat1)']},
    };

    const result = computeMissingBarData(df, ['log2(ctrl1)', 'log2(treat1)'], groups);

    expect(result.rowCount, 2);
    const pctCol = result.col('MissingPct')!;
    const grpCol = result.col('Group')!;

    // Row 0 (ctrl1): ~33.3% missing
    expect(Math.abs(pctCol.get(0) - 33.333) < 1, true);
    // Row 1 (treat1): 0% missing
    expect(Math.abs(pctCol.get(1)) < 0.01, true);

    // Group assignments
    expect(grpCol.get(0), 'Control');
    expect(grpCol.get(1), 'Treatment');
  });
});
