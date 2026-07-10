import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {createAbundanceCorrelation, correlationDescription, pearson, spearman} from '../viewers/abundance-correlation';
import {setGroups} from '../analysis/experiment-setup';
import {SEMTYPE} from '../utils/proteomics-types';

/** Candidates-shaped mock: the pre-computed AVG Group Quantity columns + a
 * `significant` column. createAbundanceCorrelation resolves abundance from the
 * quantity columns (no per-sample data needed). */
function makeCandidatesDf(num: number[], den: number[]): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('ProteinGroups', num.map((_, i) => `P${i + 1}`)),
    DG.Column.fromFloat32Array('AVG Group Quantity Numerator', new Float32Array(num)),
    DG.Column.fromFloat32Array('AVG Group Quantity Denominator', new Float32Array(den)),
  ]);
  df.col('ProteinGroups')!.semType = SEMTYPE.PROTEIN_ID;
  setGroups(df, {group1: {name: 'DMD', columns: []}, group2: {name: 'WT', columns: []}});
  return df;
}

/** Report-shaped mock: per-sample intensity columns assigned to two groups, plus
 * a `direction` column so the color binding has a target. Values are log2-scale
 * (all < 45), so the geometric-mean path rebases them as logs. */
function makeReportDf(): DG.DataFrame {
  const cols: DG.Column[] = [
    DG.Column.fromStrings('Primary Protein ID', ['P1', 'P2', 'P3', 'P4', 'P5']),
    DG.Column.fromList('double' as DG.ColumnType, 'S1', [1, 2, 3, 4, 5]),
    DG.Column.fromList('double' as DG.ColumnType, 'S2', [3, 4, 5, 6, 7]),
    DG.Column.fromList('double' as DG.ColumnType, 'S3', [5, 6, 7, 8, 9]),
    DG.Column.fromList('double' as DG.ColumnType, 'S4', [7, 8, 9, 10, 11]),
    DG.Column.fromStrings('direction',
      ['Enriched in Control', 'Enriched in Treatment', 'Not significant',
        'Enriched in Control', 'Enriched in Treatment']),
  ];
  cols[0].semType = SEMTYPE.PROTEIN_ID;
  const df = DG.DataFrame.fromColumns(cols);
  setGroups(df, {
    group1: {name: 'Control', columns: ['S1', 'S2']},
    group2: {name: 'Treatment', columns: ['S3', 'S4']},
  });
  return df;
}

category('Abundance-Correlation', () => {
  test('pearson is 1 for a perfectly correlated series', async () => {
    const {r, n} = pearson([1, 2, 3, 4, 5], [2, 4, 6, 8, 10]);
    expect(n, 5);
    expect(Math.abs(r - 1) < 1e-6, true);
  });

  test('pearson is -1 for a perfectly anti-correlated series', async () => {
    const {r} = pearson([1, 2, 3, 4], [4, 3, 2, 1]);
    expect(Math.abs(r + 1) < 1e-6, true);
  });

  test('pearson ignores rows with a null in either series', async () => {
    const {n} = pearson([1, null, 3, 4], [2, 2, null, 8]);
    expect(n, 2); // only rows 0 and 3 are pairwise-complete
  });

  test('spearman is 1 with rank ties', async () => {
    const rho = spearman([1, 2, 2, 3, 4], [1, 2, 2, 3, 4]);
    expect(Math.abs(rho - 1) < 1e-9, true);
  });

  test('builds a scatter on the log10 abundance columns (Candidates)', async () => {
    const df = makeCandidatesDf([100, 1000, 10, 50], [50, 500, 5, 25]);
    df.name = 'Abundance correlation candidates test';
    const tv = grok.shell.addTableView(df);
    try {
      const sp = createAbundanceCorrelation(df);
      expect(sp !== null, true);
      // hidden log10 columns added for both conditions
      expect(df.col('log10 abundance: DMD') !== null, true);
      expect(df.col('log10 abundance: WT') !== null, true);
      // axes: x = WT (denominator / group2), y = DMD (numerator / group1)
      expect(sp!.props.xColumnName, 'log10 abundance: WT');
      expect(sp!.props.yColumnName, 'log10 abundance: DMD');
      // Pearson/Spearman annotation goes in the on-canvas description overlay
      // (the title is not painted in a docked viewer).
      const desc = correlationDescription('DMD', 'WT', df, 'log10 abundance: WT', 'log10 abundance: DMD');
      expect(desc.includes('Pearson'), true);
      expect(desc.includes('Spearman'), true);
      expect(desc.includes('DMD vs WT'), true);
    } finally {
      tv.close();
    }
  });

  test('Report path uses the geometric mean and colors by direction', async () => {
    const df = makeReportDf();
    df.name = 'Abundance correlation report test';
    const tv = grok.shell.addTableView(df);
    try {
      const sp = createAbundanceCorrelation(df);
      expect(sp !== null, true);
      expect(sp!.props.colorColumnName, 'direction');
      // Geometric mean of a log2-scale group = log10(2) * mean(log2 values).
      // Row 0 Control {S1=1, S2=3} → mean=2 → 2*log10(2) ≈ 0.60206.
      const yCol = df.col('log10 abundance: Control')!;
      expect(Math.abs((yCol.get(0) as number) - 2 * Math.log10(2)) < 1e-5, true);
      // Row 0 Treatment {S3=5, S4=7} → mean=6 → 6*log10(2) ≈ 1.80618.
      const xCol = df.col('log10 abundance: Treatment')!;
      expect(Math.abs((xCol.get(0) as number) - 6 * Math.log10(2)) < 1e-5, true);
    } finally {
      tv.close();
    }
  });

  test('draws a y=x diagonal reference line', async () => {
    const df = makeReportDf();
    const tv = grok.shell.addTableView(df);
    try {
      createAbundanceCorrelation(df);
      const lines = df.meta.formulaLines.items;
      const diagonal = lines.find((l: any) => {
        const f = (l.formula ?? '') as string;
        return typeof f === 'string' &&
          f.includes('log10 abundance: Control') && f.includes('log10 abundance: Treatment');
      });
      expect(diagonal !== undefined, true);
    } finally {
      tv.close();
    }
  });

  test('returns null when no abundance is present', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('ProteinGroups', ['P1', 'P2']),
      DG.Column.fromFloat32Array('log2FC', new Float32Array([1, -1])),
    ]);
    expect(createAbundanceCorrelation(df), null);
  });
});
