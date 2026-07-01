import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {
  createGroupMeanCorrelation,
  computeGroupMeans,
  pearson,
  spearman,
} from '../viewers/group-mean-correlation';
import {setGroups} from '../analysis/experiment-setup';
import {SEMTYPE} from '../utils/proteomics-types';

/** Builds a 5-row, 4-sample fixture DataFrame with a `direction` column so the
 *  correlation viewer's color binding has something to point at. */
function makeFixtureDf(): DG.DataFrame {
  const cols: DG.Column[] = [];
  const pid = DG.Column.fromStrings('Primary Protein ID',
    ['P1', 'P2', 'P3', 'P4', 'P5']);
  pid.semType = SEMTYPE.PROTEIN_ID;
  cols.push(pid);

  const display = DG.Column.fromStrings('Display Name',
    ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5']);
  display.semType = SEMTYPE.DISPLAY_NAME;
  cols.push(display);

  // Intensities — known per-row mean of [g1,g2] is monotonic for groupMeanColumnsCreated.
  // Row 0: g1={1,3}→mean=2, g2={5,7}→mean=6
  // Row 1: g1={2,4}→mean=3, g2={6,8}→mean=7
  // Row 2: g1={3,5}→mean=4, g2={7,9}→mean=8
  // Row 3: g1={4,6}→mean=5, g2={8,10}→mean=9
  // Row 4: g1={5,7}→mean=6, g2={9,11}→mean=10
  cols.push(DG.Column.fromList('double' as DG.ColumnType, 'S1', [1, 2, 3, 4, 5]));
  cols.push(DG.Column.fromList('double' as DG.ColumnType, 'S2', [3, 4, 5, 6, 7]));
  cols.push(DG.Column.fromList('double' as DG.ColumnType, 'S3', [5, 6, 7, 8, 9]));
  cols.push(DG.Column.fromList('double' as DG.ColumnType, 'S4', [7, 8, 9, 10, 11]));

  // direction column — the viewer reads this for color binding (set by Plan 14-02 ensureDirectionColumn).
  cols.push(DG.Column.fromStrings('direction',
    ['Enriched in Control', 'Enriched in Treatment', 'Not significant',
      'Enriched in Control', 'Enriched in Treatment']));

  const df = DG.DataFrame.fromColumns(cols);
  setGroups(df, {
    group1: {name: 'Control', columns: ['S1', 'S2']},
    group2: {name: 'Treatment', columns: ['S3', 'S4']},
  });
  df.name = `correlation-fixture-${Math.random().toString(36).slice(2, 8)}`;
  return df;
}

category('Proteomics: 14-04', () => {
  test('groupMeanColumnsCreated', async () => {
    const df = makeFixtureDf();
    computeGroupMeans(df, {
      group1: {name: 'Control', columns: ['S1', 'S2']},
      group2: {name: 'Treatment', columns: ['S3', 'S4']},
    });
    const num = df.col('Numerator Mean');
    const den = df.col('Denominator Mean');
    expect(num !== null, true);
    expect(den !== null, true);
    expect(num!.semType, SEMTYPE.NUMERATOR_MEAN);
    expect(den!.semType, SEMTYPE.DENOMINATOR_MEAN);
    expect(num!.length, df.rowCount);
    expect(den!.length, df.rowCount);
    // Row 0 means: g1=(1+3)/2=2, g2=(5+7)/2=6.
    expect(Math.abs((num!.get(0) as number) - 2) < 1e-6, true);
    expect(Math.abs((den!.get(0) as number) - 6) < 1e-6, true);
    // Row 4 means: g1=(5+7)/2=6, g2=(9+11)/2=10.
    expect(Math.abs((num!.get(4) as number) - 6) < 1e-6, true);
    expect(Math.abs((den!.get(4) as number) - 10) < 1e-6, true);
  });

  test('groupMeanColumnsReRunSafe', async () => {
    const df = makeFixtureDf();
    const groups = {
      group1: {name: 'Control', columns: ['S1', 'S2']},
      group2: {name: 'Treatment', columns: ['S3', 'S4']},
    };
    computeGroupMeans(df, groups);
    computeGroupMeans(df, groups);
    // No duplicates after a second invocation.
    const numCount = df.columns.toList().filter((c) => c.name === 'Numerator Mean').length;
    const denCount = df.columns.toList().filter((c) => c.name === 'Denominator Mean').length;
    expect(numCount, 1);
    expect(denCount, 1);
    // Values still correct.
    const num = df.col('Numerator Mean')!;
    expect(Math.abs((num.get(0) as number) - 2) < 1e-6, true);
  });

  test('pearsonPerfectCorrelation', async () => {
    const r = pearson([1, 2, 3, 4, 5], [2, 4, 6, 8, 10]);
    expect(Math.abs(r - 1) < 1e-9, true);
  });

  test('pearsonNoCorrelation', async () => {
    const r = pearson([1, 2, 3, 4, 5], [3, 1, 4, 1, 5]);
    // Not exactly 0 — just within a reasonable tolerance away from 1.
    expect(Math.abs(r) < 0.6, true);
  });

  test('spearmanRankTies', async () => {
    const rho = spearman([1, 2, 2, 3, 4], [1, 2, 2, 3, 4]);
    expect(Math.abs(rho - 1) < 1e-9, true);
  });

  test('createGroupMeanCorrelationFactory', async () => {
    const df = makeFixtureDf();
    const sp = createGroupMeanCorrelation(df);
    expect(sp !== null, true);
    const opts: any = (sp as any).getOptions?.();
    const title: string = opts?.look?.title ?? '';
    expect(title.includes('r='), true);
    expect(title.includes('ρ='), true);
    expect(sp.props.xColumnName, 'Numerator Mean');
    expect(sp.props.yColumnName, 'Denominator Mean');
    expect(sp.props.colorColumnName, 'direction');
  });

  test('correlationDiagonalLine', async () => {
    const df = makeFixtureDf();
    createGroupMeanCorrelation(df);
    const lines = df.meta.formulaLines.items;
    const diagonal = lines.find((l: any) => {
      const f = (l.formula ?? '') as string;
      return typeof f === 'string' &&
        f.includes('Numerator Mean') && f.includes('Denominator Mean');
    });
    expect(diagonal !== undefined, true);
    expect((diagonal as any).color, '#888888');
  });

  test('correlationLabelBindsToDisplayName', async () => {
    const df = makeFixtureDf();
    const sp = createGroupMeanCorrelation(df);
    expect((sp.props.labelColumnNames as string[])?.includes('Display Name'), true);
  });
});
