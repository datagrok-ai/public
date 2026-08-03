import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseSpectronautCandidatesText} from '../parsers/spectronaut-candidates-parser';
import {dockComparisonFilterIfMultiContrast} from '../package';
import {SEMTYPE} from '../utils/proteomics-types';
import {getGroups} from '../analysis/experiment-setup';

/** Builds a minimal Spectronaut Candidates TSV string. */
function makeTsv(rows: string[][], headers: string[]): string {
  return [headers.join('\t'), ...rows.map((r) => r.join('\t'))].join('\n');
}

const BASE_HEADERS = [
  'Comparison (group1/group2)', 'PG.ProteinGroups', 'PG.Genes',
  'AVG Log2 Ratio', 'Qvalue', 'Pvalue',
];

function baseRow(comparison: string, protein: string, gene: string,
  log2fc: string, qval: string, pval: string): string[] {
  return [comparison, protein, gene, log2fc, qval, pval];
}

category('SpectronautCandidates', () => {
  test('parses standard rows', async () => {
    const tsv = makeTsv([
      baseRow('Treatment / Control', 'P04637', 'TP53', '2.5', '1e-5', '1e-6'),
      baseRow('Treatment / Control', 'P38398', 'BRCA1', '-2.3', '1e-4', '1e-5'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.rowCount, 2);
  });

  test('sets proteomics.source to spectronaut-candidates', async () => {
    const tsv = makeTsv([
      baseRow('A / B', 'P04637', 'TP53', '2.5', '1e-5', '1e-6'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.getTag('proteomics.source'), 'spectronaut-candidates');
  });

  test('sets de_complete and de_method tags', async () => {
    const tsv = makeTsv([
      baseRow('A / B', 'P04637', 'TP53', '2.5', '1e-5', '1e-6'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.getTag('proteomics.de_complete'), 'true');
    expect(df.getTag('proteomics.de_method'), 'spectronaut');
  });

  test('names contrast groups from the Comparison string (volcano legend)', async () => {
    const tsv = makeTsv([
      baseRow('DMT / WT', 'P04637', 'TP53', '2.5', '1e-5', '1e-6'),
      baseRow('DMT / WT', 'P38398', 'BRCA1', '-2.3', '1e-4', '1e-5'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    const groups = getGroups(df);
    expect(groups != null, true);
    // Assert the name SET (order depends on the parked sign-direction default).
    const names = [groups!.group1.name, groups!.group2.name].sort();
    expect(names[0], 'DMT');
    expect(names[1], 'WT');
    // Candidates carries no per-sample intensities → column lists stay empty.
    expect(groups!.group1.columns.length, 0);
    expect(groups!.group2.columns.length, 0);
  });

  test('renames AVG Log2 Ratio to log2FC', async () => {
    const tsv = makeTsv([
      baseRow('A / B', 'P04637', 'TP53', '2.5', '1e-5', '1e-6'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.col('log2FC') !== null, true);
    expect(df.col('AVG Log2 Ratio'), null);
    expect(Math.abs((df.col('log2FC')!.get(0) as number) - 2.5) < 0.01, true);
  });

  test('renames Qvalue to adj.p-value', async () => {
    const tsv = makeTsv([
      baseRow('A / B', 'P04637', 'TP53', '2.5', '0.001', '0.0001'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.col('adj.p-value') !== null, true);
    expect(df.col('Qvalue'), null);
    expect(Math.abs((df.col('adj.p-value')!.get(0) as number) - 0.001) < 1e-6, true);
  });

  test('renames Pvalue to p-value when present', async () => {
    const tsv = makeTsv([
      baseRow('A / B', 'P04637', 'TP53', '2.5', '0.001', '0.0001'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.col('p-value') !== null, true);
    expect(df.col('Pvalue'), null);
  });

  test('works without a Pvalue column', async () => {
    const headers = ['PG.ProteinGroups', 'PG.Genes', 'AVG Log2 Ratio', 'Qvalue'];
    const tsv = makeTsv([
      ['P04637', 'TP53', '2.5', '0.001'],
    ], headers);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.rowCount, 1);
    expect(df.col('p-value'), null);
    expect(df.col('adj.p-value') !== null, true);
  });

  test('assigns canonical semantic types', async () => {
    const tsv = makeTsv([
      baseRow('A / B', 'P04637', 'TP53', '2.5', '0.001', '0.0001'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.col('PG.ProteinGroups')!.semType, SEMTYPE.PROTEIN_ID);
    expect(df.col('PG.Genes')!.semType, SEMTYPE.GENE_SYMBOL);
    expect(df.col('log2FC')!.semType, SEMTYPE.LOG2FC);
    expect(df.col('adj.p-value')!.semType, SEMTYPE.P_VALUE);
    expect(df.col('p-value')!.semType, SEMTYPE.P_VALUE);
  });

  test('computes significant column from |log2FC| >= 1 and adj.p <= 0.05', async () => {
    const tsv = makeTsv([
      baseRow('A / B', 'P_sig_up',     'A', '2.5',  '0.001', '0.0001'),
      baseRow('A / B', 'P_sig_down',   'B', '-1.5', '0.001', '0.0001'),
      baseRow('A / B', 'P_small_fc',   'C', '0.3',  '0.001', '0.0001'),
      baseRow('A / B', 'P_big_qvalue', 'D', '2.0',  '0.5',   '0.4'),
      baseRow('A / B', 'P_just_inside', 'E', '1.01', '0.04', '0.03'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    const sig = df.col('significant')!;
    expect(sig.get(0), true);   // big FC, tiny q
    expect(sig.get(1), true);   // big negative FC, tiny q
    expect(sig.get(2), false);  // FC below threshold
    expect(sig.get(3), false);  // q above threshold
    expect(sig.get(4), true);   // just inside both thresholds
  });

  test('filters contam_ rows', async () => {
    const tsv = makeTsv([
      baseRow('A / B', 'P04637', 'TP53', '2.5', '0.001', '0.0001'),
      baseRow('A / B', 'contam_P99999', 'KRT1', '0.1', '0.5', '0.4'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.rowCount, 1);
    expect(df.col('PG.ProteinGroups')!.get(0), 'P04637');
  });

  test('filters rev_ decoy rows', async () => {
    const tsv = makeTsv([
      baseRow('A / B', 'P04637', 'TP53', '2.5', '0.001', '0.0001'),
      baseRow('A / B', 'rev_P12345', 'FAKE', '0.1', '0.5', '0.4'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.rowCount, 1);
  });

  test('throws on missing protein-group column', async () => {
    const headers = ['Comparison', 'AVG Log2 Ratio', 'Qvalue'];
    const tsv = makeTsv([['A / B', '2.5', '0.001']], headers);
    let threw = false;
    try {
      await parseSpectronautCandidatesText(tsv);
    } catch (e: any) {
      threw = true;
      expect(String(e.message).includes('protein-group'), true);
    }
    expect(threw, true);
  });

  test('throws on missing log2 ratio column', async () => {
    const headers = ['PG.ProteinGroups', 'PG.Genes', 'Qvalue'];
    const tsv = makeTsv([['P04637', 'TP53', '0.001']], headers);
    let threw = false;
    try {
      await parseSpectronautCandidatesText(tsv);
    } catch (e: any) {
      threw = true;
      expect(String(e.message).includes('log2 ratio'), true);
    }
    expect(threw, true);
  });

  test('throws on missing q-value column', async () => {
    const headers = ['PG.ProteinGroups', 'PG.Genes', 'AVG Log2 Ratio'];
    const tsv = makeTsv([['P04637', 'TP53', '2.5']], headers);
    let threw = false;
    try {
      await parseSpectronautCandidatesText(tsv);
    } catch (e: any) {
      threw = true;
      expect(String(e.message).includes('q-value'), true);
    }
    expect(threw, true);
  });

  test('accepts alternative column names (ProteinGroups, Log2 Ratio, Q-Value)', async () => {
    const headers = ['ProteinGroups', 'Genes', 'Log2 Ratio', 'Q-Value'];
    const tsv = makeTsv([['P04637', 'TP53', '2.5', '0.001']], headers);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.rowCount, 1);
    expect(df.col('log2FC') !== null, true);
    expect(df.col('adj.p-value') !== null, true);
  });

  test('keeps Comparison column for multi-comparison files', async () => {
    const tsv = makeTsv([
      baseRow('Treatment1 / Control', 'P04637', 'TP53',  '2.5', '0.001', '0.0001'),
      baseRow('Treatment2 / Control', 'P04637', 'TP53',  '1.8', '0.002', '0.0005'),
      baseRow('Treatment1 / Control', 'P38398', 'BRCA1', '-2.3', '0.001', '0.0001'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    expect(df.rowCount, 3);
    expect(df.col('Comparison (group1/group2)') !== null, true);
  });

  // --- R3/D-08: per-row sign normalization (13-03) ---

  const QTY_HEADERS = [
    'Comparison (group1/group2)', 'PG.ProteinGroups', 'PG.Genes',
    'AVG Log2 Ratio', 'Qvalue', 'Pvalue',
    'AVG Group Quantity Numerator', 'AVG Group Quantity Denominator',
    'Condition Numerator', 'Condition Denominator',
  ];

  test('flips reversed-comparison rows; canonical rows untouched', async () => {
    // First declared comparison ("DMD / WT") is canonical. The "WT / DMD" row
    // is its exact reverse → log2FC negated, AVG Group Quantity swapped,
    // Comparison + Condition relabeled. The canonical row is byte-identical.
    const tsv = makeTsv([
      ['DMD / WT', 'P1', 'G1', '-2.0', '0.001', '0.0001', '100', '400', 'DMD', 'WT'],
      ['WT / DMD', 'P2', 'G2', '3.0', '0.001', '0.0001', '800', '100', 'WT', 'DMD'],
    ], QTY_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    const cmp = df.col('Comparison (group1/group2)')!;
    const fc = df.col('log2FC')!;
    const num = df.col('AVG Group Quantity Numerator')!;
    const den = df.col('AVG Group Quantity Denominator')!;
    // Canonical row 0 — unchanged.
    expect(cmp.get(0), 'DMD / WT');
    expect(Math.abs((fc.get(0) as number) - (-2.0)) < 1e-6, true);
    expect(Math.abs((num.get(0) as number) - 100) < 1e-4, true);
    expect(Math.abs((den.get(0) as number) - 400) < 1e-4, true);
    // Reversed row 1 — flipped.
    expect(cmp.get(1), 'DMD / WT');
    expect(Math.abs((fc.get(1) as number) - (-3.0)) < 1e-6, true);
    expect(Math.abs((num.get(1) as number) - 100) < 1e-4, true); // swapped 800<->100
    expect(Math.abs((den.get(1) as number) - 800) < 1e-4, true);
    expect(df.col('Condition Numerator')!.get(1), 'DMD');
    expect(df.col('Condition Denominator')!.get(1), 'WT');
  });

  test('flips log2FC without AVG Group Quantity columns (no swap, no throw)', async () => {
    const tsv = makeTsv([
      baseRow('DMD / WT', 'P1', 'G1', '-2.0', '0.001', '0.0001'),
      baseRow('WT / DMD', 'P2', 'G2', '3.0', '0.001', '0.0001'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    const fc = df.col('log2FC')!;
    expect(Math.abs((fc.get(0) as number) - (-2.0)) < 1e-6, true); // canonical untouched
    expect(Math.abs((fc.get(1) as number) - (-3.0)) < 1e-6, true); // reversed flipped
    expect(df.col('Comparison (group1/group2)')!.get(1), 'DMD / WT');
    expect(df.col('AVG Group Quantity Numerator'), null); // never created
  });

  test('single-orientation file is not mirrored', async () => {
    const tsv = makeTsv([
      baseRow('Treatment / Control', 'P1', 'A', '2.5', '0.001', '0.0001'),
      baseRow('Treatment / Control', 'P2', 'B', '-1.7', '0.002', '0.0005'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    const fc = df.col('log2FC')!;
    expect(Math.abs((fc.get(0) as number) - 2.5) < 1e-6, true);
    expect(Math.abs((fc.get(1) as number) - (-1.7)) < 1e-6, true);
    expect(df.col('Comparison (group1/group2)')!.get(0), 'Treatment / Control');
  });

  test('parses CSV when delimiter is comma', async () => {
    const csv = 'PG.ProteinGroups,PG.Genes,AVG Log2 Ratio,Qvalue\n' +
      'P04637,TP53,2.5,0.001\n' +
      'P38398,BRCA1,-2.3,0.001';
    const df = await parseSpectronautCandidatesText(csv);
    expect(df.rowCount, 2);
    expect(df.col('log2FC') !== null, true);
  });

  test('creates Primary Protein ID for semicolon-delimited IDs', async () => {
    const tsv = makeTsv([
      baseRow('A / B', 'P04637;Q9Y6R7', 'TP53', '2.5', '0.001', '0.0001'),
    ], BASE_HEADERS);
    const df = await parseSpectronautCandidatesText(tsv);
    const primary = df.col('Primary Protein ID');
    expect(primary !== null, true);
    expect(primary!.get(0), 'P04637');
  });
});

/**
 * 14-03 G4 + D-05: the unified Filters viewer must scope to Comparison +
 * Display Name + Source ID, never the boolean Flags column that Phase 13
 * observed the legacy columnNames allowlist auto-attaching.
 */
function findFiltersViewer(tv: DG.TableView): DG.Viewer | null {
  for (const v of tv.viewers) {
    if (v.type === DG.VIEWER.FILTERS) return v;
  }
  return null;
}

function makeMultiContrastDf(opts: {
  withDisplayName?: boolean;
  withSourceId?: boolean;
  withFlags?: boolean;
  withProteinId?: boolean;
  comparisons?: string[];
} = {}): DG.DataFrame {
  const {
    withDisplayName = true, withSourceId = true,
    withFlags = true, withProteinId = true,
    comparisons = ['T1 / Control', 'T2 / Control', 'T1 / Control', 'T2 / Control'],
  } = opts;

  const cols: DG.Column[] = [
    DG.Column.fromStrings('Comparison (group1/group2)', comparisons),
    DG.Column.fromList('double' as DG.ColumnType, 'log2FC',
      Array.from({length: comparisons.length}, (_, i) => i - 1.5)),
  ];
  if (withProteinId) {
    const pid = DG.Column.fromStrings('Primary Protein ID',
      comparisons.map((_, i) => `P${10000 + i}`));
    pid.semType = SEMTYPE.PROTEIN_ID;
    cols.push(pid);
  }
  if (withDisplayName) {
    const dn = DG.Column.fromStrings('Display Name',
      comparisons.map((_, i) => `Gene${i + 1}`));
    dn.semType = SEMTYPE.DISPLAY_NAME;
    cols.push(dn);
  }
  if (withSourceId) {
    const sid = DG.Column.fromStrings('Source ID',
      comparisons.map((_, i) => `ENSG${10000 + i}`));
    sid.semType = SEMTYPE.SOURCE_ID;
    cols.push(sid);
  }
  if (withFlags) {
    cols.push(DG.Column.fromBitSet('Flags',
      DG.BitSet.create(comparisons.length, (i) => i % 2 === 0)));
  }
  return DG.DataFrame.fromColumns(cols);
}

category('Proteomics: 14-03', () => {
  test('filtersScopingNoFlags', async () => {
    // The dock contract is observable; the column-membership of the docked
    // viewer is NOT observable through getOptions().look — commit e527d07ba1
    // discovered the platform serializer strips both `filters[]` and
    // `columnNames` subtrees regardless of which shape we feed in. Round-3
    // human UAT (2026-06-04) confirmed the viewer renders the requested
    // columns at runtime. This test pins the docking decision; the
    // visual-membership and Flags-exclusion contracts are covered by the
    // Phase 13 round-3 HUMAN-UAT record.
    const df = makeMultiContrastDf();
    const tv = grok.shell.addTableView(df);
    try {
      const docked = dockComparisonFilterIfMultiContrast(tv, df);
      expect(docked, true);
      expect(findFiltersViewer(tv) !== null, true);
    } finally {
      tv.close();
    }
  });

  test('filtersScopingSingleContrast', async () => {
    const df = makeMultiContrastDf({comparisons: ['T1 / Control', 'T1 / Control', 'T1 / Control']});
    const tv = grok.shell.addTableView(df);
    try {
      const docked = dockComparisonFilterIfMultiContrast(tv, df);
      expect(docked, false);
      expect(findFiltersViewer(tv), null);
    } finally {
      tv.close();
    }
  });

  test('filtersScopingFallbackToProteinIdWhenNoDisplayName', async () => {
    // The DataFrame lacks Display Name + Source ID; assert the docking
    // decision still fires (Protein ID fallback path). Per the note above,
    // column-membership is not observable via getOptions().
    const df = makeMultiContrastDf({withDisplayName: false, withSourceId: false});
    const tv = grok.shell.addTableView(df);
    try {
      const docked = dockComparisonFilterIfMultiContrast(tv, df);
      expect(docked, true);
      expect(findFiltersViewer(tv) !== null, true);
    } finally {
      tv.close();
    }
  });
});
