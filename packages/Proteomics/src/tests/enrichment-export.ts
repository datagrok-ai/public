import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {SEMTYPE} from '../utils/proteomics-types';
import {buildEnrichmentInputsCsv} from '../analysis/enrichment-export';

/** DE'd DataFrame: G1 up, G2 down, G3 not significant (background only). */
function makeDeDf(): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Gene names', ['G1', 'G2', 'G3']),
    DG.Column.fromFloat32Array('log2FC', new Float32Array([2, -2, 0.1])),
    DG.Column.fromFloat32Array('adj.p-value', new Float32Array([0.01, 0.01, 0.5])),
  ]);
  df.col('Gene names')!.semType = SEMTYPE.GENE_SYMBOL;
  df.col('log2FC')!.semType = SEMTYPE.LOG2FC;
  df.col('adj.p-value')!.semType = SEMTYPE.P_VALUE;
  return df;
}

category('Enrichment Export', () => {
  test('buildEnrichmentInputsCsv: up/down/background split matches thresholds', async () => {
    const {csv, counts} = buildEnrichmentInputsCsv(makeDeDf(), 1.0, 0.05);
    expect(counts.up, 1, 'one up gene');
    expect(counts.down, 1, 'one down gene');
    expect(counts.background, 3, 'all three genes in background');

    const lines = csv.split('\n');
    expect(lines[0], 'gene,list', 'tidy long header');
    expect(lines.includes('G1,up'), true, 'G1 is up');
    expect(lines.includes('G2,down'), true, 'G2 is down');
    // Background is the full universe — every gene appears under background too.
    expect(lines.includes('G1,background') && lines.includes('G2,background') &&
      lines.includes('G3,background'), true, 'all genes listed as background');
    // G3 is not significant → never up or down.
    expect(lines.includes('G3,up') || lines.includes('G3,down'), false, 'G3 stays background-only');
  });

  test('buildEnrichmentInputsCsv: looser thresholds pull more genes in', async () => {
    // p ≤ 0.6 and |fc| ≥ 0 promotes G3 (fc 0.1, p 0.5) to up.
    const {counts} = buildEnrichmentInputsCsv(makeDeDf(), 0.0, 0.6);
    expect(counts.up, 2, 'G1 and G3 now up');
    expect(counts.down, 1, 'G2 still down');
    expect(counts.background, 3, 'background unchanged');
  });

  test('buildEnrichmentInputsCsv: quotes a gene containing a comma', async () => {
    const df = makeDeDf();
    df.col('Gene names')!.set(0, 'G1,alias');
    const {csv} = buildEnrichmentInputsCsv(df, 1.0, 0.05);
    expect(csv.includes('"G1,alias",up'), true, 'comma-bearing field is CSV-quoted');
  });

  test('buildEnrichmentInputsCsv: throws without DE columns', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Gene names', ['G1']),
    ]);
    df.col('Gene names')!.semType = SEMTYPE.GENE_SYMBOL;
    let threw = false;
    try { buildEnrichmentInputsCsv(df, 1.0, 0.05); } catch { threw = true; }
    expect(threw, true, 'missing log2FC/adj.p-value is an error');
  });
});
