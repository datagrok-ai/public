import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {
  CHART_TOP_MARKER_COL,
  NEG_LOG10_FDR_COL,
  openEnrichmentVisualization,
  prepareEnrichmentChartColumns,
  wireEnrichmentToVolcano,
} from '../viewers/enrichment-viewers';
import {SEMTYPE} from '../utils/proteomics-types';

// --- Helpers ---

/** Creates a mock enrichment DataFrame with the specified number of rows. */
function makeMockEnrichmentDf(rowCount: number, longTermName?: string): DG.DataFrame {
  const sources: string[] = [];
  const termIds: string[] = [];
  const termNames: string[] = [];
  const pValues: number[] = [];
  const fdrs: number[] = [];
  const geneCounts: number[] = [];
  const geneRatios: number[] = [];
  const intersections: string[] = [];
  const significants: boolean[] = [];

  for (let i = 0; i < rowCount; i++) {
    sources.push(i % 2 === 0 ? 'GO:BP' : 'KEGG');
    termIds.push(`GO:${String(i).padStart(7, '0')}`);
    termNames.push(longTermName && i === 0 ? longTermName : `Term ${i + 1}`);
    const fdr = 0.001 * (i + 1);
    pValues.push(fdr);
    fdrs.push(fdr);
    geneCounts.push(10 - Math.min(i, 9));
    geneRatios.push(0.2 - i * 0.005);
    intersections.push('TP53, BRCA1');
    significants.push(fdr < 0.05);
  }

  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Source', sources),
    DG.Column.fromStrings('Term ID', termIds),
    DG.Column.fromStrings('Term Name', termNames),
    DG.Column.fromFloat32Array('P-value', new Float32Array(pValues)),
    DG.Column.fromFloat32Array('FDR', new Float32Array(fdrs)),
    DG.Column.fromInt32Array('Gene Count', new Int32Array(geneCounts)),
    DG.Column.fromFloat32Array('Gene Ratio', new Float32Array(geneRatios)),
    DG.Column.fromStrings('Intersection', intersections),
    DG.Column.fromList('bool', 'Significant', significants),
  ]);
  df.setTag('proteomics.enrichment', 'true');
  return df;
}

/** Creates a mock protein DataFrame with gene symbols.
 * Includes log2FC AND adj.p-value so createVolcanoPlot's column requirements
 * (pickMetricColumn → adj.p-value) are satisfied. Tests B and C below construct
 * a co-located volcano via openEnrichmentVisualization → createVolcanoPlot;
 * without adj.p-value, the try/catch guard would swallow the dock and Tests B/C
 * would falsely report missing layout instead of asserting layout presence. */
function makeMockProteinDf(): DG.DataFrame {
  const genes = ['TP53', 'BRCA1', 'EGFR', 'MYC', 'AKT1'];
  const geneCol = DG.Column.fromStrings('Gene Name', genes);
  geneCol.semType = SEMTYPE.GENE_SYMBOL;
  const df = DG.DataFrame.fromColumns([
    geneCol,
    DG.Column.fromFloat32Array('log2FC', new Float32Array([1.5, -2.0, 0.5, 3.0, -1.0])),
    // Bulk-init via fromFloat32Array — memory feedback_dg_column_bulk_init.
    DG.Column.fromFloat32Array(
      'adj.p-value', new Float32Array([0.001, 0.01, 0.5, 0.0001, 0.04])),
  ]);
  return df;
}

/** Adds a Direction column with the supplied Up/Down assignments. The
 * directional-split branch of openEnrichmentVisualization keys on this. */
function addDirectionColumn(df: DG.DataFrame, dirs: ('Up' | 'Down')[]): void {
  df.columns.add(DG.Column.fromStrings('Direction', dirs as string[]));
}

/** Counts viewers on the supplied table view bound to `target` (same Dart-handle
 * equality the platform uses). */
function countViewersBoundTo(tv: DG.TableView, target: DG.DataFrame): number {
  let n = 0;
  for (const v of tv.viewers) {
    if ((v.dataFrame as any)?.dart === (target as any).dart) n++;
  }
  return n;
}

/** Counts all viewers on the supplied table view. tv.viewers is Iterable, not
 * an array, so we cannot read .length directly — iterate to count instead. */
function countAllViewers(tv: DG.TableView): number {
  let n = 0;
  for (const _ of tv.viewers) n++;
  return n;
}

/** Counts `true` cells in a boolean column without relying on stats semantics. */
function trueCount(col: DG.Column): number {
  let n = 0;
  for (let i = 0; i < col.length; i++) {
    if (col.get(i) === true) n++;
  }
  return n;
}

// --- Tests ---

category('Enrichment Visualization', () => {
  test('prepareEnrichmentChartColumns adds negLog10FDR on the same frame', async () => {
    const df = makeMockEnrichmentDf(20);
    prepareEnrichmentChartColumns(df, 15);
    const negLogCol = df.col(NEG_LOG10_FDR_COL);
    expect(negLogCol !== null, true);
    // Row 0 FDR = 0.001, -log10(0.001) = 3.0
    const val = negLogCol!.get(0) as number;
    expect(Math.abs(val - 3.0) < 0.01, true);
    // The charts bind to the passed frame — no cloned subset is produced.
    expect(df.rowCount, 20);
  });

  test('prepareEnrichmentChartColumns marks exactly topN terms (no Direction)', async () => {
    const df = makeMockEnrichmentDf(20);
    prepareEnrichmentChartColumns(df, 15);
    const marker = df.col(CHART_TOP_MARKER_COL)!;
    expect(trueCount(marker), 15);
  });

  test('prepareEnrichmentChartColumns marks topN per direction', async () => {
    const df = makeMockEnrichmentDf(20);
    // First 10 rows Up, last 10 Down (rows are FDR-ascending by construction).
    addDirectionColumn(df, Array.from({length: 20}, (_, i) => (i < 10 ? 'Up' : 'Down')));
    prepareEnrichmentChartColumns(df, 5);
    const marker = df.col(CHART_TOP_MARKER_COL)!;
    // 5 Up + 5 Down = 10 marked terms.
    expect(trueCount(marker), 10);
  });

  test('prepareEnrichmentChartColumns marks all when fewer rows than topN', async () => {
    const df = makeMockEnrichmentDf(5);
    prepareEnrichmentChartColumns(df, 15);
    const marker = df.col(CHART_TOP_MARKER_COL)!;
    expect(trueCount(marker), 5);
  });

  test('prepareEnrichmentChartColumns is idempotent on re-run', async () => {
    const df = makeMockEnrichmentDf(20);
    prepareEnrichmentChartColumns(df, 15);
    const colsAfterFirst = df.columns.length;
    prepareEnrichmentChartColumns(df, 15);
    // Re-run replaces (ensureFreshFloat pattern) rather than duplicating.
    expect(df.columns.length, colsAfterFirst);
  });

  test('wireEnrichmentToVolcano returns EMPTY without gene column', async () => {
    const enrichDf = makeMockEnrichmentDf(3);
    // Protein df without gene symbol semtype
    const proteinDf = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('SomeColumn', ['A', 'B', 'C']),
    ]);
    const sub = wireEnrichmentToVolcano(enrichDf, proteinDf);
    expect(sub === rxjs.Subscription.EMPTY, true);
  });

  test('wireEnrichmentToVolcano selects a term\'s member genes on the volcano', async () => {
    const enrichDf = makeMockEnrichmentDf(3);
    // Set intersection to specific genes for first row
    enrichDf.col('Intersection')!.set(0, 'TP53, BRCA1');
    const proteinDf = makeMockProteinDf();

    const sub = wireEnrichmentToVolcano(enrichDf, proteinDf);
    expect(sub !== rxjs.Subscription.EMPTY, true);

    // Select the first enrichment term (selection channel, not current row).
    enrichDf.selection.set(0, true, false);
    enrichDf.selection.fireChanged();

    // Check that TP53 (row 0) and BRCA1 (row 1) are selected
    expect(proteinDf.selection.get(0), true);  // TP53
    expect(proteinDf.selection.get(1), true);  // BRCA1
    expect(proteinDf.selection.get(2), false); // EGFR
    expect(proteinDf.selection.get(3), false); // MYC
    expect(proteinDf.selection.get(4), false); // AKT1

    sub.unsubscribe();
  });

  test('wireEnrichmentToVolcano unions member genes across multiple selected terms', async () => {
    const enrichDf = makeMockEnrichmentDf(3);
    enrichDf.col('Intersection')!.set(0, 'TP53');
    enrichDf.col('Intersection')!.set(1, 'EGFR');
    const proteinDf = makeMockProteinDf();

    const sub = wireEnrichmentToVolcano(enrichDf, proteinDf);

    // Select two terms at once (the multi-select gesture the volcano must reflect).
    enrichDf.selection.set(0, true, false);
    enrichDf.selection.set(1, true, false);
    enrichDf.selection.fireChanged();

    // Union: TP53 (row 0) from term 0, EGFR (row 2) from term 1; nothing else.
    expect(proteinDf.selection.get(0), true);  // TP53
    expect(proteinDf.selection.get(2), true);  // EGFR
    expect(proteinDf.selection.get(1), false); // BRCA1 — in neither term
    expect(proteinDf.selection.trueCount, 2);

    sub.unsubscribe();
  });

  test('wireEnrichmentToVolcano clears prior volcano selection on a new term selection', async () => {
    const enrichDf = makeMockEnrichmentDf(3);
    enrichDf.col('Intersection')!.set(0, 'TP53');
    enrichDf.col('Intersection')!.set(1, 'EGFR');
    const proteinDf = makeMockProteinDf();

    // Pre-select row 4 (AKT1)
    proteinDf.selection.set(4, true, false);
    proteinDf.selection.fireChanged();

    const sub = wireEnrichmentToVolcano(enrichDf, proteinDf);

    // Select first enrichment term (TP53 only)
    enrichDf.selection.set(0, true, false);
    enrichDf.selection.fireChanged();

    // AKT1 (row 4) should no longer be selected since selection was cleared
    expect(proteinDf.selection.get(4), false);
    // TP53 (row 0) should be selected
    expect(proteinDf.selection.get(0), true);

    sub.unsubscribe();
  });

  test('wireEnrichmentToVolcano clears the volcano selection when no terms are selected', async () => {
    const enrichDf = makeMockEnrichmentDf(3);
    enrichDf.col('Intersection')!.set(0, 'TP53');
    const proteinDf = makeMockProteinDf();

    const sub = wireEnrichmentToVolcano(enrichDf, proteinDf);

    enrichDf.selection.set(0, true, false);
    enrichDf.selection.fireChanged();
    expect(proteinDf.selection.trueCount, 1);

    // Deselecting every term clears the volcano selection.
    enrichDf.selection.setAll(false, false);
    enrichDf.selection.fireChanged();
    expect(proteinDf.selection.trueCount, 0);

    sub.unsubscribe();
  });

  test('wireEnrichmentToVolcano current term selects its proteins when no terms are selected', async () => {
    const enrichDf = makeMockEnrichmentDf(3);
    enrichDf.col('Intersection')!.set(0, 'TP53, BRCA1');
    const proteinDf = makeMockProteinDf();

    const sub = wireEnrichmentToVolcano(enrichDf, proteinDf);

    // No terms selected → making a term the current row selects its proteins
    // (the visible single-click behavior).
    enrichDf.currentRowIdx = 0;

    expect(proteinDf.selection.get(0), true);  // TP53
    expect(proteinDf.selection.get(1), true);  // BRCA1
    expect(proteinDf.selection.trueCount, 2);

    sub.unsubscribe();
  });

  test('wireEnrichmentToVolcano current term is ignored while a term selection is active', async () => {
    const enrichDf = makeMockEnrichmentDf(3);
    enrichDf.col('Intersection')!.set(0, 'TP53');       // term 0 → TP53 (row 0)
    enrichDf.col('Intersection')!.set(1, 'EGFR, MYC');  // term 1 → EGFR (row 2), MYC (row 3)
    const proteinDf = makeMockProteinDf();

    const sub = wireEnrichmentToVolcano(enrichDf, proteinDf);

    // Active selection on term 1 → EGFR + MYC selected.
    enrichDf.selection.set(1, true, false);
    enrichDf.selection.fireChanged();
    expect(proteinDf.selection.trueCount, 2);
    expect(proteinDf.selection.get(2), true); // EGFR
    expect(proteinDf.selection.get(3), true); // MYC

    // Changing the current row to term 0 must NOT clobber the active selection.
    enrichDf.currentRowIdx = 0;
    expect(proteinDf.selection.get(0), false); // TP53 not added
    expect(proteinDf.selection.get(2), true);  // EGFR still selected
    expect(proteinDf.selection.trueCount, 2);  // unchanged

    sub.unsubscribe();
  });

  // Test B (light-touch declutter, 2026-06): the directional-split Up/Down
  // dot+bar viewers now dock on the ENRICHMENT TableView, NOT the protein view —
  // the user's volcano tab stays clean. This reverses the 13-09 co-dock-on-
  // protein decision per the declutter request.
  test('openEnrichmentVisualization docks dot/bar on enrichment view (directional)', async () => {
    const enrichDf = makeMockEnrichmentDf(4);
    addDirectionColumn(enrichDf, ['Up', 'Up', 'Down', 'Down']);
    enrichDf.name = 'Enrichment (test B)';
    const proteinDf = makeMockProteinDf();
    proteinDf.name = 'Protein (test B)';

    const proteinTv = grok.shell.addTableView(proteinDf);
    const beforeViewerCount = countAllViewers(proteinTv);
    const enrichTv = openEnrichmentVisualization(enrichDf, proteinDf);
    try {
      // The protein/volcano view must NOT gain the enrichment charts.
      expect(countAllViewers(proteinTv), beforeViewerCount,
        `protein view must not gain dot/bar viewers; gained ${countAllViewers(proteinTv) - beforeViewerCount}`);

      // The 4 charts (Up dot + Up bar + Down dot + Down bar) live on the
      // enrichment view alongside its grid (grid + 4 charts ⇒ ≥4).
      expect(countAllViewers(enrichTv) >= 4, true,
        `enrichment view should carry ≥4 dot/bar viewers; got ${countAllViewers(enrichTv)}`);

      // No proteinDf-bound viewer should live on the enrichment view.
      expect(countViewersBoundTo(enrichTv, proteinDf), 0,
        'enrichment view must not carry a proteinDf-bound volcano');

      // Data-layer cross-link still fires: selecting a term marks its proteins.
      enrichDf.selection.set(0, true, false);
      enrichDf.selection.fireChanged();
      expect(proteinDf.selection.trueCount >= 1, true,
        'selecting an enrichment term should mark matching protein rows');
    } finally {
      enrichTv.close();
      proteinTv.close();
    }
  });

  // Test C (declutter): single-direction fallback path also docks on the enrichment view.
  test('openEnrichmentVisualization docks dot/bar on enrichment view (single-direction)', async () => {
    const enrichDf = makeMockEnrichmentDf(3);
    // Intentionally NO Direction column — exercises the fallback branch.
    enrichDf.name = 'Enrichment (test C)';
    const proteinDf = makeMockProteinDf();
    proteinDf.name = 'Protein (test C)';

    const proteinTv = grok.shell.addTableView(proteinDf);
    const beforeViewerCount = countAllViewers(proteinTv);
    const enrichTv = openEnrichmentVisualization(enrichDf, proteinDf);
    try {
      // Protein/volcano view stays clean.
      expect(countAllViewers(proteinTv), beforeViewerCount,
        `protein view must not gain dot/bar viewers; gained ${countAllViewers(proteinTv) - beforeViewerCount}`);
      // Single-direction fallback adds 2 charts (dot + bar) on the enrichment view.
      expect(countAllViewers(enrichTv) >= 2, true,
        `enrichment view should carry ≥2 dot/bar viewers; got ${countAllViewers(enrichTv)}`);

      expect(countViewersBoundTo(enrichTv, proteinDf), 0,
        'enrichment view must not carry a proteinDf-bound volcano');
    } finally {
      enrichTv.close();
      proteinTv.close();
    }
  });

  // Test D (declutter): focus now lands on the ENRICHMENT view (its dot/bar
  // charts), leaving the volcano untouched on the protein tab.
  test('openEnrichmentVisualization switches focus to enrichment TableView', async () => {
    const enrichDf = makeMockEnrichmentDf(4);
    addDirectionColumn(enrichDf, ['Up', 'Up', 'Down', 'Down']);
    enrichDf.name = 'Enrichment (test D)';
    const proteinDf = makeMockProteinDf();
    proteinDf.name = 'Protein (test D)';

    const proteinTv = grok.shell.addTableView(proteinDf);
    const enrichTv = openEnrichmentVisualization(enrichDf, proteinDf);
    try {
      // Duck-type via Dart-handle identity — toJs returns fresh JS wrappers.
      // grok.shell.v is a ViewBase (no dataFrame field on the base type), so
      // cast to any before reaching for it.
      const focusedDart = ((grok.shell.v as any)?.dataFrame as any)?.dart;
      expect(focusedDart === (enrichDf as any).dart, true,
        'focused view should be the enrichment TableView after enrichment runs');
    } finally {
      enrichTv.close();
      proteinTv.close();
    }
  });

  // --- Smart pathway filter banner (14-05 D-14) ---

  test('enrichmentBannerRendersWhenTagSet', async () => {
    const enrichDf = makeMockEnrichmentDf(3);
    enrichDf.name = 'Enrichment (banner-on)';
    enrichDf.setTag('proteomics.enrichment_smart_filtered', 'true');
    enrichDf.setTag('proteomics.enrichment_smart_filtered_kept', '20');
    enrichDf.setTag('proteomics.enrichment_smart_filtered_total', '47');
    enrichDf.setTag('proteomics.enrichment_smart_filtered_dropped_parents', '3');
    enrichDf.setTag('proteomics.enrichment_smart_filtered_cap', '15');
    const proteinDf = makeMockProteinDf();
    const tv = openEnrichmentVisualization(enrichDf, proteinDf);
    try {
      const banners = tv.root.querySelectorAll('[data-smart-filter-banner="true"]');
      expect(banners.length >= 1, true, 'banner should be docked when tag is set');
      const txt = (banners[0] as HTMLElement).textContent ?? '';
      expect(txt.includes('Smart pathway filter active'), true);
      expect(txt.includes('showing 20 of 47 terms'), true);
      expect(txt.includes('dropped 3 generic parents'), true);
      expect(txt.includes('capped at 15 per source'), true);
      expect(txt.includes('Re-run with the filter off'), true);
    } finally {
      tv.close();
    }
  });

  test('enrichmentBannerAbsentWhenTagUnset', async () => {
    const enrichDf = makeMockEnrichmentDf(3);
    enrichDf.name = 'Enrichment (banner-off)';
    // No smart-filter tags set.
    const proteinDf = makeMockProteinDf();
    const tv = openEnrichmentVisualization(enrichDf, proteinDf);
    try {
      const banners = tv.root.querySelectorAll('[data-smart-filter-banner="true"]');
      expect(banners.length, 0);
    } finally {
      tv.close();
    }
  });
});
