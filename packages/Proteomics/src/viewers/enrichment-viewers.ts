import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {findColumn} from '../utils/column-detection';
import {SEMTYPE} from '../utils/proteomics-types';

/** Module-level subscriptions for cleanup on re-open. */
let activeSubscriptions: rxjs.Subscription[] = [];

/**
 * Creates a top-N enrichment DataFrame sorted by FDR ascending.
 * Adds a -log10(FDR) column for visualization, truncates long term names.
 */
export function createTopNEnrichmentDf(enrichDf: DG.DataFrame, topN: number = 15): DG.DataFrame {
  const fdrCol = enrichDf.col('FDR')!;
  const indexed: {idx: number; fdr: number}[] = [];
  for (let i = 0; i < enrichDf.rowCount; i++) {
    if (!fdrCol.isNone(i))
      indexed.push({idx: i, fdr: fdrCol.get(i) as number});
  }
  indexed.sort((a, b) => a.fdr - b.fdr);

  const count = Math.min(topN, indexed.length);
  const mask = DG.BitSet.create(enrichDf.rowCount);
  for (let i = 0; i < count; i++)
    mask.set(indexed[i].idx, true);

  const topDf = enrichDf.clone(mask);
  topDf.name = `Top ${count} Enriched Terms`;

  // Truncate long term names in place via init (re-reads via getRawData wouldn't
  // help for string columns; reading get(i) here is the simplest correct form).
  const termCol = topDf.col('Term Name');
  if (termCol) {
    const orig: (string | null)[] = [];
    for (let i = 0; i < topDf.rowCount; i++) orig.push(termCol.get(i) as string | null);
    termCol.init((i) => {
      const val = orig[i];
      return val && val.length > 50 ? val.substring(0, 50) + '...' : val;
    });
  }

  // Add -log10(FDR) column for color mapping. FDR=0 (underflow / extreme
  // enrichment) plots at the float-underflow ceiling so it stays visible but
  // doesn't collide with real values like FDR=1e-10.
  const UNDERFLOW_NEGLOG10 = -Math.log10(Number.MIN_VALUE);
  const topFdrRaw = topDf.col('FDR')!.getRawData() as Float32Array | Float64Array;
  const negLogCol = topDf.columns.addNewFloat('negLog10FDR');
  negLogCol.init((i) => topFdrRaw[i] > 0 ? -Math.log10(topFdrRaw[i]) : UNDERFLOW_NEGLOG10);

  return topDf;
}

/**
 * Creates an enrichment dot plot using a scatter plot viewer.
 * X = Gene Ratio, Y = Term Name, size = Gene Count, color = -log10(FDR).
 */
export function createEnrichmentDotPlot(topDf: DG.DataFrame): DG.ScatterPlotViewer {
  const sp = DG.Viewer.scatterPlot(topDf, {
    x: 'Gene Ratio',
    y: 'Term Name',
    sizeColumnName: 'Gene Count',
    colorColumnName: 'negLog10FDR',
    markerMinSize: 5,
    markerMaxSize: 25,
    title: 'Enrichment Dot Plot',
  });
  return sp;
}

/**
 * Creates an enrichment bar chart showing top terms ranked by -log10(FDR).
 */
export function createEnrichmentBarChart(topDf: DG.DataFrame): DG.Viewer {
  const bar = DG.Viewer.barChart(topDf, {
    splitColumnName: 'Term Name',
    valueColumnName: 'negLog10FDR',
    valueAggrType: 'avg',
    barSortType: 'by value',
    barSortOrder: 'desc',
    orientation: 'Horizontal',
    title: 'Top Enriched Terms',
  } as any);
  return bar;
}

/**
 * Wires enrichment row changes to protein DataFrame selection.
 * When a row is selected in the enrichment table, matching protein rows
 * (by gene symbol from the Intersection column) are highlighted.
 */
export function wireEnrichmentToVolcano(
  enrichDf: DG.DataFrame,
  proteinDf: DG.DataFrame,
): rxjs.Subscription {
  const geneCol = findColumn(proteinDf, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol']);
  if (!geneCol) return rxjs.Subscription.EMPTY;

  // Build gene-to-row index
  const geneToRows = new Map<string, number[]>();
  for (let i = 0; i < proteinDf.rowCount; i++) {
    if (!geneCol.isNone(i)) {
      const gene = (geneCol.get(i) as string).trim();
      if (!geneToRows.has(gene)) geneToRows.set(gene, []);
      geneToRows.get(gene)!.push(i);
    }
  }

  return enrichDf.onCurrentRowChanged.subscribe(() => {
    const rowIdx = enrichDf.currentRowIdx;
    if (rowIdx < 0) return;

    const intersectionCol = enrichDf.col('Intersection');
    if (!intersectionCol) return;

    const memberGenesStr = intersectionCol.get(rowIdx) as string;
    if (!memberGenesStr) return;

    const memberGenes = memberGenesStr.split(',').map((g) => g.trim()).filter(Boolean);

    // Clear previous selection
    proteinDf.selection.setAll(false, false);

    // Set selection for matching genes
    for (const gene of memberGenes) {
      const rows = geneToRows.get(gene);
      if (rows) {
        for (const row of rows)
          proteinDf.selection.set(row, true, false);
      }
    }
    proteinDf.selection.fireChanged();
  });
}

/**
 * Clones the enrichment DataFrame down to one direction via the same
 * clone-by-mask idiom as createTopNEnrichmentDf. Returns null when the
 * Direction column is absent or the direction has no rows.
 */
function filterByDirection(enrichDf: DG.DataFrame, direction: 'Up' | 'Down'): DG.DataFrame | null {
  const dirCol = enrichDf.col('Direction');
  if (!dirCol) return null;
  const mask = DG.BitSet.create(enrichDf.rowCount);
  let any = false;
  for (let i = 0; i < enrichDf.rowCount; i++) {
    if (dirCol.get(i) === direction) { mask.set(i, true); any = true; }
  }
  if (!any) return null;
  const sub = enrichDf.clone(mask);
  sub.name = `${direction} Enrichment`;
  return sub;
}

/**
 * Opens a full enrichment visualization dashboard with dot plot, bar chart,
 * and cross-DataFrame selection wiring to the protein table. When the result
 * carries an R2 Direction column with both Up and Down rows, the Up and Down
 * dot+bar charts are docked side-by-side; otherwise the original single
 * dot+bar layout is used. The Phase-9 volcano cross-link is unchanged.
 */
export function openEnrichmentVisualization(
  enrichDf: DG.DataFrame,
  proteinDf: DG.DataFrame,
): DG.TableView {
  // Clean up previous subscriptions
  for (const sub of activeSubscriptions)
    sub.unsubscribe();
  activeSubscriptions = [];

  // Find or reuse the enrichment results table view. Match by stable .dart
  // identity, NOT `=== enrichDf`: grok.shell.tableViews returns a fresh toJs
  // wrapper per access, so strict-equality misses an already-open view and
  // spawns a duplicate "Enrichment Results (2)" tab (precedent: the protein-TV
  // lookup that used to live here; 13-07 SUMMARY).
  let existing: DG.TableView | undefined;
  for (const v of grok.shell.tableViews) {
    if ((v.dataFrame as any)?.dart === (enrichDf as any).dart) { existing = v; break; }
  }
  const tv = existing ?? grok.shell.addTableView(enrichDf);

  // Dock the dot/bar charts onto the ENRICHMENT results view — not the protein
  // view where the volcano lives — so the volcano tab stays uncluttered
  // (light-touch declutter). Surface this view FIRST: docking onto a
  // backgrounded host leaves zero-dimension dock splitters and crashes the Dart
  // scatter plot's smart-labels with
  //   "Failed to execute 'drawImage'... canvas with width or height of 0".
  grok.shell.v = tv;
  const chartHost = tv;

  // D-14 — surface the smart-filter transform stats above the grid when it
  // actually dropped rows. Tag is set only by runEnrichmentPipeline when
  // kept < total; the banner copy is locked in 14-UI-SPEC.md §"Copywriting".
  if (enrichDf.getTag('proteomics.enrichment_smart_filtered') === 'true') {
    const kept = enrichDf.getTag('proteomics.enrichment_smart_filtered_kept') ?? '?';
    const total = enrichDf.getTag('proteomics.enrichment_smart_filtered_total') ?? '?';
    const parents = enrichDf.getTag('proteomics.enrichment_smart_filtered_dropped_parents') ?? '?';
    const cap = enrichDf.getTag('proteomics.enrichment_smart_filtered_cap') ?? '15';
    const banner = ui.divText(
      `Smart pathway filter active: showing ${kept} of ${total} terms ` +
      `(dropped ${parents} generic parents, capped at ${cap} per source). ` +
      `Re-run with the filter off to see all.`);
    banner.dataset.smartFilterBanner = 'true';
    tv.dockManager.dock(banner, DG.DOCK_TYPE.TOP, null, 'Smart Filter Info', 0.1);
  }

  // Dock the dot/bar charts (extracted so the publish path reuses the exact
  // 2×2) and cross-link the full frame plus each chart-backing subset. All keep
  // the Intersection column the Phase-9 wiring keys on.
  const subsetDfs = dockEnrichmentCharts(chartHost, enrichDf);
  activeSubscriptions.push(wireEnrichmentToVolcano(enrichDf, proteinDf));
  for (const sub of subsetDfs)
    activeSubscriptions.push(wireEnrichmentToVolcano(sub, proteinDf));

  // The dot/bar charts, the merged enrichment grid, and the smart-filter banner
  // all share this enrichment tab; the volcano stays on its own (protein) tab.
  // The cross-link keeps the two tabs in sync on selection.
  return tv;
}

/** Docks the directional (Up|Down 2×2) or single Dot+Bar enrichment charts onto
 * `host` and returns the top-N subset DataFrames they bind to. Pure layout — no
 * subscriptions, no focus changes — so both the live viewer and the publish path
 * (which bundles the returned subsets as project children) reuse one code path. */
export function dockEnrichmentCharts(host: DG.TableView, enrichDf: DG.DataFrame): DG.DataFrame[] {
  const upSub = filterByDirection(enrichDf, 'Up');
  const downSub = filterByDirection(enrichDf, 'Down');

  if (upSub && downSub) {
    const upTop = createTopNEnrichmentDf(upSub, 15);
    const downTop = createTopNEnrichmentDf(downSub, 15);
    // Clear names for the publish tree (the live path doesn't surface these).
    upTop.name = 'Enrichment — Up (top terms)';
    downTop.name = 'Enrichment — Down (top terms)';

    // Balanced 2×2: [Up Dot | Down Dot] over [Up Bar | Down Bar]. Split the
    // chart region into the two columns FIRST, THEN drop each bar under its dot
    // — docking the bars before the split carved up only the top-left cell and
    // produced a lopsided layout.
    const upDotNode = host.dockManager.dock(
      createEnrichmentDotPlot(upTop), DG.DOCK_TYPE.RIGHT, null, 'Up — Dot Plot', 0.6);
    const downDotNode = host.dockManager.dock(
      createEnrichmentDotPlot(downTop), DG.DOCK_TYPE.RIGHT, upDotNode, 'Down — Dot Plot', 0.5);
    host.dockManager.dock(
      createEnrichmentBarChart(upTop), DG.DOCK_TYPE.DOWN, upDotNode, 'Up — Bar Chart', 0.5);
    host.dockManager.dock(
      createEnrichmentBarChart(downTop), DG.DOCK_TYPE.DOWN, downDotNode, 'Down — Bar Chart', 0.5);
    return [upTop, downTop];
  }

  // Single-direction or no Direction column → one Dot + Bar.
  const topDf = createTopNEnrichmentDf(enrichDf, 15);
  topDf.name = 'Enrichment — Top terms';
  const dotNode = host.dockManager.dock(
    createEnrichmentDotPlot(topDf), DG.DOCK_TYPE.RIGHT, null, 'Dot Plot', 0.5);
  host.dockManager.dock(
    createEnrichmentBarChart(topDf), DG.DOCK_TYPE.DOWN, dotNode, 'Bar Chart', 0.5);
  return [topDf];
}
