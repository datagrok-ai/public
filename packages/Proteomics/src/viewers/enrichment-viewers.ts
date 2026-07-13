import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {findColumn} from '../utils/column-detection';
import {SEMTYPE} from '../utils/proteomics-types';

/** Module-level subscriptions for cleanup on re-open. */
let activeSubscriptions: rxjs.Subscription[] = [];

/** Default number of terms each dot/bar chart shows (per direction). */
export const CHART_TOP_N = 15;

/** Value/color column the dot + bar charts read. Kept as a plain (non-`~`)
 * name because it is a viewer VALUE binding — a missing value column hard-fails
 * a bar chart on reopen (the original "Column negLog10FDR does not exist" bug). */
export const NEG_LOG10_FDR_COL = 'negLog10FDR';
/** Hidden marker column (auto-hidden via `~` prefix, ignored by the publish
 * shape contract) flagging the top-N terms per direction. Used only inside the
 * per-viewer formula filter, where a parse failure degrades to "show all"
 * rather than crashing — so the `~` name is safe here but not for the value. */
export const CHART_TOP_MARKER_COL = '~enrichChartTop';

/**
 * Adds the two derived columns the enrichment dot/bar charts bind to, directly
 * onto `enrichDf` (idempotent — ensureFreshFloat pattern):
 *   - `negLog10FDR` — -log10(FDR); dot color + bar value. FDR=0 (underflow /
 *     extreme enrichment) maps to the float-underflow ceiling so it stays
 *     visible without colliding with real values like FDR=1e-10.
 *   - `~enrichChartTop` — true for the top-N terms by FDR (ascending) within
 *     each Direction, or overall when there is no Direction column.
 *
 * Both charts bind to `enrichDf` itself rather than to cloned subsets: docked
 * viewers rebind to their host view's primary frame on project reopen, so the
 * columns they read must live on that frame (else the bar chart fatal-errors).
 * Keeping everything on one frame also stops the subset tables from leaking
 * into a published project's table tree.
 */
export function prepareEnrichmentChartColumns(enrichDf: DG.DataFrame, topN: number = CHART_TOP_N): void {
  const fdrCol = enrichDf.col('FDR');
  if (!fdrCol) return;

  const UNDERFLOW_NEGLOG10 = -Math.log10(Number.MIN_VALUE);
  const fdrRaw = fdrCol.getRawData() as Float32Array | Float64Array;
  if (enrichDf.columns.contains(NEG_LOG10_FDR_COL)) enrichDf.columns.remove(NEG_LOG10_FDR_COL);
  const negLogCol = enrichDf.columns.addNewFloat(NEG_LOG10_FDR_COL);
  negLogCol.init((i) => fdrRaw[i] > 0 ? -Math.log10(fdrRaw[i]) : UNDERFLOW_NEGLOG10);
  // Keep it out of the reviewer's grid (best-effort; `~` isn't usable on a
  // viewer value binding). If the tag is stripped on round-trip the only cost
  // is one extra, meaningful numeric column in the grid.
  negLogCol.setTag('.hidden', 'true');

  // Rank terms by FDR (ascending) within each direction and mark the top N.
  const dirCol = enrichDf.col('Direction');
  const groups = new Map<string, {idx: number; fdr: number}[]>();
  for (let i = 0; i < enrichDf.rowCount; i++) {
    if (fdrCol.isNone(i)) continue;
    const key = dirCol ? String(dirCol.get(i)) : '';
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key)!.push({idx: i, fdr: fdrCol.get(i) as number});
  }
  const topMask = new Array<boolean>(enrichDf.rowCount).fill(false);
  for (const arr of groups.values()) {
    arr.sort((a, b) => a.fdr - b.fdr);
    for (let k = 0; k < Math.min(topN, arr.length); k++)
      topMask[arr[k].idx] = true;
  }
  if (enrichDf.columns.contains(CHART_TOP_MARKER_COL)) enrichDf.columns.remove(CHART_TOP_MARKER_COL);
  enrichDf.columns.addNewBool(CHART_TOP_MARKER_COL).init((i) => topMask[i]);
}

/** Per-viewer formula filter: the top-N marked terms, optionally narrowed to one
 * regulation direction. Degrades to "show all rows" if the platform drops the
 * formula on reopen — never a missing-column crash. */
function chartFilter(direction?: 'Up' | 'Down'): string {
  const clauses = [`\${${CHART_TOP_MARKER_COL}}`];
  if (direction) clauses.push(`\${Direction} == "${direction}"`);
  return clauses.join(' and ');
}

/**
 * Creates an enrichment dot plot bound to the shared enrichment frame.
 * X = Gene Ratio, Y = Term Name, size = Gene Count, color = -log10(FDR).
 * Pass `direction` to split into Up/Down via a per-viewer formula filter.
 */
export function createEnrichmentDotPlot(enrichDf: DG.DataFrame, direction?: 'Up' | 'Down'): DG.ScatterPlotViewer {
  const sp = DG.Viewer.scatterPlot(enrichDf, {
    x: 'Gene Ratio',
    y: 'Term Name',
    sizeColumnName: 'Gene Count',
    colorColumnName: NEG_LOG10_FDR_COL,
    markerMinSize: 5,
    markerMaxSize: 25,
    filter: chartFilter(direction),
    title: direction ? `Enrichment — ${direction}` : 'Enrichment Dot Plot',
  } as any);
  return sp;
}

/**
 * Creates an enrichment bar chart (top terms ranked by -log10(FDR)) bound to the
 * shared enrichment frame. Pass `direction` to split into Up/Down.
 */
export function createEnrichmentBarChart(enrichDf: DG.DataFrame, direction?: 'Up' | 'Down'): DG.Viewer {
  const bar = DG.Viewer.barChart(enrichDf, {
    splitColumnName: 'Term Name',
    valueColumnName: NEG_LOG10_FDR_COL,
    valueAggrType: 'avg',
    barSortType: 'by value',
    barSortOrder: 'desc',
    orientation: 'Horizontal',
    filter: chartFilter(direction),
    title: direction ? `Top Enriched — ${direction}` : 'Top Enriched Terms',
  } as any);
  return bar;
}

/**
 * Wires enrichment term interactions to the protein DataFrame, keyed on each
 * term's comma-separated `Intersection` member genes → protein rows. All three
 * enrichment charts (grid, dot plot, bar chart) share `enrichDf`, so this single
 * full-frame wiring covers interactions in any of them.
 *
 * Everything drives the volcano's `selection` — the only channel the scatterplot
 * paints visibly for a *set* of proteins (verified live: programmatic
 * `rows.highlight` / `mouseOverRowFunc` do not render). Two channels feed that
 * one selection, with the multi-term selection winning:
 *  - (1)+(3) enrichment **selection** (dot-plot rubber-band, shift/ctrl-click
 *    bars, ctrl-click grid rows) → the UNION of member proteins across every
 *    selected term. An empty enrichment selection clears the volcano selection.
 *  - (2) enrichment **current row** (a single clicked term) → that term's member
 *    proteins, but ONLY when no terms are selected. While a multi-term selection
 *    is active it governs, and current-row changes are ignored — so a trailing
 *    current-row event (fired alongside a select gesture) can't clobber the union.
 *
 * The enrichment **filter** is intentionally NOT reflected (4) — filtering terms
 * declutters the enrichment view; it must not repaint the volcano.
 */
export function wireEnrichmentToVolcano(
  enrichDf: DG.DataFrame,
  proteinDf: DG.DataFrame,
): rxjs.Subscription {
  const geneCol = findColumn(proteinDf, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol']);
  if (!geneCol) return rxjs.Subscription.EMPTY;

  // Build gene -> protein-row index once.
  const geneToRows = new Map<string, number[]>();
  for (let i = 0; i < proteinDf.rowCount; i++) {
    if (!geneCol.isNone(i)) {
      const gene = (geneCol.get(i) as string).trim();
      if (!geneToRows.has(gene)) geneToRows.set(gene, []);
      geneToRows.get(gene)!.push(i);
    }
  }

  // Union of protein rows whose gene is a member of ANY of the given enrichment
  // term rows (via each term's comma-separated Intersection string).
  const proteinRowsForTerms = (termRows: Int32Array | number[]): Set<number> => {
    const out = new Set<number>();
    const intersectionCol = enrichDf.col('Intersection');
    if (!intersectionCol) return out;
    for (const r of termRows) {
      if (r < 0 || intersectionCol.isNone(r)) continue;
      const memberGenesStr = intersectionCol.get(r) as string;
      if (!memberGenesStr) continue;
      for (const gene of memberGenesStr.split(',').map((g) => g.trim()).filter(Boolean)) {
        const rows = geneToRows.get(gene);
        if (rows) for (const row of rows) out.add(row);
      }
    }
    return out;
  };

  // Replace the volcano selection with exactly `rows` (empty set → cleared).
  const selectProteins = (rows: Set<number>): void => {
    proteinDf.selection.setAll(false, false);
    for (const row of rows) proteinDf.selection.set(row, true, false);
    proteinDf.selection.fireChanged();
  };

  // (1)+(3) Selected terms → union of member proteins; no selection → cleared.
  const applySelection = (): void => {
    const termRows = enrichDf.selection.getSelectedIndexes();
    selectProteins(termRows.length > 0 ? proteinRowsForTerms(termRows) : new Set());
  };

  // (2) Single current term → that term's member proteins, but only while no
  // terms are selected; an active multi-term selection governs and wins.
  const applyCurrentRow = (): void => {
    if (enrichDf.selection.trueCount > 0) return;
    const idx = enrichDf.currentRowIdx;
    selectProteins(idx >= 0 ? proteinRowsForTerms([idx]) : new Set());
  };

  const sub = new rxjs.Subscription();
  sub.add(enrichDf.onSelectionChanged.subscribe(() => applySelection()));
  sub.add(enrichDf.onCurrentRowChanged.subscribe(() => applyCurrentRow()));
  return sub;
}

/** True when `enrichDf` carries a Direction column with at least one row in
 * `direction`. Drives the Up/Down 2×2 vs single-chart layout choice. */
function hasDirection(enrichDf: DG.DataFrame, direction: 'Up' | 'Down'): boolean {
  const dirCol = enrichDf.col('Direction');
  if (!dirCol) return false;
  for (let i = 0; i < enrichDf.rowCount; i++) {
    if (dirCol.get(i) === direction) return true;
  }
  return false;
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
  // 2×2) and cross-link the frame. Every chart binds to `enrichDf` itself, so
  // the single full-frame wiring covers all of them (and keeps the Intersection
  // column the Phase-9 wiring keys on).
  dockEnrichmentCharts(chartHost, enrichDf);
  activeSubscriptions.push(wireEnrichmentToVolcano(enrichDf, proteinDf));

  // The dot/bar charts, the merged enrichment grid, and the smart-filter banner
  // all share this enrichment tab; the volcano stays on its own (protein) tab.
  // The cross-link keeps the two tabs in sync on selection.
  return tv;
}

/** Docks the directional (Up|Down 2×2) or single Dot+Bar enrichment charts onto
 * `host`. Every chart binds to `enrichDf` itself — the derived `negLog10FDR` /
 * `~enrichChartTop` columns are added to that frame and the Up/Down split is a
 * per-viewer formula filter. No cloned subsets, so nothing leaks into a
 * published project's table tree and the charts survive a project round-trip
 * (docked viewers rebind to their host frame on reopen). Pure layout — no
 * subscriptions, no focus changes — shared by the live and publish paths. */
export function dockEnrichmentCharts(host: DG.TableView, enrichDf: DG.DataFrame): void {
  prepareEnrichmentChartColumns(enrichDf);

  if (hasDirection(enrichDf, 'Up') && hasDirection(enrichDf, 'Down')) {
    // Balanced 2×2: [Up Dot | Down Dot] over [Up Bar | Down Bar]. Split the
    // chart region into the two columns FIRST, THEN drop each bar under its dot
    // — docking the bars before the split carved up only the top-left cell and
    // produced a lopsided layout.
    const upDotNode = host.dockManager.dock(
      createEnrichmentDotPlot(enrichDf, 'Up'), DG.DOCK_TYPE.RIGHT, null, 'Up — Dot Plot', 0.6);
    const downDotNode = host.dockManager.dock(
      createEnrichmentDotPlot(enrichDf, 'Down'), DG.DOCK_TYPE.RIGHT, upDotNode, 'Down — Dot Plot', 0.5);
    host.dockManager.dock(
      createEnrichmentBarChart(enrichDf, 'Up'), DG.DOCK_TYPE.DOWN, upDotNode, 'Up — Bar Chart', 0.5);
    host.dockManager.dock(
      createEnrichmentBarChart(enrichDf, 'Down'), DG.DOCK_TYPE.DOWN, downDotNode, 'Down — Bar Chart', 0.5);
    return;
  }

  // Single-direction or no Direction column → one Dot + Bar.
  const dotNode = host.dockManager.dock(
    createEnrichmentDotPlot(enrichDf), DG.DOCK_TYPE.RIGHT, null, 'Dot Plot', 0.5);
  host.dockManager.dock(
    createEnrichmentBarChart(enrichDf), DG.DOCK_TYPE.DOWN, dotNode, 'Bar Chart', 0.5);
}
