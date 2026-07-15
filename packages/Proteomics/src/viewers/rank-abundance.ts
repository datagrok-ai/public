import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {AbundanceByCondition, ConditionAbundance, findAbundanceByCondition} from '../utils/abundance-detection';

/** Percentiles driving the LOD line and the dynamic-range band. */
const LOD_PCT = 0.05;
const BAND_LO_PCT = 0.25;
const BAND_HI_PCT = 0.75;

/** Column-name prefixes for the derived per-condition plot columns. The log10
 * columns are viewer VALUE bindings so they stay plain (a `~` name can break a
 * value binding); both are `.hidden`-tagged to keep the protein grid clean. */
const LOG10_PREFIX = 'log10 intensity: ';
const RANK_PREFIX = 'rank: ';

function ensureFresh(df: DG.DataFrame, name: string, type: DG.ColumnType): DG.Column {
  if (df.columns.contains(name)) df.columns.remove(name);
  return type === DG.COLUMN_TYPE.INT ? df.columns.addNewInt(name) : df.columns.addNewFloat(name);
}

/** rank 1 = highest abundance; unquantified rows (null) sort last and get 0. */
function ranksByAbundance(log10: (number | null)[]): Int32Array {
  const idx = log10.map((_, i) => i)
    .filter((i) => log10[i] != null)
    .sort((a, b) => (log10[b] as number) - (log10[a] as number));
  const ranks = new Int32Array(log10.length); // 0 = not quantified
  idx.forEach((origIdx, pos) => { ranks[origIdx] = pos + 1; });
  return ranks;
}

/** nearest-rank percentile of the non-null values. */
function percentile(values: (number | null)[], p: number): number | null {
  const sorted = values.filter((v): v is number => v != null).sort((a, b) => a - b);
  if (sorted.length === 0) return null;
  return sorted[Math.min(sorted.length - 1, Math.floor(p * (sorted.length - 1)))];
}

interface ConditionCols {
  name: string;
  rankCol: string;
  log10Col: string;
  lod: number | null;
  bandLo: number | null;
  bandHi: number | null;
}

/**
 * Adds the derived per-condition `log10 intensity` and `rank` columns onto `df`
 * (idempotent), and returns the plot-column names + LOD/band thresholds. The
 * scatters bind to `df` itself so their current row is shared with the volcano
 * — clicking a protein anywhere highlights it in the abundance plot for free.
 */
function prepareConditionColumns(df: DG.DataFrame, cond: ConditionAbundance): ConditionCols {
  const log10Name = LOG10_PREFIX + cond.name;
  const rankName = RANK_PREFIX + cond.name;

  const log10Col = ensureFresh(df, log10Name, DG.COLUMN_TYPE.FLOAT);
  log10Col.init((i) => cond.log10Intensity[i]);
  log10Col.setTag('.hidden', 'true');

  const ranks = ranksByAbundance(cond.log10Intensity);
  const rankCol = ensureFresh(df, rankName, DG.COLUMN_TYPE.INT);
  rankCol.init((i) => ranks[i]);
  rankCol.setTag('.hidden', 'true');

  return {
    name: cond.name,
    rankCol: rankName,
    log10Col: log10Name,
    lod: percentile(cond.log10Intensity, LOD_PCT),
    bandLo: percentile(cond.log10Intensity, BAND_LO_PCT),
    bandHi: percentile(cond.log10Intensity, BAND_HI_PCT),
  };
}

/** Removes any rank-abundance formula lines this module previously added (keyed
 * by the `log10 intensity:` column reference) so re-runs don't stack lines. */
function clearRankAbundanceLines(df: DG.DataFrame): void {
  const fl: any = (df as any).meta?.formulaLines;
  if (!fl) return;
  try {
    const items: any[] = Array.isArray(fl.items) ? fl.items : [];
    fl.items = items.filter((line) =>
      !(typeof line?.formula === 'string' && line.formula.includes(`\${${LOG10_PREFIX}`)));
  } catch { /* best effort */ }
}

function addLine(df: DG.DataFrame, yCol: string, value: number, color: string, style: string): void {
  try {
    (df as any).meta.formulaLines.addLine(
      {formula: `\${${yCol}} = ${value}`, color, width: 1, style, visible: true});
  } catch { /* best effort */ }
}

/** Creates one condition's rank-abundance scatter (x = rank, y = log10 intensity)
 * with the LOD line + dynamic-range band, current-row highlighting on. */
export function createRankAbundanceScatter(df: DG.DataFrame, c: ConditionCols): DG.ScatterPlotViewer {
  if (c.lod != null) addLine(df, c.log10Col, c.lod, '#d62728', 'dashed'); // LOD
  if (c.bandLo != null) addLine(df, c.log10Col, c.bandLo, '#2ca02c', 'dotted');
  if (c.bandHi != null) addLine(df, c.log10Col, c.bandHi, '#2ca02c', 'dotted');

  return DG.Viewer.scatterPlot(df, {
    x: c.rankCol,
    y: c.log10Col,
    showViewerFormulaLines: true,
    showCurrentRowMarker: true,
    markerType: 'circle',
    markerDefaultSize: 5,
    xAxisType: 'linear',
    showRegressionLine: false,
    title: `Rank–abundance — ${c.name}`,
  } as any);
}

/**
 * Docks the two per-condition rank-abundance (dynamic-range) S-curves side by
 * side onto `view`, both bound to `df`. Returns true when charts were docked,
 * false when the frame carries no resolvable abundance (caller warns).
 */
export function dockRankAbundanceCharts(view: DG.TableView, df: DG.DataFrame): boolean {
  const abundance: AbundanceByCondition | null = findAbundanceByCondition(df);
  if (!abundance) return false;

  clearRankAbundanceLines(df);
  const c1 = prepareConditionColumns(df, abundance.group1);
  const c2 = prepareConditionColumns(df, abundance.group2);

  const node1 = view.dockManager.dock(
    createRankAbundanceScatter(df, c1), DG.DOCK_TYPE.RIGHT, null, `Abundance — ${c1.name}`, 0.5);
  view.dockManager.dock(
    createRankAbundanceScatter(df, c2), DG.DOCK_TYPE.DOWN, node1, `Abundance — ${c2.name}`, 0.5);
  return true;
}

/**
 * Menu entry point: docks the rank-abundance charts on the active protein view,
 * or warns when there is no abundance data (Report without groups / older
 * Candidates export). Mirrors the enrichment-charts handler shape.
 */
export function openRankAbundance(df: DG.DataFrame): void {
  const view = grok.shell.tv;
  if (!view || (view.dataFrame as any)?.dart !== (df as any).dart) {
    grok.shell.warning('Open the protein table first, then run Rank–Abundance.');
    return;
  }
  if (!dockRankAbundanceCharts(view, df)) {
    grok.shell.warning(
      'No per-condition abundance found. Rank–Abundance needs the Spectronaut ' +
      'Candidates "AVG Group Quantity" columns, or an annotated Report (Annotate ' +
      'Experiment) with per-sample intensities.');
  }
}
