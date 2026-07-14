import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {debounceTime} from 'rxjs/operators';

import {findAbundanceByCondition} from '../utils/abundance-detection';
import {findColumn} from '../utils/column-detection';
import {ensureDirectionColumn} from './volcano';
import {SEMTYPE, DEFAULT_FC_THRESHOLD, DEFAULT_P_THRESHOLD} from '../utils/proteomics-types';

/**
 * Abundance Correlation viewer — one dot per protein, its abundance in group 1
 * (y) against group 2 (x), on geometric-mean log10 axes. Answers: how does each
 * protein's abundance compare between the two conditions, and how tightly do the
 * conditions track overall? Proteins on the y=x diagonal are unchanged; those off
 * it are differentially abundant. Consolidates the former Group-Mean Correlation
 * and Correlation Plot viewers.
 *
 * Works on both pipeline shapes via `findAbundanceByCondition`: a Spectronaut
 * Candidates file (pre-computed AVG Group Quantity) or an annotated Report
 * (per-sample intensities collapsed by geometric mean). Points are colored by the
 * volcano's `direction` (reused when present so the two plots never disagree); the
 * title carries Pearson r and Spearman rho recomputed live over the current
 * filter. No regression line or fold-change guides — the y=x diagonal plus the
 * `direction` color already carry effect size and the significance call.
 */

/** Hidden log10-abundance axis columns. Deliberately distinct from
 * rank-abundance's `log10 intensity:` prefix — that viewer uses the arithmetic
 * mean, this one the geometric, so they must not share a column and overwrite
 * each other when both are opened on the same table. */
const LOG10_PREFIX = 'log10 abundance: ';
const DIAGONAL_COLOR = '#888888';

/** Bulk-init a hidden float column, replacing any prior one so re-runs don't
 * duplicate. */
function ensureFreshHiddenFloat(df: DG.DataFrame, name: string, values: (number | null)[]): string {
  if (df.columns.contains(name)) df.columns.remove(name);
  const c = df.columns.addNewFloat(name);
  c.init((i) => values[i]);
  c.setTag('.hidden', 'true');
  return name;
}

/** Pearson r over the pairwise-complete rows and the pair count. */
export function pearson(x: (number | null)[], y: (number | null)[]): {r: number; n: number} {
  let n = 0; let sx = 0; let sy = 0; let sxx = 0; let syy = 0; let sxy = 0;
  for (let i = 0; i < x.length; i++) {
    const a = x[i]; const b = y[i];
    if (a == null || b == null || !Number.isFinite(a) || !Number.isFinite(b)) continue;
    n++; sx += a; sy += b; sxx += a * a; syy += b * b; sxy += a * b;
  }
  if (n < 2) return {r: NaN, n};
  const cov = sxy - (sx * sy) / n;
  const vx = sxx - (sx * sx) / n;
  const vy = syy - (sy * sy) / n;
  const denom = Math.sqrt(vx * vy);
  return {r: denom > 0 ? cov / denom : NaN, n};
}

/** Fractional ranks (1-based), ties averaged. Nulls keep their slot as null so
 * the caller can drop the pair. */
function fractionalRanks(v: (number | null)[]): (number | null)[] {
  const idx = v.map((_, i) => i).filter((i) => v[i] != null && Number.isFinite(v[i] as number));
  idx.sort((a, b) => (v[a] as number) - (v[b] as number));
  const out: (number | null)[] = new Array(v.length).fill(null);
  let i = 0;
  while (i < idx.length) {
    let j = i;
    while (j + 1 < idx.length && (v[idx[j + 1]] as number) === (v[idx[i]] as number)) j++;
    const avg = (i + j) / 2 + 1;
    for (let k = i; k <= j; k++) out[idx[k]] = avg;
    i = j + 1;
  }
  return out;
}

/** Spearman rho = Pearson over fractional ranks of the pairwise-complete rows. */
export function spearman(x: (number | null)[], y: (number | null)[]): number {
  return pearson(fractionalRanks(x), fractionalRanks(y)).r;
}

/** Reads two columns over the current filter into aligned (null-preserving) arrays. */
function filteredPairs(df: DG.DataFrame, xName: string, yName: string): {xs: (number | null)[]; ys: (number | null)[]} {
  const xc = df.col(xName); const yc = df.col(yName);
  const xs: (number | null)[] = []; const ys: (number | null)[] = [];
  if (!xc || !yc) return {xs, ys};
  for (let i = 0; i < df.rowCount; i++) {
    if (!df.filter.get(i)) continue;
    xs.push(xc.isNone(i) ? null : (xc.get(i) as number));
    ys.push(yc.isNone(i) ? null : (yc.get(i) as number));
  }
  return {xs, ys};
}

/** The on-canvas description (group comparison + correlation stats) over the
 * current filter. Rendered via the scatterplot's `description` overlay — the
 * `title` property is not painted in a docked viewer (the dock tab supplies the
 * header), so the stats live here where they are actually visible. */
export function correlationDescription(
  g1Name: string, g2Name: string, df: DG.DataFrame, xName: string, yName: string,
): string {
  const {xs, ys} = filteredPairs(df, xName, yName);
  const {r, n} = pearson(xs, ys);
  const rho = spearman(xs, ys);
  const rStr = Number.isFinite(r) ? r.toFixed(2) : 'NA';
  const rhoStr = Number.isFinite(rho) ? rho.toFixed(2) : 'NA';
  return `${g1Name} vs ${g2Name} — Pearson r = ${rStr} · Spearman ρ = ${rhoStr} · n = ${n}`;
}

/**
 * Builds the abundance correlation scatter, or returns null when the frame
 * carries no resolvable per-condition abundance (caller warns). Group naming
 * follows the package (group1 = numerator = y).
 */
export function createAbundanceCorrelation(df: DG.DataFrame): DG.ScatterPlotViewer | null {
  const ab = findAbundanceByCondition(df, {reportMean: 'geometric'});
  if (!ab) return null;

  const yName = ensureFreshHiddenFloat(df, LOG10_PREFIX + ab.group1.name, ab.group1.log10Intensity);
  const xName = ensureFreshHiddenFloat(df, LOG10_PREFIX + ab.group2.name, ab.group2.log10Intensity);

  // Reuse the volcano's existing coloring when present so the two plots agree;
  // otherwise classify with the same thresholds. Best-effort — a frame without
  // log2FC (e.g. abundance-only) simply renders uncolored.
  let colorCol: string | undefined;
  try {
    colorCol = df.col('direction') ? 'direction'
      : ensureDirectionColumn(df, DEFAULT_FC_THRESHOLD, DEFAULT_P_THRESHOLD);
  } catch { colorCol = undefined; }

  // y = x reference diagonal — replace any prior one so re-runs don't stack.
  const diagFormula = `\${${yName}} = \${${xName}}`;
  const fl: any = (df as any).meta?.formulaLines;
  if (fl) {
    try {
      const items: any[] = Array.isArray(fl.items) ? fl.items : [];
      fl.items = items.filter((l) => !(typeof l?.formula === 'string' &&
        l.formula.includes(`\${${xName}}`) && l.formula.includes(`\${${yName}}`)));
      fl.addLine({formula: diagFormula, color: DIAGONAL_COLOR, width: 1, visible: true});
      // addLine defaults an empty title back to the formula, so blank it afterward
      // to suppress the raw `${y} = ${x}` label Datagrok otherwise prints along the line.
      const cur: any[] = Array.isArray(fl.items) ? fl.items : [];
      const added = cur.find((l) => l?.formula === diagFormula);
      if (added) { added.title = ''; fl.items = cur; }
    } catch { /* best effort */ }
  }

  const sp = DG.Viewer.scatterPlot(df, {
    x: xName,
    y: yName,
    ...(colorCol ? {color: colorCol} : {}),
    showRegressionLine: false,
    showViewerFormulaLines: true,
    showFilteredOutPoints: true,
    markerType: 'circle',
    markerDefaultSize: 5,
    title: 'Abundance Correlation',
    // Stats go in the always-on description overlay (the title isn't painted).
    description: correlationDescription(ab.group1.name, ab.group2.name, df, xName, yName),
    descriptionVisibilityMode: 'Always',
  } as any);

  // Label with Display Name, gene fallback (matches Group-Mean / Plan 02).
  const labelCol = findColumn(df, SEMTYPE.DISPLAY_NAME, ['display name']) ??
    findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol']);
  if (labelCol) {
    sp.props.labelColumnNames = [labelCol.name];
    sp.props.showLabelsFor = 'MouseOverRow';
  }

  // Filter-live stats: recompute r/rho on filter changes; tear the subscription
  // down when this viewer closes so it doesn't leak or write to a dead viewer.
  const sub = df.onFilterChanged.pipe(debounceTime(100)).subscribe(() => {
    try {
      sp.setOptions({description: correlationDescription(ab.group1.name, ab.group2.name, df, xName, yName)});
    } catch { /* viewer closed */ }
  });
  const closeSub = grok.events.onViewerClosed.subscribe((args) => {
    const v: any = (args as any)?.args?.viewer;
    if (v && (v as any).dart === (sp as any).dart) {
      sub.unsubscribe();
      closeSub.unsubscribe();
    }
  });

  return sp;
}

/**
 * Menu entry point: docks the abundance correlation plot on the active protein
 * view, or warns when there is no per-condition abundance.
 */
export function openAbundanceCorrelation(df: DG.DataFrame): void {
  const view = grok.shell.tv;
  if (!view || (view.dataFrame as any)?.dart !== (df as any).dart) {
    grok.shell.warning('Open the protein table first, then run Abundance Correlation.');
    return;
  }
  const sp = createAbundanceCorrelation(df);
  if (!sp) {
    grok.shell.warning(
      'No per-condition abundance found. Abundance Correlation needs the Spectronaut ' +
      'Candidates "AVG Group Quantity" columns, or an annotated Report (Annotate ' +
      'Experiment) with per-sample intensities.');
    return;
  }
  view.dockManager.dock(sp, DG.DOCK_TYPE.RIGHT, null, 'Abundance Correlation', 0.5);
}
