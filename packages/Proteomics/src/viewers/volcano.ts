import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';
import {findColumn} from '../utils/column-detection';
import {SEMTYPE} from '../utils/proteomics-types';
import {parseAccession} from '../panels/uniprot-panel';
import {getGroups} from '../analysis/experiment-setup';
import {
  getSubcellularLocations,
  LOCATION_COLORS,
  ProgressCb,
} from '../analysis/subcellular-location';

/** -log10 of the smallest representable positive double, used as a sentinel for
 * p-values that underflow to literal 0 (e.g. from t-test extreme inputs).
 * Plotting at this height keeps the point visible but clearly off-scale,
 * rather than colliding with real values like p = 1e-10. */
const UNDERFLOW_NEGLOG10 = -Math.log10(Number.MIN_VALUE);

/** Significance metric the volcano Y axis / classification is computed from. */
export type MetricKind = 'adj.p-value' | 'p-value';

/** Active color dimension on the volcano (significance categories vs subcellular location). */
export type ColorDim = 'significance' | 'location';

/** Stable, generic column names. The Y/color bindings keep these names across
 * a metric/color toggle; only the values change in place (Pitfall 5 — never
 * early-return on existence, never remove+re-add a bound column). */
const NEG_LOG_COL = 'negLog10P';
const DIRECTION_COL = 'direction';
const LOCATION_COL = 'Subcellular Location';

/** Persists the active metric on the DataFrame so the Volcano Options dialog
 * (G2) and the axis-label overlay (G1) can derive Y-axis text without
 * re-reading scatter props. recomputeVolcano updates this on every toggle. */
export const VOLCANO_METRIC_TAG = 'proteomics.volcano_metric';

/** Module-level subscription store mirrors `enrichment-viewers.ts:9` —
 * disposed at the start of every createVolcanoPlot to prevent stacked
 * overlays / frozen counters on viewer re-entry (Pitfall 6). */
let activeVolcanoSubscriptions: rxjs.Subscription[] = [];

/** Resolves the significance column for the requested metric. `p-value` falls
 * back to `adj.p-value` when the raw p-value column is absent (D-06 guard — no
 * throw, stays on adj.p-value). Throws only when adj.p-value itself is missing
 * (genuinely no DE result). */
function pickMetricColumn(df: DG.DataFrame, metric: MetricKind): DG.Column {
  if (metric === 'p-value') {
    const p = df.col('p-value');
    if (p) return p;
  }
  const adj = df.col('adj.p-value');
  if (!adj)
    throw new Error('adj.p-value column not found');
  return adj;
}

/** Ensures the -log10(metric) column exists and reflects `metric`. The column
 * name is FIXED so the scatter Y binding is stable; values are ALWAYS re-init
 * in place so a metric toggle updates the axis without rebinding. */
export function ensureNegLog10Column(
  df: DG.DataFrame, metric: MetricKind = 'adj.p-value',
): string {
  const pCol = pickMetricColumn(df, metric);
  const col = df.col(NEG_LOG_COL) ?? df.columns.addNewFloat(NEG_LOG_COL);
  col.init((i) => {
    if (pCol.isNone(i)) return DG.FLOAT_NULL;
    const p = pCol.get(i) as number;
    return p > 0 ? -Math.log10(p) : UNDERFLOW_NEGLOG10;
  });
  return NEG_LOG_COL;
}

/** LOCKED color contract (D-04, CONTEXT.md §"Specifics"): magenta = numerator
 * (group1), cyan = denominator (group2), gray = NS. Keys are derived per call
 * from `getGroups(df)`; the ARGB ints here are the invariant base. */
export const DIRECTION_COLORS_BASE = {
  enrichedG1: 0xFFFF00FF, // magenta (ARGB)
  enrichedG2: 0xFF00FFFF, // cyan (ARGB)
  notSig: 0xFFAAAAAA, // gray (ARGB)
} as const;

/** Ensures the direction column, classified against the chosen `metric` (same
 * column the Y axis uses, so dots and coloring never disagree). Updated in
 * place — never remove+re-add (preserves viewer/grid bindings). D-04 migration:
 * category strings derive from `getGroups(df)` (`Enriched in {g1}` / `Enriched
 * in {g2}` / `Not significant`), color map is the LOCKED magenta/cyan/gray. */
export function ensureDirectionColumn(
  df: DG.DataFrame,
  fcThreshold: number,
  pThreshold: number,
  metric: MetricKind = 'adj.p-value',
): string {
  const fcCol = df.col('log2FC');
  if (!fcCol)
    throw new Error('log2FC column not found');
  const pCol = pickMetricColumn(df, metric);

  const groups = getGroups(df);
  const g1Label = groups ? `Enriched in ${groups.group1.name}` : 'Enriched in group1';
  const g2Label = groups ? `Enriched in ${groups.group2.name}` : 'Enriched in group2';
  const nsLabel = 'Not significant';

  const n = df.rowCount;
  const vals: string[] = new Array(n);
  for (let i = 0; i < n; i++) {
    if (fcCol.isNone(i) || pCol.isNone(i)) { vals[i] = nsLabel; continue; }
    const fc = fcCol.get(i) as number;
    const p = pCol.get(i) as number;
    if (p <= pThreshold && fc > fcThreshold) vals[i] = g1Label;
    else if (p <= pThreshold && fc < -fcThreshold) vals[i] = g2Label;
    else vals[i] = nsLabel;
  }

  const col = df.col(DIRECTION_COL) ?? df.columns.addNewString(DIRECTION_COL);
  col.init((i) => vals[i]);
  col.meta.colors.setCategorical({
    [g1Label]: DIRECTION_COLORS_BASE.enrichedG1,
    [g2Label]: DIRECTION_COLORS_BASE.enrichedG2,
    [nsLabel]: DIRECTION_COLORS_BASE.notSig,
  });
  return DIRECTION_COL;
}

/** FNV-1a 32-bit hex hash. Deterministic, dependency-free, fast enough for the
 * ~8k-accession stable-sorted join that ensureLocationColumn feeds it. Used
 * only as a cache-invalidation key on the Subcellular Location column. */
function fnv1aHex(s: string): string {
  let h = 0x811c9dc5;
  for (let i = 0; i < s.length; i++) {
    h ^= s.charCodeAt(i);
    h = (h + ((h << 1) + (h << 4) + (h << 7) + (h << 8) + (h << 24))) >>> 0;
  }
  return h.toString(16).padStart(8, '0');
}

/** Cache-invalidation key on the Subcellular Location COLUMN (not the
 * DataFrame). This is the first column-level tag in the package — the
 * DataFrame tag convention covers whole-pipeline workflow state, while
 * per-column cache-invalidation metadata belongs on the column. Documented
 * in 13-08 SUMMARY.md as a new convention. */
const LOCATION_HASH_TAG = 'proteomics.location_acc_hash';

/** Ensures the 'Subcellular Location' categorical column using the shared
 * verbatim-CK-omics service (13-04). Bulk-init string column, the LOCKED 11+1
 * palette via setCategorical, semType = SUBCELLULAR_LOCATION. In place.
 *
 * Short-circuits when the column already exists AND its accession-set hash
 * tag matches the current accession set — skipping the entire fetch+init
 * pipeline. The palette and semType survive on the column object across
 * the short-circuit because `col.meta.colors.setCategorical` and
 * `col.semType` are set on the same column object as last run.
 *
 * `progress` is forwarded to getSubcellularLocations for the fetch phases
 * AND fired once for the 'init-column' phase so the dialog's pi can advance
 * past the post-fetch bulk-init. */
export async function ensureLocationColumn(
  df: DG.DataFrame, progress?: ProgressCb,
): Promise<string> {
  const idCol = findColumn(df, SEMTYPE.PROTEIN_ID,
    ['primary protein id', 'protein id', 'protein ids', 'accession', 'uniprot']) ??
    df.col('Primary Protein ID') ?? df.col('Protein ID');

  const n = df.rowCount;
  const accForRow: (string | null)[] = new Array(n).fill(null);
  const accessions = new Set<string>();
  if (idCol) {
    for (let i = 0; i < n; i++) {
      const raw = idCol.get(i) as string | null;
      const acc = raw ? parseAccession(raw) : null;
      accForRow[i] = acc;
      if (acc) accessions.add(acc);
    }
  }

  // Stable-sorted join → deterministic across calls regardless of row order.
  const currentHash = fnv1aHex([...accessions].sort().join('|'));

  const existing = df.col(LOCATION_COL);
  if (existing && existing.getTag(LOCATION_HASH_TAG) === currentHash) {
    // Short-circuit: column already populated for this exact accession set.
    // Fire 'init-column' so the dialog's pi can advance to 100% and close
    // cleanly without sitting at 0%.
    progress?.(1, 1, 'init-column');
    return LOCATION_COL;
  }

  const locByAcc = await getSubcellularLocations([...accessions], undefined, progress);
  const vals: string[] = new Array(n);
  for (let i = 0; i < n; i++) {
    const acc = accForRow[i];
    vals[i] = (acc && locByAcc.get(acc)) || 'Unknown';
  }

  const col = existing ?? df.columns.addNewString(LOCATION_COL);
  col.init((i) => vals[i]);
  col.meta.colors.setCategorical(LOCATION_COLORS);
  col.semType = SEMTYPE.SUBCELLULAR_LOCATION;
  // Column-level cache-invalidation tag — see plan 13-08 patterns-established.
  col.setTag(LOCATION_HASH_TAG, currentHash);
  // Final progress tick after bulk-init so the dialog's pi advances past
  // the network phase into 'init-column' done.
  progress?.(1, 1, 'init-column');
  return LOCATION_COL;
}

/** Replaces (never stacks) the volcano threshold lines: horizontal at the
 * significance cutoff and the ±fold-change verticals. */
function applyThresholdLines(
  df: DG.DataFrame, yColName: string, fcThreshold: number, pThreshold: number,
): void {
  const hLine = -Math.log10(pThreshold);
  const yFormulaPrefix = `\${${yColName}}`;
  const fcFormulaPrefix = `\${log2FC}`;
  df.meta.formulaLines.items = df.meta.formulaLines.items.filter((line) => {
    const f = line.formula ?? '';
    return !(typeof f === 'string' &&
      (f.startsWith(yFormulaPrefix) || f.startsWith(fcFormulaPrefix)));
  });
  df.meta.formulaLines.addLine(
    {formula: `${yFormulaPrefix} = ${hLine}`, color: '#888888', width: 1, visible: true});
  df.meta.formulaLines.addLine(
    {formula: `${fcFormulaPrefix} = ${fcThreshold}`, color: '#888888', width: 1, visible: true});
  df.meta.formulaLines.addLine(
    {formula: `${fcFormulaPrefix} = ${-fcThreshold}`, color: '#888888', width: 1, visible: true});
}

/** UI-SPEC §"Copywriting": Y-axis label depends on the active metric;
 * defaults to adj.p-value mode when the tag is absent. */
function yAxisLabelForMetric(metric: MetricKind | null | undefined): string {
  return metric === 'p-value' ? '-Log10(p-value)' : '-Log10(Q-value)';
}

/** G2 dialog-preload contract: every input in the Volcano Options dialog is
 * seeded from {@link readVolcanoState}. Source-of-truth precedence:
 *   1. `sp.getOptions().look.colorColumnName` → colorDim
 *   2. `df.tag(VOLCANO_METRIC_TAG)` → metric (set by recomputeVolcano)
 * Pitfall 2 — the dialog snapshots at open time; OK uses input values, never
 * a re-read of sp.getOptions. */
export interface VolcanoState {
  metric: MetricKind;
  colorDim: ColorDim;
}

export function readVolcanoState(df: DG.DataFrame, sp: DG.ScatterPlotViewer): VolcanoState {
  const opts = sp.getOptions() as any;
  const look = opts?.look ?? {};
  const colorDim: ColorDim = look.colorColumnName === LOCATION_COL ? 'location' : 'significance';
  const metric = (df.getTag(VOLCANO_METRIC_TAG) as MetricKind | null) ?? 'adj.p-value';
  return {metric, colorDim};
}

/** Volcano title synthesis from `proteomics.groups` per UI-SPEC §"Copywriting":
 * `Volcano Plot: {g1} vs {g2}`, or `Volcano Plot` when groups are unset. */
function synthesizeVolcanoTitle(df: DG.DataFrame): string {
  const groups = getGroups(df);
  return groups ?
    `Volcano Plot: ${groups.group1.name} vs ${groups.group2.name}` :
    'Volcano Plot';
}

/** Disposes every subscription in `activeVolcanoSubscriptions` and removes any
 * axis-label / counter overlays attached to `sp.root`. Called at the start of
 * every createVolcanoPlot so re-entry never stacks DOM or strands subscribers
 * pointing at a detached viewer (Pitfall 6). */
function disposeVolcanoAttachments(sp?: DG.ScatterPlotViewer): void {
  for (const s of activeVolcanoSubscriptions) s.unsubscribe();
  activeVolcanoSubscriptions = [];
  if (sp) {
    for (const sel of [
      '[data-volcano-axis-x]', '[data-volcano-axis-y]', '[data-volcano-counter]',
      '[data-volcano-busy]',
    ])
      sp.root.querySelectorAll(sel).forEach((el) => el.remove());
  }
}

/** Re-reads the Y-axis label DOM (if present) from the current metric tag.
 * Called by recomputeVolcano after the tag is persisted so the label reflects
 * the new metric immediately (sp.onPropertyChanged doesn't fire when only the
 * tag changes — the yColumnName stays stable). */
function refreshVolcanoAxisLabel(sp: DG.ScatterPlotViewer, df: DG.DataFrame): void {
  const yLabel = sp.root.querySelector<HTMLElement>('[data-volcano-axis-y]');
  if (yLabel) {
    const metric = df.getTag(VOLCANO_METRIC_TAG) as MetricKind | null;
    yLabel.textContent = yAxisLabelForMetric(metric);
  }
}

/** D-06 live counter overlay: floating div anchored bottom-right of the
 * volcano. Heading + Total + per-category counts (significance categories
 * when color = direction, subcellular location categories when color =
 * Subcellular Location). Counts iterate `df.filter` only (selection toggles
 * trigger a recompute defensively but never change the displayed counts —
 * CONTEXT D-06). Subscriptions are pushed into `activeVolcanoSubscriptions`
 * so disposeVolcanoAttachments tears them down on re-entry (Pitfall 6). */
function attachCounterOverlay(sp: DG.ScatterPlotViewer, df: DG.DataFrame): void {
  const host = sp.root;
  host.style.position = host.style.position || 'relative';

  const overlay = ui.divV([]);
  overlay.dataset.volcanoCounter = 'true';
  Object.assign(overlay.style, {
    position: 'absolute', bottom: '8px', right: '8px',
    padding: '8px 16px', background: 'rgba(255,255,255,0.85)',
    borderRadius: '4px', fontSize: '0.9em',
    pointerEvents: 'none', zIndex: '5',
  } as Partial<CSSStyleDeclaration>);
  host.appendChild(overlay);

  const recompute = (): void => {
    const total = df.filter.trueCount;
    const colorColName = sp.props.colorColumnName;
    const colorCol = colorColName ? df.col(colorColName) : null;

    // Iterate categories from the column's own category list so zero-count
    // categories still render (D-06: empty filter still shows zeros, not
    // just the categories present in the filtered subset).
    const counts = new Map<string, number>();
    let perCategoryEnabled = false;
    if (colorCol && (colorColName === DIRECTION_COL || colorColName === LOCATION_COL)) {
      perCategoryEnabled = true;
      for (const cat of colorCol.categories) counts.set(cat, 0);
      for (let i = 0; i < df.rowCount; i++) {
        if (!df.filter.get(i)) continue;
        if (colorCol.isNone(i)) continue;
        const cat = colorCol.get(i) as string;
        counts.set(cat, (counts.get(cat) ?? 0) + 1);
      }
    }

    overlay.replaceChildren();
    const heading = ui.divText('Visible Proteins');
    heading.style.fontWeight = 'bold';
    overlay.appendChild(heading);
    overlay.appendChild(ui.divText(`Total: ${total.toLocaleString()}`));
    if (perCategoryEnabled) {
      for (const [cat, n] of counts)
        overlay.appendChild(ui.divText(`${cat}: ${n.toLocaleString()}`));
    }
  };

  recompute();
  activeVolcanoSubscriptions.push(
    df.onFilterChanged.pipe(debounceTime(50)).subscribe(recompute),
    df.onSelectionChanged.pipe(debounceTime(50)).subscribe(recompute),
    sp.onPropertyValueChanged.pipe(debounceTime(50)).subscribe(recompute),
  );
}

/** Attaches DOM overlays for X and Y axis labels (the scatter has no
 * xAxisCustomTitle / yAxisCustomTitle prop). X label is invariant; Y label is
 * driven from VOLCANO_METRIC_TAG and re-rendered by refreshVolcanoAxisLabel. */
function attachAxisLabels(sp: DG.ScatterPlotViewer, df: DG.DataFrame): void {
  const host = sp.root;
  host.style.position = host.style.position || 'relative';

  const xLabel = ui.divText('Log2 Fold Change');
  xLabel.dataset.volcanoAxisX = 'true';
  Object.assign(xLabel.style, {
    position: 'absolute', bottom: '4px', left: '50%',
    transform: 'translateX(-50%)', fontSize: '0.85em',
    pointerEvents: 'none', background: 'rgba(255,255,255,0.7)',
    padding: '0 4px', zIndex: '4',
  } as Partial<CSSStyleDeclaration>);
  host.appendChild(xLabel);

  const metric = df.getTag(VOLCANO_METRIC_TAG) as MetricKind | null;
  const yLabel = ui.divText(yAxisLabelForMetric(metric));
  yLabel.dataset.volcanoAxisY = 'true';
  Object.assign(yLabel.style, {
    position: 'absolute', left: '4px', top: '50%',
    transform: 'translateY(-50%) rotate(-90deg)',
    transformOrigin: 'left center', fontSize: '0.85em',
    pointerEvents: 'none', background: 'rgba(255,255,255,0.7)',
    padding: '0 4px', zIndex: '4',
  } as Partial<CSSStyleDeclaration>);
  host.appendChild(yLabel);

  // Defensive: refresh on any property change as a secondary trigger
  // (recomputeVolcano calls refreshVolcanoAxisLabel directly — this catches
  // any external prop edit, e.g. user changes color column via property panel).
  activeVolcanoSubscriptions.push(
    sp.onPropertyValueChanged.pipe(debounceTime(50))
      .subscribe(() => refreshVolcanoAxisLabel(sp, df)),
  );
}

/** 13-10 in-volcano busy overlay. Centered card on `sp.root` showing a phase
 * label and an optional detail line; updated per ProgressCb tick from the
 * caller (e.g. volcanoOptions during a first-time subcellular-location fetch).
 *
 * The TaskBarProgressIndicator at the bottom of the platform shell is too easy
 * to miss while the user stares at a stale-looking volcano (post-13-08 UAT
 * Test 3). This overlay anchors progress on the surface the user is actually
 * watching. Tagged with `data-volcano-busy` so disposeVolcanoAttachments
 * sweeps it on createVolcanoPlot re-entry. */
export function showVolcanoBusy(sp: DG.ScatterPlotViewer, label: string): void {
  // Idempotent: replace any pre-existing overlay rather than stacking.
  hideVolcanoBusy(sp);

  const host = sp.root;
  host.style.position = host.style.position || 'relative';

  const overlay = ui.divV([]);
  overlay.dataset.volcanoBusy = 'true';
  Object.assign(overlay.style, {
    position: 'absolute', top: '50%', left: '50%',
    transform: 'translate(-50%, -50%)',
    padding: '12px 20px', background: 'rgba(255,255,255,0.95)',
    border: '1px solid #ccc', borderRadius: '6px',
    boxShadow: '0 2px 8px rgba(0,0,0,0.15)',
    fontSize: '0.95em', textAlign: 'center',
    pointerEvents: 'none', zIndex: '6',
    maxWidth: '320px',
  } as Partial<CSSStyleDeclaration>);

  const labelDiv = ui.divText(label);
  labelDiv.style.fontWeight = 'bold';
  labelDiv.style.marginBottom = '4px';
  const detailDiv = ui.divText('');
  detailDiv.style.fontSize = '0.85em';
  detailDiv.style.color = '#666';
  overlay.appendChild(labelDiv);
  overlay.appendChild(detailDiv);
  host.appendChild(overlay);
}

export function updateVolcanoBusy(
  sp: DG.ScatterPlotViewer, label: string, detail?: string,
): void {
  const overlay = sp.root.querySelector<HTMLElement>('[data-volcano-busy]');
  if (!overlay) return; // already torn down — caller raced detach
  const [labelDiv, detailDiv] = Array.from(overlay.children) as HTMLElement[];
  if (labelDiv) labelDiv.textContent = label;
  if (detailDiv) detailDiv.textContent = detail ?? '';
}

export function hideVolcanoBusy(sp: DG.ScatterPlotViewer): void {
  sp.root.querySelectorAll('[data-volcano-busy]').forEach((el) => el.remove());
}

/** Column that drives the volcano's point labels. Sparse: it holds the chosen
 * display name only for the labeled rows and '' everywhere else. `~`-prefixed so
 * it stays out of the grid. Labels bind to THIS column (showLabelsFor='All')
 * instead of to df.selection, so labeling no longer hijacks the user's selection
 * — that coupling surfaced as the confusing orange "N selected" on open and let
 * the enrichment cross-link wipe the labels. */
export const VOLCANO_LABEL_COL = '~Volcano label';
/** Persisted top-N label count so the options dialog + search path can read it. */
export const VOLCANO_TOPN_TAG = 'proteomics.volcano_top_n';

/** Display Name (Plan 14-01) is canonical; Gene name is the pre-14-01 fallback. */
function findLabelSource(df: DG.DataFrame): DG.Column | null {
  return findColumn(df, SEMTYPE.DISPLAY_NAME, ['display name']) ??
    findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol']);
}

/** Ranks filtered-in rows by significance (−log10p desc, |log2FC| desc) and
 * returns the top-N row indices. Excludes FLOAT_NULL on either axis. */
function topNBySignificance(df: DG.DataFrame, n: number): number[] {
  if (n <= 0) return [];
  const yCol = df.col(NEG_LOG_COL);
  const fcCol = df.col('log2FC');
  if (!yCol || !fcCol) return [];
  const yRaw = yCol.getRawData() as Float32Array | Float64Array;
  const fcRaw = fcCol.getRawData() as Float32Array | Float64Array;

  const candidates: number[] = [];
  for (let i = 0; i < df.rowCount; i++) {
    if (!df.filter.get(i)) continue;
    if (yRaw[i] === DG.FLOAT_NULL || fcRaw[i] === DG.FLOAT_NULL) continue;
    candidates.push(i);
  }
  candidates.sort((a, b) => {
    const dy = yRaw[b] - yRaw[a]; // higher -log10(p) first
    if (dy !== 0) return dy;
    return Math.abs(fcRaw[b]) - Math.abs(fcRaw[a]);
  });
  return candidates.slice(0, n);
}

/** Reads the persisted top-N label count (default 15). */
export function getVolcanoTopN(df: DG.DataFrame): number {
  const n = parseInt(df.getTag(VOLCANO_TOPN_TAG) ?? '', 10);
  return Number.isFinite(n) && n >= 0 ? n : 15;
}

/** D-03 top-N point labels — decoupled from df.selection (Pitfall: selection-as-
 * labels conflated the orange highlight with the label set and let the
 * enrichment cross-link wipe it). Populates the sparse {@link VOLCANO_LABEL_COL}
 * with the top-N proteins by significance (plus any `extraRows`, e.g. search
 * matches), '' everywhere else, and binds it as the scatter's label source. The
 * user's selection is never touched. `n=0` with no `extraRows` clears labels. */
export function applyTopNLabels(
  df: DG.DataFrame,
  sp: DG.ScatterPlotViewer,
  n: number = 15,
  extraRows?: Iterable<number>,
): void {
  df.setTag(VOLCANO_TOPN_TAG, String(n));
  const src = findLabelSource(df);
  if (!src) return; // nothing to label with → don't add an empty technical column

  // ensureFreshFloat-style idempotency: reuse the column if present, else add.
  const labelCol = df.col(VOLCANO_LABEL_COL) ?? df.columns.addNewString(VOLCANO_LABEL_COL);

  const show = new Set<number>(topNBySignificance(df, n));
  if (extraRows) for (const i of extraRows) show.add(i);

  // Bulk init (memory feedback_dg_column_bulk_init — never per-row set).
  labelCol.init((i) => show.has(i) ? String(src.get(i) ?? '') : '');

  sp.props.labelColumnNames = [VOLCANO_LABEL_COL];
  sp.props.showLabelsFor = 'All';
}

/** Creates a volcano plot (ScatterPlotViewer). Defaults match the client
 * CK-omics figure: metric = adj.p-value, color = significance (magenta/cyan/
 * gray per D-04). G1 title is synthesized from `proteomics.groups`, D-03
 * top-N labels are seeded into `df.selection`, label binding prefers
 * `Display Name` (Plan 14-01 contract) with `Gene name` as a graceful fallback
 * for DataFrames predating Plan 14-01. */
export function createVolcanoPlot(
  df: DG.DataFrame,
  options?: {fcThreshold?: number; pThreshold?: number; topNLabels?: number; title?: string},
): DG.ScatterPlotViewer {
  const fcThreshold = options?.fcThreshold ?? 1.0;
  const pThreshold = options?.pThreshold ?? 0.05;

  // Re-entry contract — dispose every overlay / subscription from a prior
  // createVolcanoPlot call BEFORE attaching new ones (Pitfall 6). We do this
  // again AFTER the new sp is built (with sp passed in) so any orphan DOM on
  // the new sp.root from a previous render is cleared too.
  disposeVolcanoAttachments();

  const initialMetric: MetricKind = (df.getTag(VOLCANO_METRIC_TAG) as MetricKind | null) ?? 'adj.p-value';
  df.setTag(VOLCANO_METRIC_TAG, initialMetric);

  const yColName = ensureNegLog10Column(df, initialMetric);
  const colorColName = ensureDirectionColumn(df, fcThreshold, pThreshold, initialMetric);

  const sp = df.plot.scatter({
    x: 'log2FC',
    y: yColName,
    color: colorColName,
  });

  applyThresholdLines(df, yColName, fcThreshold, pThreshold);
  sp.props.showViewerFormulaLines = true;

  // G1 title synthesis — runs AFTER any caller-passed options.title so the
  // contract is the synthesized string. Callers that need an explicit title
  // can read the synthesized value back via sp.getOptions().look.title.
  const title = options?.title ?? synthesizeVolcanoTitle(df);
  sp.setOptions({title});

  attachAxisLabels(sp, df);
  attachCounterOverlay(sp, df);

  // D-03 (decoupled): label the top-N by significance via the sparse label
  // column — NOT via df.selection. topNLabels=0 opts out (no labels).
  applyTopNLabels(df, sp, options?.topNLabels ?? 15);

  return sp;
}

/** Single synchronized recompute: re-runs the metric-parameterized Y and
 * direction helpers and the threshold lines together, then points the scatter
 * at the chosen color dimension. After a Q↔P toggle the Y values, up/down/NS
 * classes and the threshold line all move together — they never disagree.
 * `colorDim`: 'significance' → direction column, 'location' → the locked
 * Subcellular-Location palette. `metric` 'p-value' is guarded (falls back to
 * adj.p-value when the column is absent). */
export async function recomputeVolcano(
  df: DG.DataFrame,
  sp: DG.ScatterPlotViewer,
  metric: MetricKind,
  colorDim: ColorDim,
  fcThreshold = 1.0,
  pThreshold = 0.05,
  progress?: ProgressCb,
): Promise<void> {
  // Persist BEFORE the column-recompute path so the axis-label refresh below
  // sees the new metric. The tag is the single source of truth for G2.
  df.setTag(VOLCANO_METRIC_TAG, metric);

  const yName = ensureNegLog10Column(df, metric);
  const dirName = ensureDirectionColumn(df, fcThreshold, pThreshold, metric);
  applyThresholdLines(df, yName, fcThreshold, pThreshold);

  sp.props.xColumnName = 'log2FC';
  sp.props.yColumnName = yName;
  if (colorDim === 'location')
    sp.props.colorColumnName = await ensureLocationColumn(df, progress);
  else
    sp.props.colorColumnName = dirName;

  refreshVolcanoAxisLabel(sp, df);
}
