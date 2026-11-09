/**
 * Composite-mean pooling for destructive / multi-occasion NONMEM TK designs.
 *
 * @see src/tests/reference-nonmem-tests.ts `poolByCell` — the validated
 *   pooling this promotes from the test harness into production.
 *
 * **Why.** `normalizeEventDrivenObservations` emits one row per OBSERVATION,
 * keyed by an occasion-scoped pseudo-subject (`<id>_Day28`). For a destructive-
 * composite TK study (`06_rat_tk28d_nonmem` — 360 animals sampled once each)
 * that means ~1 observation per pseudo-subject, so the dense per-subject pipeline
 * is degenerate (360 single-point profiles). The discriminating PKNCA card
 * (`rat_tk28d_nca_per_cell.csv`) is instead a composite-mean per
 * dose × sex × occasion CELL: pool the animals sharing a cell at each nominal
 * TAD, average the concentration, and run the EXISTING dense pipeline on the
 * resulting cell-keyed composite profiles. This mirrors the shipped pre-pooled
 * convention (`06_tk_serial_demo.csv`: `Subject = "D100_F_Day1"`).
 *
 * AD-1 is preserved — the transform still emits per-observation rows; pooling is
 * a separate, auto-detected adapter step the pane applies only when the
 * per-pseudo-subject profiles are not dense-viable. The NCA math, TK
 * accumulation and grouping are untouched: they consume the pooled frame
 * unchanged. Pure `DataFrame → DataFrame`.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import type {nca} from '@datagrok-libraries/sci-comp';

import {f64Column} from './adapter-utils';
import {fnv1a} from './study-schema';
import {compositeDoseWeight} from './unit-conversion';

import {SUBJECT_COL, TAD_COL, CONC_COL, DOSE_COL, OCCASION_COL} from './event-driven';
import {DERIVED_ROUTE_CODE_COL} from './route-normalize';
import {detectMixedCompoundGroups} from './column-autodetect';
import type {MixedCompoundGroup} from './column-autodetect';

/** The canonical columns the pooler produces / consumes by fixed name — never
 *  carried as a generic covariate (they are rebuilt per cell). */
const CORE_COLS: ReadonlySet<string> = new Set([
  SUBJECT_COL, TAD_COL, CONC_COL, DOSE_COL, DERIVED_ROUTE_CODE_COL, OCCASION_COL,
]);

/** Derived study-day column the pooler emits from {@link OCCASION_COL}. Named to
 *  match the TK aggregator's `/^study_?day$/i` repeat-dose gate (toxicokinetics
 *  AD-1) so Day-1 vs Day-N accumulation lights up on the pooled frame. */
export const STUDY_DAY_COL = 'StudyDay';

/** The cell's per-kg body weight (nca-cl-bw-normalised-units AD-5b): mean BW over the members
 *  with both a dose and a weight, adjusted so that per kg = mean dose / mean BW / AUC over those
 *  SAME members (`compositeDoseWeight`). Mapped as the pooled frame's `bodyWeight`; the per-animal
 *  BW column keeps its honest null on a non-constant cell. In the source BW column's unit. */
export const CELL_BW_COL = '~CELL_MEAN_BW';

/** Resolved group-key column names for the cell key. Both nullable so a design
 *  missing a covariate degrades gracefully (key omits the absent dimension). */
export interface CompositePoolColumns {
  /** Per-kg (nominal) dose-group column — the cell-key dose dimension AND the
   *  stable cross-cell grouping key (e.g. `DOSE_MGKG`). Null ⇒ fall back to the
   *  absolute {@link DOSE_COL} rounded (less canonical, but never silent). */
  readonly doseGroupCol: string | null;
  /** Sex column for the cell key (e.g. `Sex`). Null ⇒ the key omits sex. */
  readonly sexCol: string | null;
  /** Per-animal body-weight column; when set, the pooled frame carries {@link CELL_BW_COL}. */
  readonly bodyWeightCol?: string | null;
}

export interface CompositePoolResult {
  /** The cell-keyed composite-mean observations frame (one row per cell × TAD).
   *  Carries the canonical contract columns + `StudyDay` + every cell-constant
   *  covariate (Sex / DOSE_MGKG / NOAEL / Human* …) so the existing dense
   *  pipeline, `appendGroupingColumns` and TK aggregator consume it unchanged. */
  readonly pooled: DG.DataFrame;
  /** Distinct composite cells (= row count of the resulting dense results grid). */
  readonly cellCount: number;
  /** Per cell key: an FNV-1a hash of the cell's MEMBER observations (source subject id, TAD,
   *  concentration — canonical strings, source order). Feeds the review-decision fingerprint
   *  (nca-workflow-stages AC-S1.4.3, R1/A1): reassigning an animal between cells changes the
   *  hash even when the cell-id string does not. Built in the Pass-2 loop — no new math. */
  readonly memberHashByCell: ReadonlyMap<string, string>;
  /** Cells that pool animals of ≥ 2 compounds (a column named like a compound id varies within
   *  the cell — BUG-07 pattern). The cell key has no compound dimension, so such a composite
   *  averages different compounds; the caller must REFUSE the pooled frame when non-empty. */
  readonly mixedCompound: ReadonlyArray<MixedCompoundGroup>;
  /** Per cell key: the per-kg notes (`COMPOSITE_BW_PARTIAL`, `COMPOSITE_NOT_WEIGHT_BASED`) — the
   *  caller merges them into that cell's profile warnings. Empty without a body-weight column. */
  readonly warningsByCell: ReadonlyMap<string, ReadonlyArray<nca.ParameterWarning>>;
}

/** Per-cell aggregate scratch — accumulated over all the cell's source rows. */
interface CellAgg {
  readonly cell: string;
  doseSum: number; // Σ finite absolute dose (mean → the composite "mean animal" dose)
  doseN: number;
  readonly route: number; // first route code seen in the cell
  readonly occasion: number;
  readonly firstRow: number; // any cell row — covariates are cell-constant
}

/** Read a numeric cell as JS number, mapping platform nulls → NaN. */
function readNum(col: DG.Column, row: number): number {
  const raw = col.getRawData();
  const v = raw[row];
  if (raw instanceof Float64Array || raw instanceof Float32Array)
    return v === DGNS.FLOAT_NULL ? NaN : v;
  if (raw instanceof Int32Array)
    return v === DGNS.INT_NULL ? NaN : v;
  return typeof v === 'number' ? v : NaN;
}

/**
 * Pool a normalized event-driven frame into composite-mean profiles, one per
 * (dose × sex × occasion) cell. Emits `Subject` = the cell key
 * (`D{dose}_{sex}_Day{occasion}`), the per-cell-TAD mean concentration, the
 * cell-mean absolute dose, the cell route, a derived `StudyDay`, and every
 * cell-constant covariate carried verbatim.
 */
export function poolCompositeByCell(
  norm: DG.DataFrame,
  cols: CompositePoolColumns,
): CompositePoolResult {
  const tad = mustCol(norm, TAD_COL);
  const conc = mustCol(norm, CONC_COL);
  const dose = mustCol(norm, DOSE_COL);
  const route = mustCol(norm, DERIVED_ROUTE_CODE_COL);
  const occ = mustCol(norm, OCCASION_COL);
  const doseGroup = cols.doseGroupCol !== null ? norm.col(cols.doseGroupCol) : null;
  const sex = cols.sexCol !== null ? norm.col(cols.sexCol) : null;
  const bw = cols.bodyWeightCol != null ? norm.col(cols.bodyWeightCol) : null;

  const n = norm.rowCount;

  // Pass 1 — cell-level aggregates (mean dose, route, occasion, an anchor row).
  const cells = new Map<string, CellAgg>();
  const cellMembers = new Map<string, Array<{dose: number; bw: number}>>();
  const cellKeyOf = (i: number): string => {
    const d = doseGroup !== null ? Math.round(readNum(doseGroup, i)) : Math.round(readNum(dose, i));
    const s = sex !== null && !sex.isNone(i) ? `_${sex.getString(i)}` : '';
    return `D${d}${s}_Day${readNum(occ, i)}`;
  };
  for (let i = 0; i < n; i++) {
    const cell = cellKeyOf(i);
    let g = cells.get(cell);
    if (g === undefined) {
      g = {cell, doseSum: 0, doseN: 0, route: route.get(i) as number, occasion: readNum(occ, i), firstRow: i};
      cells.set(cell, g);
      cellMembers.set(cell, []);
    }
    const d = readNum(dose, i);
    if (Number.isFinite(d)) {g.doseSum += d; g.doseN++;}
    if (bw !== null) cellMembers.get(cell)!.push({dose: d, bw: readNum(bw, i)});
  }
  // The cell's per-kg weight + notes, from the ONE helper shared with the sparse pooler (AD-5b).
  const perKgByCell = new Map<string, ReturnType<typeof compositeDoseWeight>>();
  if (bw !== null) for (const [cell, mem] of cellMembers) perKgByCell.set(cell, compositeDoseWeight(mem));

  // Pass 2 — per (cell × TAD) mean concentration. Keyed map preserves a stable
  // emission order (insertion) so the dense pipeline's lex sort is the only
  // reordering downstream.
  const points = new Map<string, {sum: number; nReal: number; tad: number; cell: string}>();
  const memberSubj = norm.col(SUBJECT_COL);
  const members = new Map<string, string[]>();
  for (let i = 0; i < n; i++) {
    const cell = cellKeyOf(i);
    let mem = members.get(cell);
    if (mem === undefined) {mem = []; members.set(cell, mem);}
    mem.push(`${memberSubj === null || memberSubj.isNone(i) ? '' : memberSubj.getString(i)}|` +
      `${tad.isNone(i) ? '' : tad.getString(i)}|${conc.isNone(i) ? '' : conc.getString(i)}`);
    const t = readNum(tad, i);
    const key = `${cell}|${t}`;
    let p = points.get(key);
    if (p === undefined) {p = {sum: 0, nReal: 0, tad: t, cell}; points.set(key, p);}
    const c = readNum(conc, i);
    if (Number.isFinite(c)) {p.sum += c; p.nReal++;}
  }

  // Emit one output row per (cell × TAD) point.
  const m = points.size;
  const outSubject: string[] = new Array(m);
  const outTad = new Float64Array(m);
  const outConc = new Float64Array(m);
  const outDose = new Float64Array(m);
  const outRoute = new Int32Array(m);
  const outStudyDay = new Int32Array(m);
  const outBw = new Float64Array(m);
  /** Source anchor row per output row — for cell-constant covariate reads. */
  const anchorRow = new Int32Array(m);

  let r = 0;
  for (const p of points.values()) {
    const g = cells.get(p.cell)!;
    outSubject[r] = p.cell;
    outTad[r] = p.tad;
    outConc[r] = p.nReal > 0 ? p.sum / p.nReal : NaN;
    outDose[r] = g.doseN > 0 ? g.doseSum / g.doseN : NaN;
    outRoute[r] = g.route;
    outStudyDay[r] = Number.isFinite(g.occasion) ? Math.round(g.occasion) : DGNS.INT_NULL;
    outBw[r] = perKgByCell.get(p.cell)?.perKgDivisor ?? NaN;
    anchorRow[r] = g.firstRow;
    r++;
  }

  const outCols: DG.Column[] = [
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, SUBJECT_COL, outSubject),
    f64Column(TAD_COL, outTad),
    f64Column(CONC_COL, outConc),
    f64Column(DOSE_COL, outDose),
    DGNS.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, outRoute),
    DGNS.Column.fromInt32Array(STUDY_DAY_COL, outStudyDay),
    ...(bw !== null ? [f64Column(CELL_BW_COL, outBw)] : []),
  ];

  // Carry every covariate (Sex / DOSE_MGKG / ROUTE string / NOAEL / Human* / LLOQ …)
  // at each output row's anchor row — BUT only where the value is constant within
  // the cell. A destructive composite cell pools distinct animals, so a per-animal
  // column (e.g. body weight) has no single cell value: emitting an arbitrary
  // anchor-animal value would be silently wrong (a composite "mean animal" has no
  // real BW). Such columns are emitted as null per non-constant cell — the honest
  // representation. Cell-constant study/cell metadata (Sex, dose level, NOAEL,
  // human-projection constants) is carried verbatim. (Independent of the
  // downstream `appendGroupingColumns` constancy gate, which only protects
  // group-key SELECTION; this protects the column VALUES for any direct reader.)
  const carriedCols = [...(norm.columns as Iterable<DG.Column>)]
    .filter((c) => !CORE_COLS.has(c.name) && c.name !== STUDY_DAY_COL);
  const nonConstantCells = computeNonConstantCells(carriedCols, cellKeyOf, n);
  for (const src of carriedCols)
    outCols.push(buildCarriedColumn(src, anchorRow, outSubject, nonConstantCells.get(src.name)!, m));

  const memberHashByCell = new Map<string, string>();
  for (const [cell, mem] of members) memberHashByCell.set(cell, fnv1a(mem.join('\n')));
  // Checked with the SAME cell key the pooling used (one source for "which rows share a cell").
  const mixedCompound = detectMixedCompoundGroups(norm, cellKeyOf);
  const warningsByCell = new Map<string, ReadonlyArray<nca.ParameterWarning>>();
  for (const [cell, w] of perKgByCell) if (w.warnings.length > 0) warningsByCell.set(cell, w.warnings);
  return {pooled: DGNS.DataFrame.fromColumns(outCols), cellCount: cells.size, memberHashByCell, mixedCompound,
    warningsByCell};
}

function mustCol(df: DG.DataFrame, name: string): DG.Column {
  const col = df.col(name);
  if (col === null)
    throw new Error(`poolCompositeByCell: required column '${name}' missing from the normalized frame.`);
  return col;
}

/** For each carried column, the set of cell keys where the value is NOT constant
 *  across the cell's source rows (single pass, string-normalized comparison). A
 *  cell in the set gets a null carried value (it has no single representative). */
function computeNonConstantCells(
  carriedCols: ReadonlyArray<DG.Column>, cellKeyOf: (i: number) => string, n: number,
): Map<string, Set<string>> {
  const nonConstant = new Map<string, Set<string>>();
  const firstVal = new Map<string, Map<string, string | null>>();
  for (const c of carriedCols) {nonConstant.set(c.name, new Set()); firstVal.set(c.name, new Map());}
  for (let i = 0; i < n; i++) {
    const cell = cellKeyOf(i);
    for (const c of carriedCols) {
      const v = c.isNone(i) ? null : c.getString(i);
      const fm = firstVal.get(c.name)!;
      if (!fm.has(cell)) fm.set(cell, v);
      else if (fm.get(cell) !== v) nonConstant.get(c.name)!.add(cell);
    }
  }
  return nonConstant;
}

/** Build a carried covariate column by reading `src` at each output row's anchor
 *  row, preserving the source column type. Emits null where the column is not
 *  constant within the output row's cell (`nonConstant` set). */
function buildCarriedColumn(
  src: DG.Column, anchorRow: Int32Array, outSubject: ReadonlyArray<string>,
  nonConstant: ReadonlySet<string>, m: number,
): DG.Column {
  const isNull = (i: number): boolean => nonConstant.has(outSubject[i]) || src.isNone(anchorRow[i]);
  switch (src.type) {
  case DGNS.COLUMN_TYPE.INT: {
    const out = new Int32Array(m);
    for (let i = 0; i < m; i++) out[i] = isNull(i) ? DGNS.INT_NULL : (src.get(anchorRow[i]) as number);
    return DGNS.Column.fromInt32Array(src.name, out);
  }
  case DGNS.COLUMN_TYPE.FLOAT: {
    const out = new Float64Array(m);
    for (let i = 0; i < m; i++) out[i] = isNull(i) ? DGNS.FLOAT_NULL : (src.get(anchorRow[i]) as number);
    return f64Column(src.name, out);
  }
  case DGNS.COLUMN_TYPE.BOOL: {
    const out = new Array<boolean | null>(m);
    for (let i = 0; i < m; i++) out[i] = isNull(i) ? null : (src.get(anchorRow[i]) as boolean);
    return DGNS.Column.fromList(DGNS.COLUMN_TYPE.BOOL, src.name, out);
  }
  default: {
    const out = new Array<string | null>(m);
    for (let i = 0; i < m; i++) out[i] = isNull(i) ? null : src.getString(anchorRow[i]);
    return DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, src.name, out);
  }
  }
}
