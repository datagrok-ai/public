/**
 * Results writer — turn an array of `nca.ComputeResult` (one per subject) into
 * a Datagrok DataFrame following the 36-column UPPERCASE_UNDERSCORE schema
 * (22 base + `C0_EXTRAP` + per-row `CL_BASIS` + per-row `DOSE_NORM` + the 6
 * FR-200 moment / lag params `AUMC_LAST, AUMC_INF, MRT, VSS, TLAG,
 * PCT_AUMC_EXTRAP`, sci-comp 0.7.0 + the λz `SPAN_RATIO` diagnostic,
 * sci-comp 0.10.0 + the absolute companions `CL_ABS, VZ_ABS, VSS_ABS` and
 * `BW_KG`).
 *
 * **Float64 invariant.** Every numeric output column is created via
 * `adapter-utils.f64Column(name, array)` — which wraps
 * `DG.Column.fromFloat64Array`. Never `Column.float(...)`, never
 * `addNewFloat(...)`: both alternatives yield Float32 storage and break PKNCA
 * bit-identity (NFR-04).
 *
 * **Missing-value sentinel — `DG.FLOAT_NULL`, written at construction**
 * (nca-dg-native-missing-values, user decision D1, 2026-08-03). This REVERSES
 * the former NaN convention that this header used to document, and the reversal
 * was anticipated here ("if it doesn't, this schema reverts to `DG.FLOAT_NULL`").
 *
 * sci-comp returns IEEE-754 NaN for an unfit value (e.g. `ParameterValues.aucInf`
 * when λz didn't fit) and that contract is immutable (rule 5, BFIELD-01..15). But
 * a raw NaN copied verbatim into a DG column is NOT `none` to the platform:
 * `isNone` returns false, so grid sort, column filter, `col.stats`, group-by and
 * CSV/`.d42` export all read an unfit parameter as a PRESENT value — and one such
 * cell collapses min/max/avg/sum/stdev on the whole column. `f64Column` performs
 * the NaN → `FLOAT_NULL` translation at the DataFrame boundary, so the conversion
 * is part of construction rather than a post-pass someone can forget.
 *
 * The `.d42` question this header used to defer is CLOSED and favourable:
 * `FLOAT_NULL` survives `toByteArray`/`fromByteArray` BIT-EXACT and still reads
 * `isNone` (probe PR-6, pinned in `adapter-results-writer-tests.ts`). No schema
 * change, and no parallel `*_VALID` bool-mask columns are needed.
 *
 * `LAMBDAZ_NPOINTS` is Int32 and uses `DG.INT_NULL`, spelled literally below —
 * it needs no helper because it is the only Int32 column with a missing case.
 *
 * **Per-row CL / VZ basis (per-row-route AD-4).** A single column holds the
 * numeric `cl` for every subject — the engine doesn't know about route. Its
 * physical interpretation (`CL` for intravascular vs `CL/F` for extravascular,
 * and `Vz` vs `Vz/F`) is **per subject**, following that subject's route, and
 * is written to a dedicated per-row `CL_BASIS` STRING column (`'CL'` | `'CL/F'`;
 * the same basis governs Vz). This replaces the former all-or-nothing
 * `allIntravenous` study flag, which mislabeled every subject of a mixed-route
 * frame. The `CL` / `VZ` numeric columns therefore carry the **neutral**
 * `tags.label` (`'CL'` / `'Vz'`); read `CL_BASIS` per row for the F-corrected
 * interpretation. Single-route frames have one `CL_BASIS` value on every row
 * (unchanged appearance).
 *
 * **Per-row engine dose basis (nca-cl-vz-convention AD-3).** A second per-row
 * STRING column, `DOSE_NORM` (`'absolute'` | `'per_kg'`), records how the ENGINE
 * dose was formed: absolute mg, or per kg (a `/kg` dose with no body weight). It
 * is read by the summary-cache key, dose-proportionality and the TK guards — it
 * is NOT the reporting basis of CL / Vz / Vss.
 *
 * **Reporting basis = the column identity (nca-cl-bw-normalised-units AD-1).**
 * `CL` / `VZ` / `VSS` are per kg on EVERY row (mL/min/kg, L/kg, L/kg), formed by
 * `normaliseClVolumes` from the engine value, the engine dose basis and the
 * profile's body weight; NaN (+ `PER_KG_NOT_FORMED`) when there is no weight to
 * normalise an absolute dose by. The absolute values ride `CL_ABS` / `VZ_ABS` /
 * `VSS_ABS` (L/h, L, L) and the weight `BW_KG`. One meaning per column, so every
 * column's `units` tag comes from `CANONICAL_UNIT_TAGS` — no frame-uniformity rule.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import {f64Column} from './adapter-utils';
import {nca} from '@datagrok-libraries/sci-comp';

import type {ProfileIndex} from './profile-index';
import {CANONICAL_UNIT_TAGS, NO_BODY_WEIGHT, normaliseClVolumes} from './unit-conversion';
import type {BodyWeightContext} from './unit-conversion';

/** Inverse of `ROUTE_STRING_TO_CODE` in profile-inputs.ts, keyed by RouteCode
 *  index 0..5. This is the CANONICAL route label the results `ROUTE` column carries
 *  (and the study grid / tier show); exported so any other consumer that must display
 *  the SAME route labels (e.g. the Mean±SD chart grouping by route) uses one source. */
export const ROUTE_CODE_TO_STRING: ReadonlyArray<string> = [
  'IV-bolus', // 0  nca.ROUTE_IV_BOLUS
  'IV-infusion', // 1  nca.ROUTE_IV_INFUSION
  'PO', // 2  nca.ROUTE_PO
  'SC', // 3  nca.ROUTE_SC
  'IM', // 4  nca.ROUTE_IM
  'other', // 5  nca.ROUTE_OTHER
];

/**
 * Write the 36-column results DataFrame.
 *
 * Pre-allocates a `Float64Array` for each numeric column and an `Int32Array`
 * for `LAMBDAZ_NPOINTS`, then makes a single pass over the parallel inputs
 * arrays. Variable-length provenance (lambdaZ.pointsUsed, blqApplied.excluded)
 * is serialized into `PROVENANCE_JSON`; the rest of `provenance.warnings` is
 * serialized into `WARNINGS_JSON`.
 *
 * Invariant: `keys.length === routes.length === results.length`. Caller (the
 * Phase 4 controller) is responsible for keeping the three arrays parallel.
 */
/** The FIXED per-subject results-frame column names this writer emits (structural +
 *  every NCA parameter) — the single source of truth for "reserved results-frame column".
 *  Consumers (e.g. the v10 passthrough candidate gate `RESERVED_RESULT_COLS`) derive their
 *  collision set from THIS constant instead of re-typing the list, so a new result column can
 *  never silently become an offerable passthrough name (architect finding, 2026-07-14). The
 *  `cols` array in {@link writeResultsToDataFrame} asserts its own names match this list at
 *  build time (drift guard) — add a new column here AND there together, or the assertion trips. */
export const FIXED_RESULT_COLUMNS: readonly string[] = [
  'SUBJECT', 'ROUTE', 'CL_BASIS', 'DOSE_NORM', 'STATUS', 'CMAX', 'C0_EXTRAP', 'TMAX',
  'AUC_LAST', 'AUC_INF', 'PCT_EXTRAP', 'LAMBDA_Z', 'HALF_LIFE', 'CL', 'VZ',
  'AUMC_LAST', 'AUMC_INF', 'MRT', 'VSS', 'TLAG', 'PCT_AUMC_EXTRAP',
  'CL_ABS', 'VZ_ABS', 'VSS_ABS', 'BW_KG',
  'LAMBDAZ_R2', 'LAMBDAZ_ADJ_R2', 'LAMBDAZ_NPOINTS', 'SPAN_RATIO',
  'LAMBDAZ_TFIRST', 'LAMBDAZ_TLAST',
  'LAMBDAZ_INTERCEPT', 'AUC_METHOD', 'COMPENSATED', 'WARNINGS_JSON', 'PROVENANCE_JSON',
];

/**
 * Every per-subject scalar the results frame carries for one row. Computed once
 * by {@link computeResultRowValues} and consumed by BOTH the full-build writer
 * loop and the in-place {@link patchResultRow} — so the two paths cannot drift
 * (AC-1.6). `lzNPoints === DG.INT_NULL` and any `NaN` numeric mark a missing cell;
 * the consumers apply the DG sentinels.
 */
export interface ResultRowValues {
  readonly subject: string;
  readonly route: string;
  readonly clBasis: string;
  readonly doseNorm: string;
  readonly status: string;
  readonly aucMethod: string;
  readonly compensated: boolean;
  readonly warningsJson: string;
  readonly provenanceJson: string;
  readonly cmax: number; readonly c0: number; readonly tmax: number;
  readonly aucLast: number; readonly aucInf: number; readonly pctExtrap: number;
  // `cl` / `vz` / `vss` are per kg (mL/min/kg, L/kg, L/kg) — `normaliseClVolumes`.
  readonly lambdaZ: number; readonly halfLife: number; readonly cl: number; readonly vz: number;
  readonly aumcLast: number; readonly aumcInf: number; readonly mrt: number;
  readonly vss: number; readonly tlag: number; readonly pctAumcExtrap: number;
  /** Absolute companions (L/h, L, L) and the body weight (kg) the per-kg values used. */
  readonly clAbs: number; readonly vzAbs: number; readonly vssAbs: number; readonly bwKg: number;
  readonly lzR2: number; readonly lzAdjR2: number; readonly lzTFirst: number;
  readonly lzTLast: number; readonly lzIntercept: number; readonly spanRatio: number;
  /** Point count, or `DG.INT_NULL` when λz did not fit. */
  readonly lzNPoints: number;
}

/**
 * Compute one subject's full result-row values from the engine output. The
 * SINGLE per-row deriver shared by {@link writeResultsToDataFrame} and
 * {@link patchResultRow}. NaN numerics pass through verbatim (the DG-boundary
 * sentinel translation is the consumer's job — `f64Column` at build, per-cell at
 * patch).
 */
export function computeResultRowValues(
  k: ProfileIndex, r: nca.RouteCode, res: nca.ComputeResult, c0: number,
  doseNorm: 'absolute' | 'per_kg', bw: BodyWeightContext = NO_BODY_WEIGHT,
): ResultRowValues {
  const prov = res.provenance;
  const lz = prov.lambdaZ;
  // The per-kg rule's warnings are merged into the profile's provenance by the controller (one
  // stream for every reader); here only the values are taken.
  const norm = normaliseClVolumes(res.values, doseNorm, bw);
  return {
    subject: k.subjectKey,
    route: ROUTE_CODE_TO_STRING[r],
    clBasis: (r === nca.ROUTE_IV_BOLUS || r === nca.ROUTE_IV_INFUSION) ? 'CL' : 'CL/F',
    doseNorm,
    status: res.status,
    aucMethod: prov.aucMethod,
    compensated: prov.compensated,
    warningsJson: JSON.stringify(prov.warnings),
    provenanceJson: JSON.stringify({
      lambdaZ: lz === null ? null : {pointsUsed: Array.from(lz.pointsUsed)},
      blqApplied: {excluded: Array.from(prov.blqApplied.excluded)},
      // c0 provenance (sci-comp ≥ 0.11.0, BUG-12 residual): WHICH chain method
      // produced the dose-time concentration, whether a BLQ / non-positive
      // pre-dose row was replaced, and the back-extrapolated share of AUC₀–∞.
      // `null` for non-IV-bolus routes and on a failed profile; absent (older
      // engine / test doubles) → null. Additive; not hashed by any cache.
      c0: prov.c0 == null ? null : {
        method: prov.c0.method,
        replacedDoseTimeRow: prov.c0.replacedDoseTimeRow,
        pctAucBackExtrap: Number.isNaN(prov.c0.pctAucBackExtrap) ? null : prov.c0.pctAucBackExtrap,
      },
    }),
    cmax: res.values.cmax, c0, tmax: res.values.tmax,
    aucLast: res.values.aucLast, aucInf: res.values.aucInf, pctExtrap: res.values.pctExtrap,
    lambdaZ: res.values.lambdaZ, halfLife: res.values.halfLife,
    cl: norm.cl, vz: norm.vz,
    aumcLast: res.values.aumcLast, aumcInf: res.values.aumcInf, mrt: res.values.mrt,
    vss: norm.vss, tlag: res.values.tlag, pctAumcExtrap: res.values.pctExtrapAumc,
    clAbs: norm.clAbs, vzAbs: norm.vzAbs, vssAbs: norm.vssAbs, bwKg: norm.bwKg,
    lzR2: lz === null ? NaN : lz.rSquared,
    lzAdjR2: lz === null ? NaN : lz.adjRSquared,
    lzTFirst: lz === null ? NaN : lz.tStart,
    lzTLast: lz === null ? NaN : lz.tEnd,
    lzIntercept: lz === null ? NaN : lz.intercept,
    spanRatio: lz === null ? NaN : lz.spanRatio,
    lzNPoints: lz === null ? DGNS.INT_NULL : lz.pointsUsed.length,
  };
}

export function writeResultsToDataFrame(
  keys: ReadonlyArray<ProfileIndex>,
  routes: ReadonlyArray<nca.RouteCode>,
  results: ReadonlyArray<nca.ComputeResult>,
  c0s?: Float64Array,
  doseNorms?: ReadonlyArray<'absolute' | 'per_kg'>,
  bodyWeights?: ReadonlyArray<BodyWeightContext>,
): DG.DataFrame {
  const n = keys.length;
  if (routes.length !== n || results.length !== n) {
    throw new Error(
      `writeResultsToDataFrame: parallel-array length mismatch — keys=${n}, ` +
      `routes=${routes.length}, results=${results.length}.`,
    );
  }
  if (c0s !== undefined && c0s.length !== n) {
    throw new Error(
      `writeResultsToDataFrame: c0s length mismatch — keys=${n}, c0s=${c0s.length}.`,
    );
  }
  if (doseNorms !== undefined && doseNorms.length !== n) {
    throw new Error(
      `writeResultsToDataFrame: doseNorms length mismatch — keys=${n}, ` +
      `doseNorms=${doseNorms.length}.`,
    );
  }
  if (bodyWeights !== undefined && bodyWeights.length !== n) {
    throw new Error(
      `writeResultsToDataFrame: bodyWeights length mismatch — keys=${n}, ` +
      `bodyWeights=${bodyWeights.length}.`,
    );
  }
  // Default to all-absolute when the caller didn't supply a basis (AD-3):
  // pre-AD-2 call sites and single-compound absolute frames stay unchanged
  // (DOSE_NORM = 'absolute' on every row, CL/VZ tags = L/h / L).
  const doseNormEffective: ReadonlyArray<'absolute' | 'per_kg'> =
    doseNorms ?? new Array<'absolute' | 'per_kg'>(n).fill('absolute');
  // Default to NaN-filled when the caller didn't supply c0s: non-IV-bolus
  // profiles (and IV-bolus profiles where estimateC0 declined to return a
  // value) are reported as NaN per the writer's no-fit sentinel convention.
  const c0Effective = c0s ?? new Float64Array(n).fill(NaN);


  const subject = new Array<string>(n);
  const route = new Array<string>(n);
  const clBasis = new Array<string>(n);
  const doseNorm = new Array<string>(n);
  const status = new Array<string>(n);
  const aucMethod = new Array<string>(n);
  const compensated = new Array<boolean>(n);
  const warningsJson = new Array<string>(n);
  const provenanceJson = new Array<string>(n);

  const cmax = new Float64Array(n);
  const tmax = new Float64Array(n);
  const aucLast = new Float64Array(n);
  const aucInf = new Float64Array(n);
  const pctExtrap = new Float64Array(n);
  const lambdaZ = new Float64Array(n);
  const halfLife = new Float64Array(n);
  const cl = new Float64Array(n);
  const vz = new Float64Array(n);
  const aumcLast = new Float64Array(n);
  const aumcInf = new Float64Array(n);
  const mrt = new Float64Array(n);
  const vss = new Float64Array(n);
  const tlag = new Float64Array(n);
  const pctAumcExtrap = new Float64Array(n);
  const clAbs = new Float64Array(n);
  const vzAbs = new Float64Array(n);
  const vssAbs = new Float64Array(n);
  const bwKg = new Float64Array(n);
  const lzR2 = new Float64Array(n);
  const lzAdjR2 = new Float64Array(n);
  const lzTFirst = new Float64Array(n);
  const lzTLast = new Float64Array(n);
  const lzIntercept = new Float64Array(n);
  const spanRatio = new Float64Array(n);
  const lzNPoints = new Int32Array(n);

  for (let i = 0; i < n; i++) {
    // Single per-row deriver, shared with patchResultRow (AC-1.6) so the two
    // paths cannot drift. C0_EXTRAP is emitted from `c0Effective` directly below,
    // so `rv.c0` is unused here.
    const rv = computeResultRowValues(keys[i], routes[i], results[i], c0Effective[i], doseNormEffective[i],
      bodyWeights?.[i] ?? NO_BODY_WEIGHT);
    subject[i] = rv.subject;
    route[i] = rv.route;
    clBasis[i] = rv.clBasis;
    doseNorm[i] = rv.doseNorm;
    status[i] = rv.status;
    aucMethod[i] = rv.aucMethod;
    compensated[i] = rv.compensated;
    warningsJson[i] = rv.warningsJson;
    provenanceJson[i] = rv.provenanceJson;
    cmax[i] = rv.cmax;
    tmax[i] = rv.tmax;
    aucLast[i] = rv.aucLast;
    aucInf[i] = rv.aucInf;
    pctExtrap[i] = rv.pctExtrap;
    lambdaZ[i] = rv.lambdaZ;
    halfLife[i] = rv.halfLife;
    cl[i] = rv.cl;
    vz[i] = rv.vz;
    aumcLast[i] = rv.aumcLast;
    aumcInf[i] = rv.aumcInf;
    mrt[i] = rv.mrt;
    vss[i] = rv.vss;
    tlag[i] = rv.tlag;
    pctAumcExtrap[i] = rv.pctAumcExtrap;
    clAbs[i] = rv.clAbs;
    vzAbs[i] = rv.vzAbs;
    vssAbs[i] = rv.vssAbs;
    bwKg[i] = rv.bwKg;
    lzR2[i] = rv.lzR2;
    lzAdjR2[i] = rv.lzAdjR2;
    lzTFirst[i] = rv.lzTFirst;
    lzTLast[i] = rv.lzTLast;
    lzIntercept[i] = rv.lzIntercept;
    spanRatio[i] = rv.spanRatio;
    lzNPoints[i] = rv.lzNPoints;
  }

  // AD-4: neutral parameter-name labels; the F-corrected interpretation is the
  // per-row CL_BASIS column, not a frame-level tag.
  const clLabel = 'CL';
  const vzLabel = 'Vz';

  const cols: DG.Column[] = [
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'SUBJECT', subject),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'ROUTE', route),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'CL_BASIS', clBasis),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'DOSE_NORM', doseNorm),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'STATUS', status),
    f64Column('CMAX', cmax),
    // C0_EXTRAP — back-extrapolated conc at dose time for IV bolus (sci-comp
    // `nca.estimateC0`, faithful PKNCA `pk.calc.c0` port).
    // NaN for non-IV-bolus profiles ONLY. It is NOT NaN for a "failed" IV-bolus
    // profile: the default chain ('c0','logslope','c1','cmin','set0') ends in
    // `set0`, which returns 0 unconditionally, so `estimateC0` never returns
    // null here — an all-BLQ IV-bolus profile reports C0 = 0, not blank
    // (MEASURED against sci-comp c0.ts, 2026-07-29; sci-comp c0.test.ts:135,141
    // pins both halves). That failed-row 0 is WHY the controller's separate
    // `estimateC0` call stays the source rather than `provenance.c0.value`
    // (the engine reports `c0: null` on a failed profile — no chain ran); on
    // every non-failed IV-bolus row the two are asserted EXACTLY equal
    // (reference-suite-tests, drift guard). The method that produced the value
    // IS recorded per row since sci-comp 0.11.0 — `PROVENANCE_JSON.c0.method`.
    // Sibling to CMAX (which stays observed-only) — Ed validation obs #4.
    f64Column('C0_EXTRAP', c0Effective),
    f64Column('TMAX', tmax),
    f64Column('AUC_LAST', aucLast),
    f64Column('AUC_INF', aucInf),
    f64Column('PCT_EXTRAP', pctExtrap),
    f64Column('LAMBDA_Z', lambdaZ),
    f64Column('HALF_LIFE', halfLife),
    f64Column('CL', cl).setTag('label', clLabel),
    f64Column('VZ', vz).setTag('label', vzLabel),
    // FR-200 moment / lag params. AUMC_LAST/AUMC_INF mirror the AUC_LAST/AUC_INF
    // pair; VSS (like VZ) carries the neutral label.
    f64Column('AUMC_LAST', aumcLast),
    f64Column('AUMC_INF', aumcInf),
    f64Column('MRT', mrt),
    f64Column('VSS', vss).setTag('label', 'Vss'),
    f64Column('TLAG', tlag),
    f64Column('PCT_AUMC_EXTRAP', pctAumcExtrap),
    // Absolute companions + the weight the per-kg values used (AD-3). Neutral labels like
    // their per-kg siblings: CL_BASIS says CL vs CL/F for both.
    f64Column('CL_ABS', clAbs).setTag('label', 'CL abs.'),
    f64Column('VZ_ABS', vzAbs).setTag('label', 'Vz abs.'),
    f64Column('VSS_ABS', vssAbs).setTag('label', 'Vss abs.'),
    f64Column('BW_KG', bwKg),
    f64Column('LAMBDAZ_R2', lzR2),
    f64Column('LAMBDAZ_ADJ_R2', lzAdjR2),
    DGNS.Column.fromInt32Array('LAMBDAZ_NPOINTS', lzNPoints),
    // SPAN_RATIO — the λz window's span in half-lives, sci-comp 0.10.0. Sits
    // beside its N-points sibling: both answer "is there enough terminal phase
    // here to trust λz?", and a scientist reads them together. `f64Column` turns
    // the no-fit NaN into `DG.FLOAT_NULL` at construction, so an unfit row is
    // `isNone` rather than a finite sentinel.
    f64Column('SPAN_RATIO', spanRatio),
    f64Column('LAMBDAZ_TFIRST', lzTFirst),
    f64Column('LAMBDAZ_TLAST', lzTLast),
    f64Column('LAMBDAZ_INTERCEPT', lzIntercept),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'AUC_METHOD', aucMethod),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.BOOL, 'COMPENSATED', compensated),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'WARNINGS_JSON', warningsJson),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'PROVENANCE_JSON', provenanceJson),
  ];

  // Apply canonical `column.tags.units` per AD-11 / E11. Single source of
  // truth in `unit-conversion.ts:CANONICAL_UNIT_TAGS` — `results-writer.ts`
  // walks the map by column name; string/boolean/count/JSON columns are
  // intentionally absent from the map and keep `tags.units === undefined`.
  // CL / VZ / VSS are one basis on every row (per kg), so they are tagged from
  // the map like every other column (nca-cl-bw-normalised-units AD-1).
  for (const col of cols) {
    const tag = CANONICAL_UNIT_TAGS[col.name];
    if (tag !== undefined) col.setTag('units', tag);
  }

  // Drift guard: the exported FIXED_RESULT_COLUMNS (consumed by the passthrough candidate gate)
  // must stay in lockstep with the columns actually emitted here. A mismatch means a result
  // column was added/removed without updating the constant — fail loudly rather than let a stale
  // reserved-name set silently offer a colliding passthrough column.
  const actualNames = cols.map((c) => c.name);
  if (actualNames.length !== FIXED_RESULT_COLUMNS.length ||
      actualNames.some((n, i) => n !== FIXED_RESULT_COLUMNS[i])) {
    throw new Error(
      `writeResultsToDataFrame: FIXED_RESULT_COLUMNS out of sync with emitted columns — ` +
      `expected [${FIXED_RESULT_COLUMNS.join(', ')}], got [${actualNames.join(', ')}]`);
  }

  const df = DGNS.DataFrame.fromColumns(cols);
  df.name = 'NCAResults';
  return df;
}

/**
 * Patch a SINGLE subject's row in an EXISTING results frame, in place — for the
 * interactive λz single-subject recompute (`NcaController.recomputeProfile`,
 * AC-1.2/AC-1.6). Honors the interface §9.4 contract: does NOT create a new
 * DataFrame, so every viewer/grid/selection binding on the frame is preserved;
 * the patched row equals what a full {@link writeResultsToDataFrame} rebuild would
 * produce for that subject.
 *
 * The row is located by `SUBJECT === index.subjectKey`. NaN numerics are written
 * as `DG.FLOAT_NULL` PER CELL (N-2): `f64Column` applies the sentinel only at
 * array construction, so a `failed`/`partial` recompute would otherwise write a
 * raw NaN and reintroduce the missing-value defect. Unit tags are fixed per
 * column (CANONICAL_UNIT_TAGS) and are intentionally left untouched.
 */
export function patchResultRow(
  resultsDf: DG.DataFrame,
  index: ProfileIndex,
  route: nca.RouteCode,
  result: nca.ComputeResult,
  c0: number,
  doseNorm: 'absolute' | 'per_kg',
  bw: BodyWeightContext = NO_BODY_WEIGHT,
): void {
  const subjectCol = resultsDf.col('SUBJECT');
  if (subjectCol === null)
    throw new Error('patchResultRow: results frame has no \'SUBJECT\' column.');
  let row = -1;
  for (let i = 0; i < subjectCol.length; i++)
    if (subjectCol.getString(i) === index.subjectKey) {row = i; break;}

  if (row < 0)
    throw new Error(`patchResultRow: subject '${index.subjectKey}' not found in results frame.`);

  const rv = computeResultRowValues(index, route, result, c0, doseNorm, bw);

  const setStr = (name: string, v: string): void => resultsDf.col(name)?.set(row, v);
  const setNum = (name: string, v: number): void => {
    // Per-cell NaN → DG.FLOAT_NULL (N-2). ±Infinity passes through, mirroring
    // f64Column's AC-3 (only NaN is rewritten).
    resultsDf.col(name)?.set(row, Number.isNaN(v) ? DGNS.FLOAT_NULL : v);
  };

  setStr('SUBJECT', rv.subject);
  setStr('ROUTE', rv.route);
  setStr('CL_BASIS', rv.clBasis);
  setStr('DOSE_NORM', rv.doseNorm);
  setStr('STATUS', rv.status);
  setStr('AUC_METHOD', rv.aucMethod);
  resultsDf.col('COMPENSATED')?.set(row, rv.compensated);
  setStr('WARNINGS_JSON', rv.warningsJson);
  setStr('PROVENANCE_JSON', rv.provenanceJson);

  setNum('CMAX', rv.cmax);
  setNum('C0_EXTRAP', rv.c0);
  setNum('TMAX', rv.tmax);
  setNum('AUC_LAST', rv.aucLast);
  setNum('AUC_INF', rv.aucInf);
  setNum('PCT_EXTRAP', rv.pctExtrap);
  setNum('LAMBDA_Z', rv.lambdaZ);
  setNum('HALF_LIFE', rv.halfLife);
  setNum('CL', rv.cl);
  setNum('VZ', rv.vz);
  setNum('AUMC_LAST', rv.aumcLast);
  setNum('AUMC_INF', rv.aumcInf);
  setNum('MRT', rv.mrt);
  setNum('VSS', rv.vss);
  setNum('TLAG', rv.tlag);
  setNum('PCT_AUMC_EXTRAP', rv.pctAumcExtrap);
  setNum('CL_ABS', rv.clAbs);
  setNum('VZ_ABS', rv.vzAbs);
  setNum('VSS_ABS', rv.vssAbs);
  setNum('BW_KG', rv.bwKg);
  setNum('LAMBDAZ_R2', rv.lzR2);
  setNum('LAMBDAZ_ADJ_R2', rv.lzAdjR2);
  setNum('SPAN_RATIO', rv.spanRatio);
  setNum('LAMBDAZ_TFIRST', rv.lzTFirst);
  setNum('LAMBDAZ_TLAST', rv.lzTLast);
  setNum('LAMBDAZ_INTERCEPT', rv.lzIntercept);
  // LAMBDAZ_NPOINTS is Int32 — rv.lzNPoints is already DG.INT_NULL when λz did
  // not fit (computeResultRowValues), so set it verbatim.
  resultsDf.col('LAMBDAZ_NPOINTS')?.set(row, rv.lzNPoints);
}
