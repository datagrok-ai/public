/**
 * Results aggregator — SPIKE (topic: per-compound-aggregation).
 *
 * Group-by aggregation of the per-subject NCA results DataFrame (the
 * `writeResultsToDataFrame` output) into a per-group summary table: arithmetic
 * mean, sample SD, CV%, and N per NCA parameter, grouped by a chosen column.
 *
 * **Placement (CLAUDE.md rule 5 / Architecture Gotchas).** This is pure
 * post-`computeNca` post-processing on a DG.DataFrame — NOT a numerical
 * primitive. It uses the native `DataFrame.groupBy().aggregate()` builder and
 * adds no math sci-comp already owns. No re-derivation, no upstream concern.
 *
 * **Arithmetic mean default (VAL-01 parity).** Ed's reference
 * `multicompound_nca_per_compound.csv` is computed with R `mean()` (arithmetic);
 * matching it keeps the rule-18 PKNCA-parity gate green. Geometric mean +
 * geo-CV% (the DMPK exposure-reporting convention) is an additive follow-up,
 * NOT wired here — see the gap note appended to TODO.md by this spike.
 *
 * **SPIKE SCOPE.** Assumes `results` already carries `groupColName`. Production
 * must propagate a grouping column into the results table — the per-subject
 * value is recoverable via `perProfile[i].index.observationRows[0]` indexing
 * back into the source DataFrame (the grouping metadata, e.g. CompoundID, is
 * constant per subject). Wiring that (controller / results-writer carry-through
 * or a join helper) + the UI "Group by" control + docking the grid are
 * deliberately out of spike scope.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import {compositeGroupKey} from './composite-group-key';
import {f64Column} from './adapter-utils';

/**
 * Numeric NCA parameter columns aggregated by default. Mirrors the
 * `writeResultsToDataFrame` numeric schema; the headline DMPK screening set
 * (Cmax / AUC / t½ / CL / Vz / %extrap) plus Tmax and λz for completeness, and
 * the FR-200 moment headlines (AUMC_INF / MRT / Vss).
 *
 * `Vss` pools IV rows only for free: a non-IV row's Vss is the platform null —
 * written that way by the results writer at construction (`f64Column`), where it
 * used to be converted by a shadow clone in this module — so DG `avg` skips it.
 * `CL` / `VZ` / `VSS` are per kg on every row (nca-cl-bw-normalised-units AD-1),
 * so a group mixing `absolute` and `per_kg` engine doses pools ONE basis; the
 * absolute companions `CL_ABS` / `VZ_ABS` / `VSS_ABS` are null on per-kg rows, so
 * they average the rows that have them, with their own `N_`. AUMC_LAST, TLAG and
 * PCT_AUMC_EXTRAP are NOT aggregated by default — AUMC_INF is the headline moment,
 * lag is categorical-ish, %extrap is a per-row flag (AD-13). Columns absent from a
 * given results frame are skipped (see filter below).
 */
export const DEFAULT_AGG_PARAMS: ReadonlyArray<string> = [
  'CMAX', 'TMAX', 'AUC_LAST', 'AUC_INF', 'PCT_EXTRAP',
  'LAMBDA_Z', 'HALF_LIFE', 'CL', 'VZ',
  'AUMC_INF', 'MRT', 'VSS',
  'CL_ABS', 'VZ_ABS', 'VSS_ABS',
];

export interface AggregateOptions {
  /** Parameter columns to aggregate. Defaults to {@link DEFAULT_AGG_PARAMS}. */
  readonly params?: ReadonlyArray<string>;
}

/** Per-row-route AD-4 column holding the per-subject `'CL'` / `'CL/F'` basis. */
const CL_BASIS_COL = 'CL_BASIS';

/**
 * The SET of composite group-keys (over the `selectedKeys` tuple) whose rows pool
 * more than one distinct `stratifierCol` value — i.e. the aggregate cells that
 * genuinely mix bases. The single source of truth for BOTH the mixed-basis banner
 * (fires iff `.size > 0`) and the {@link blankMixedBasisFamily} post-pass (which
 * withholds the corrupted family cells for exactly these keys). Deriving the banner
 * and the blank from ONE scan is what keeps them coupled (they can never desync —
 * synthesis architecture spine, Open Q2).
 *
 * **Per composite-group, not frame-global (SF-2, multi-column-group-by).** Buckets
 * rows by the FULL `(k₁, k₂, …)` tuple via {@link compositeGroupKey}, then collects
 * every bucket that spans ≥2 stratifier values. Under a single-column selection this
 * reduces to the original per-group test. Under a `{Sex, Dose}` selection it asks
 * "does this *cell* mix CL basis" — NOT "does Sex alone mix basis" (which
 * over-triggers when Sex spans routes even where each Sex×Dose cell is single-basis,
 * the gap-18 noisy-warning lesson). A key that 1:1-determines the stratifier (e.g.
 * `ROUTE` → `CL_BASIS`: IV ⇒ CL, PO/SC ⇒ CL/F) never pools within a cell, so the set
 * stays empty. A `selectedKeys` member equal to `stratifierCol` makes the stratifier
 * constant within each bucket → empty set (redundant split).
 *
 * The keys are the SAME byte strings a caller recomputes over the post-aggregation
 * frame's group columns (DG `groupBy` preserves the group columns' type, so
 * `compositeGroupKey`'s `getString` yields identical strings on both frames) — the
 * match {@link blankMixedBasisFamily} relies on (AC-F2.3 numeric-key guard proves it).
 */
function mixedGroupKeys(
  results: DG.DataFrame,
  selectedKeys: ReadonlyArray<string>,
  stratifierCol: string,
): Set<string> {
  const mixed = new Set<string>();
  const col = results.col(stratifierCol);
  if (col === null) return mixed;
  const keyCols: DG.Column[] = [];
  for (const k of selectedKeys) {
    const kc = results.col(k);
    if (kc === null) return mixed; // a missing key column → cannot bucket → no split
    keyCols.push(kc);
  }
  if (keyCols.length === 0) return mixed;
  const firstSeen = new Map<string, string>();
  const n = results.rowCount;
  for (let i = 0; i < n; i++) {
    if (col.isNone(i)) continue;
    if (keyCols.some((kc) => kc.isNone(i))) continue;
    const g = compositeGroupKey(keyCols, i);
    const v = col.getString(i);
    const prev = firstSeen.get(g);
    if (prev === undefined) firstSeen.set(g, v);
    else if (prev !== v) mixed.add(g); // this composite cell spans ≥2 stratifier values
  }
  return mixed;
}

/**
 * Resolve the group-by key list for an aggregate and detect the mixed-basis cells.
 *
 * **The keys are EXACTLY the user's deduped selection — no auto-append**
 * (group-by-cl-basis-explicit-control F1, reversing per-row-route AC-C4). Silently
 * appending `CL_BASIS`/`DOSE_NORM` split one requested `Sex` row into four opaque
 * rows the user never asked for. Instead the aggregate honours the user's grouping
 * and the caller *withholds* (NE) the corrupted family cells in the pooled rows (see
 * {@link blankMixedBasisFamily}), driven by the key sets returned here.
 *
 * **`CL_BASIS` detection (per-row-route AC-C4 science, preserved).** Pooling true
 * `CL` (IV) and apparent `CL/F` (EV) into one `MEAN_CL` is a quantitative science
 * error (~1.67× CL mis-statement in a 50/50 IV:PO split at F = 25%). `mixedBasis` is
 * true when some composite cell spans both; `mixedBasisKeys` is exactly those cells.
 *
 * **No `DOSE_NORM` detection any more (nca-cl-bw-normalised-units AD-4).** `CL` /
 * `VZ` / `VSS` are per kg on every row whatever the engine dose basis, so a group
 * mixing `absolute` and `per_kg` engine doses pools one unit; the former
 * `mixedDoseNorm` withholding had no member left and was removed.
 *
 * `mixedBasis` is a convenience boolean (`.size > 0`) for the `MIXED_BASIS_AGGREGATE`
 * banner guard; the key set feeds the blank pass. Banner and blank share ONE
 * detection, so they cannot desync (Open Q2).
 */
export function resolveAggregateGroupKeys(
  results: DG.DataFrame,
  selectedKeys: string | ReadonlyArray<string>,
): {
  readonly keys: string[];
  readonly mixedBasis: boolean;
  readonly mixedBasisKeys: Set<string>;
} {
  // Widened single→multi seam (multi-column-group-by F3b): a composite selection
  // passes the full tuple; the mixed-basis detection buckets by the tuple, not by
  // one key. A single string keeps the original single-key behaviour. Dedupe so a
  // repeated key never double-groups.
  const selected = typeof selectedKeys === 'string' ? [selectedKeys] : [...selectedKeys];
  // The detection buckets by the USER'S tuple only — "does the user's grouping pool CL
  // and CL/F" in some cell.
  const userKeys = Array.from(new Set(selected));
  const mixedBasisKeys = mixedGroupKeys(results, userKeys, CL_BASIS_COL);
  return {keys: userKeys, mixedBasis: mixedBasisKeys.size > 0, mixedBasisKeys};
}

/**
 * The IV↔EV-incomparable parameter family (user-ratified 2026-07-10). When a pooled
 * aggregate cell mixes true CL (IV) and apparent CL/F (EV), the MEAN/SD/CV of these
 * params are not comparable across bases and are withheld (NE): `CL`/`VZ` because
 * `÷F` corrupts the units; `MRT`/`AUMC_INF` because they carry the absorption-moment
 * (MAT) term (`AUMC_INF` is `MRT`'s numerator input) — a *meaning* corruption, not
 * merely an exposure-level pool; `VSS` for conservative consistency (it is IV-only, so
 * a mixed cell's `MEAN_VSS` was already IV-only — blanking it hides an unrequested
 * sub-set). `AUC_INF`/`CMAX` are NOT here (identical units across route — a separate,
 * deferred cross-route *exposure* pool, RG-2). `HALF_LIFE`/`LAMBDA_Z`/`PCT_EXTRAP` are
 * NOT here — populated + carried by the banner flip-flop caveat (Decision A′). The absolute
 * companions `CL_ABS`/`VZ_ABS`/`VSS_ABS` carry the same ÷F as their per-kg siblings, so
 * they are withheld with them (nca-cl-bw-normalised-units AC-3.2).
 */
export const IV_EV_INCOMPARABLE_FAMILY: ReadonlyArray<string> =
  ['CL', 'VZ', 'VSS', 'MRT', 'AUMC_INF', 'CL_ABS', 'VZ_ABS', 'VSS_ABS'];

/**
 * Post-aggregation blank pass: withhold (write `FLOAT_NULL` → the grid renders an
 * empty NE cell) the `MEAN_/SD_/CV_PCT_<P>` cells of a scientific `family` for every
 * output row whose composite group key is in `mixedKeys` — the pooled cells that mix
 * incomparable bases. `N_<P>` and group `N` are left INTACT (the cohort still exists;
 * only the pooled central tendency is withheld — R2-5).
 *
 * **Where the blank lands (synthesis architecture spine).** `aggregateNcaResults` is
 * byte-untouched (AC-C4/AD-6 PKNCA 60/60 parity pin); this runs on the RETURNED `agg`
 * frame in the caller. The per-row key is recomputed over the same `userKeys` columns
 * `groupBy` carried through — DG preserves the group columns' type, so
 * {@link compositeGroupKey}'s `getString` yields the identical byte string on the
 * post-aggregation frame as on the pre-aggregation frame that built `mixedKeys` (the
 * coupling AC-F2.3's numeric-key fixture proves, not merely asserts).
 *
 * Mutates `agg` in place (the freshly built aggregate — never the source frame).
 * A no-op when `mixedKeys` is empty (the parity-pin fast path).
 */
export function blankMixedBasisFamily(
  agg: DG.DataFrame,
  mixedKeys: ReadonlySet<string>,
  userKeys: ReadonlyArray<string>,
  family: ReadonlyArray<string>,
): void {
  if (mixedKeys.size === 0) return;
  const keyCols: DG.Column[] = [];
  for (const k of userKeys) {
    const kc = agg.col(k);
    if (kc === null) return; // key column not carried through → cannot match → no-op
    keyCols.push(kc);
  }
  if (keyCols.length === 0) return;
  // The MEAN/SD/CV stat columns present in this aggregate for the family params. N_<P>
  // is deliberately excluded — the contributor count survives the withhold.
  const targets: DG.Column[] = [];
  for (const p of family) {
    for (const stat of [`MEAN_${p}`, `SD_${p}`, `CV_PCT_${p}`]) {
      const c = agg.col(stat);
      if (c !== null) targets.push(c);
    }
  }
  if (targets.length === 0) return;
  const n = agg.rowCount;
  for (let i = 0; i < n; i++) {
    if (!mixedKeys.has(compositeGroupKey(keyCols, i))) continue;
    for (const c of targets) c.set(i, DGNS.FLOAT_NULL);
  }
}

/**
 * Aggregate a per-subject NCA results DataFrame into a per-group summary.
 *
 * For each requested parameter `P` present in `results`, emits `MEAN_P`
 * (arithmetic), `SD_P` (sample), and `CV_PCT_P` (= 100·SD/MEAN). Adds one
 * `N` row-count column and carries the grouping column(s) through. Row order is
 * whatever the native groupBy yields (not contractually sorted in the spike).
 *
 * `groupBy` accepts a single column name (the original spike contract) or a
 * list of names — the multi-key form backs the per-row-route mixed-basis guard
 * (see {@link resolveAggregateGroupKeys}). Single-key behaviour is unchanged.
 *
 * @throws if any group column is not a column of `results`.
 */
export function aggregateNcaResults(
  results: DG.DataFrame,
  groupBy: string | ReadonlyArray<string>,
  options: AggregateOptions = {},
): DG.DataFrame {
  const groupCols = typeof groupBy === 'string' ? [groupBy] : [...groupBy];
  for (const g of groupCols) {
    if (results.col(g) === null)
      throw new Error(`aggregateNcaResults: group column '${g}' not found in results.`);
  }

  const params = (options.params ?? DEFAULT_AGG_PARAMS)
    .filter((p) => results.col(p) !== null);

  // ---------------------------------------------------------------------
  // Aggregate the SOURCE frame directly (na.rm semantics).
  //
  // **A NaN → FLOAT_NULL pre-pass used to live here and has been DELETED**
  // (nca-dg-native-missing-values F3, 2026-08-03). It cloned every parameter
  // column so that unfit values reached `groupBy` as `FLOAT_NULL` rather than
  // as IEEE NaN, because native DG `avg` does NOT skip a raw NaN — one NaN in a
  // group nulls that group's whole mean, while `FLOAT_NULL` is excluded cleanly
  // (spike-verified here; the probe later extended the finding from `avg` to the
  // entire stats surface — min/max/avg/sum/stdev, PR-3).
  //
  // That clone was a COMPENSATION FOR A DEFECT ELSEWHERE: the result writers
  // emitted the wrong sentinel. They now emit `DG.FLOAT_NULL` at construction
  // (`adapter-utils.f64Column`), so every parameter column arrives already
  // correct and the clone is dead code. **This is not a change to `groupBy`
  // SEMANTICS** — the N-schema (`N` = cohort row count, `N_<P>` = per-param
  // `valueCount` of contributors), the unit/label carry-through and the %CV
  // guard are all untouched. Three prior blueprints marked this module
  // DO-NOT-MODIFY on exactly those grounds; none of them is affected.
  //
  // **Accepted regime — unbounded `na.rm`.** DG `avg`/`stdev` apply no
  // minimum-n bound, and none is added here. PKNCA — the parity reference
  // NFR-04/05 binds us to — DOES bound it (`PKNCA.options("max.missing")`
  // returns `NA` when too many values are missing). So a group with 7 of 8
  // subjects unfit yields a one-animal "group mean" that renders identically to
  // an eight-animal one. `N_<P>` discloses it HERE; the newly-exposed DG-native
  // surfaces (ad-hoc grid group-by, `col.stats`, filters, export) carry no such
  // companion. Recorded as a stated limitation with its own filed gap, not a
  // silent choice — bounding it would change `MEAN_*`/`SD_*` output and is
  // outside the user's fixed decision, which is about the REPRESENTATION of
  // absence.
  // ---------------------------------------------------------------------
  const gb = results.groupBy(groupCols);
  for (const p of params) {
    gb.avg(p, `MEAN_${p}`);
    gb.stdev(p, `SD_${p}`);
    gb.valueCount(p, `N_${p}`);
  }
  gb.count('N');
  const agg = gb.aggregate();

  // Per-param post-pass: carry-through `tags.units` + CL/F-vs-CL `tags.label`
  // from the source results column, then compute the derived `CV_PCT_<P>`.
  //
  // Units: MEAN_<P> and SD_<P> inherit the source unit (mean / SD of mg/L is
  // still mg/L); CV_PCT_<P> is always '%'.
  // Label: the prefixed convention ("Mean CL/F", "SD CL/F", "%CV CL/F")
  // preserves the IV-vs-extravascular distinction in the aggregated grid
  // header so a reader doesn't confuse `MEAN_CL` (intravenous) with
  // `MEAN_CL/F` (oral/SC bioavailability-uncorrected).
  // CV guard: a null (FLOAT_NULL) / zero group mean yields NaN %CV rather
  // than ±Inf or a garbage ratio. `isNone` catches FLOAT_NULL; `Number.isFinite`
  // would NOT, because FLOAT_NULL (2.68e-34) is finite.
  for (const p of params) {
    const meanCol = agg.col(`MEAN_${p}`);
    const sdCol = agg.col(`SD_${p}`);
    if (meanCol === null || sdCol === null) continue;

    const sourceCol = results.col(p)!;
    const sourceUnits = sourceCol.tags['units'];
    const sourceLabel = sourceCol.tags['label'];
    // Carry the source unit through ONLY when it is a real unit string. An
    // UNSET DG tag reads back as `null` (not `undefined`) — e.g. a hand-built
    // frame. Guard on both so an untagged source leaves MEAN/SD untagged (no
    // fallback unit) rather than carrying `null`.
    if (sourceUnits !== undefined && sourceUnits !== null) {
      meanCol.setTag('units', sourceUnits);
      sdCol.setTag('units', sourceUnits);
    }
    // Guard BOTH undefined and null: an UNSET DG tag reads back as `null` (not
    // `undefined`) — same as the units guard above. Only CL/VZ/VSS carry a
    // per-subject `label`; every other param has none, so without the null check
    // we would write the literal string "Mean null" / "SD null" as the label.
    if (sourceLabel !== undefined && sourceLabel !== null) {
      meanCol.setTag('label', `Mean ${sourceLabel}`);
      sdCol.setTag('label', `SD ${sourceLabel}`);
    }

    const cv = new Float64Array(agg.rowCount);
    for (let i = 0; i < agg.rowCount; i++) {
      if (meanCol.isNone(i) || sdCol.isNone(i)) {
        cv[i] = NaN;
        continue;
      }
      const m = meanCol.get(i) as number;
      const s = sdCol.get(i) as number;
      cv[i] = m === 0 ? NaN : (100 * s) / m;
    }
    // F1 scope: `CV_PCT_<P>` is a Float64 column of a RESULT (aggregated) frame,
    // built three statements below the workaround F3 deleted. Without `f64Column`
    // here the migration would ship its own trigger symptom — the literal text
    // "NaN" in a grid cell — on the aggregated grid.
    const cvCol = f64Column(`CV_PCT_${p}`, cv);
    cvCol.setTag('units', '%');
    if (sourceLabel !== undefined && sourceLabel !== null)
      cvCol.setTag('label', `%CV ${sourceLabel}`);
    agg.columns.add(cvCol);
  }

  agg.name = 'NCAResultsByGroup';
  return agg;
}
