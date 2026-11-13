/**
 * Internal adapter helpers shared across modules in `src/adapter/`.
 *
 * **Re-export policy (updated 2026-08-03, nca-dg-native-missing-values F6).**
 * This file used to declare itself package-private ("NOT re-exported from
 * `./index.ts`"). That is no longer true and the header would otherwise deny
 * what the file does: {@link numAt} and {@link finiteMinMax} are the sanctioned
 * DG-column READ accessors and are re-exported from `./index.ts` so the `views/`
 * layer consumes the same spelling the adapter does. {@link f64Column} is
 * re-exported as well — `views/nca-mapping-pane.ts` genuinely builds one RESULT
 * column (`TLAST`), and withholding the helper is precisely what let that site
 * keep writing a raw NaN until the Arm A inter-rater pass found it (2026-08-03).
 * Only {@link stringifyKey} stays adapter-internal (imported by path).
 *
 * `views/` must never host these helpers: that layer carries a pre-existing
 * import cycle (`nca-launch → nca-mapping-pane → nca-add-to-study →
 * nca-analysis-open → nca-launch`), and a top-level binding there that reads
 * another `views/` module at init time takes the whole webpack bundle down.
 * This file's only import is `datagrok-api/dg`, which cannot participate.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

/**
 * Stringify a cell value (subjectId, route, etc.) for use as a Map key or
 * error-message identifier. `null` / `undefined` → empty string; everything
 * else via `String(v)`. Shared between validators.ts (DUPLICATE_SUBJECT_TIME
 * scan) and profile-index.ts (groupBy bucket key).
 */
export function stringifyKey(v: unknown): string {
  if (v === null || v === undefined) return '';
  return String(v);
}

// ---------------------------------------------------------------------------
// The DG missing-value boundary (nca-dg-native-missing-values, user decision D1).
//
// THE TWO SENTINELS ARE NOT THE SAME BITS AND NOT THE SAME PREDICATE:
//
//   in-package / sci-comp  →  IEEE-754 NaN.        `Number.isNaN` catches it.
//   Datagrok column        →  DG.FLOAT_NULL        `Column.isNone` catches it.
//                             = 2.6789344063684636e-34 (const.d.ts:51)
//
// MEASURED (probe 2026-07-30, `probe_discharged` in the cycle state):
//   - `isNone(rawNaN) === false`. A raw NaN in a DG column is a PRESENT VALUE to
//     every DG-native consumer (grid sort, filter, `col.stats`, group-by, export)
//     and POISONS every statistic on the column (PR-3: min/max/avg/sum/stdev all
//     NaN from a single cell). `FLOAT_NULL` is excluded cleanly instead.
//   - `Column.get()` / `getNumber()` return the RAW sentinel even when
//     `isNone(i)` is true — the datagrok-api JSDoc is wrong (PR-2). Only
//     `isNone()` is trustworthy on the Column API. `getString()` and
//     `GridCell.value` DO honour null (PR-9).
//
// THE RISK IS BIDIRECTIONAL. Per IEEE 754-2019 §5.11 a raw NaN is fail-SAFE
// under every relational operator (`NaN > 0`, `NaN < t`, `NaN <= n` are all
// false) and under `isFinite`. `DG.FLOAT_NULL` is FINITE and POSITIVE, hence
// fail-DANGEROUS under `> 0`, `< threshold`, `<= floor`, `!== 0` and `isFinite`.
// So this boundary makes DG-native reads safer and in-package raw reads MORE
// dangerous. That is why {@link numAt} / {@link finiteMinMax} exist: any read of
// a result column must go through an `isNone`-first accessor, never a bare
// `get()` + `isFinite`.
// ---------------------------------------------------------------------------

/**
 * Build a Float64 DG column from a numeric array, translating the in-package /
 * sci-comp IEEE-754 `NaN` sentinel to the platform's `DG.FLOAT_NULL`.
 *
 * **This is the ONLY sanctioned way to build a numeric result column.** Every
 * `Column.fromFloat64Array` call in a result / derived frame routes through
 * here, so the conversion is part of CONSTRUCTION rather than a post-pass an
 * author can forget (the half-migrated state this replaced had four piecemeal
 * fixes in three different spellings).
 *
 * Semantics (each is an acceptance criterion of the migration):
 *  - **AC-1** `NaN` → `DG.FLOAT_NULL`; the cell reads `isNone`.
 *  - **AC-2** a value that is ALREADY `DG.FLOAT_NULL` passes through unchanged
 *    and still reads `isNone` (the dual-check becomes structural).
 *  - **AC-3** `±Infinity` passes through UNCHANGED — the predicate is
 *    `Number.isNaN`, **not** `!Number.isFinite`. Translating it would be an
 *    undeclared behaviour change: `bioavailability.ts` can emit `+Infinity` for
 *    a zero-exposure IV reference, that value is aggregated as PRESENT today,
 *    and nulling it here would both change `N_F_PCT` and manufacture a false
 *    "no IV reference" disclosure downstream. The `ivRef === 0` defect is filed
 *    as its own gap; this helper does not silently absorb it. (User decision,
 *    2026-08-03.)
 *  - **AC-4** COPIES into a fresh `Float64Array`; never mutates `values`.
 *    `fromFloat64Array` ADOPTS the array it is given as zero-copy storage, and
 *    at least one call site passes a caller-owned array (`results-writer.ts`
 *    `c0Effective` is the controller's `c0s`).
 *  - **AC-4b** a `null` / `undefined` / non-numeric element is a PROGRAMMING
 *    ERROR, not a missing value, and throws. `Float64Array.from([null])` yields
 *    `0`, which on a semi-log concentration axis is a real point at zero rather
 *    than a gap — this helper must not silently decide that question.
 *  - **AC-6** real values round-trip BIT-EXACT (no perturbation of a
 *    PKNCA-parity value — NFR-04).
 */
export function f64Column(name: string, values: ArrayLike<number>): DG.Column {
  const n = values.length;
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const v = values[i];
    // AC-4b — throw rather than coerce. `typeof NaN === 'number'`, so a genuine
    // missing value still takes the branch below; only null/undefined/string/
    // object elements land here.
    if (typeof v !== 'number') {
      throw new Error(
        `f64Column('${name}'): element ${i} is ${v === null ? 'null' : typeof v}, not a number. ` +
        `A missing value must reach this helper as NaN or DG.FLOAT_NULL — null/undefined would ` +
        `coerce to 0 and read as a real datum.`);
    }
    // AC-1 / AC-2 / AC-3 / AC-6 — NaN is the only value rewritten.
    out[i] = Number.isNaN(v) ? DGNS.FLOAT_NULL : v;
  }
  return DGNS.Column.fromFloat64Array(name, out);
}

/**
 * Read a numeric DG cell as a JS number, mapping the platform's null to `NaN`
 * — i.e. back across the boundary {@link f64Column} writes.
 *
 * **`isNone`-FIRST, and that is the whole point.** `Column.get()` returns the
 * RAW sentinel even for a null cell (PR-2), so a bare `get()` + `Number.isFinite`
 * admits `DG.FLOAT_NULL` (2.68e-34) as a real value. Every read of a result
 * column goes through here.
 *
 * Hoisted verbatim from three BYTE-IDENTICAL private copies
 * (`dose-proportionality.ts`, `toxicokinetics.ts`, `sparse-nca.ts`) plus an
 * inline fourth in `bioavailability.ts` — not a new abstraction, three copies
 * collapsed to one (nca-dg-native-missing-values F1).
 */
export function numAt(col: DG.Column, row: number): number {
  return col.isNone(row) ? NaN : Number(col.get(row));
}

/**
 * Normalise EITHER missing-value spelling onto the platform null, for a value
 * being copied into a chart-side numeric array.
 *
 * The `Number.isNaN` half is LOAD-BEARING: `Column.isNone()` returns false for a
 * raw IEEE NaN, so an unguarded NaN reaches the grid as the literal text "NaN"
 * and collapses `min`/`max`/`avg`/`sum`/`stdev` on the whole column. The
 * `=== FLOAT_NULL` half is belt-and-braces — `Math.fround(FLOAT_NULL) ===
 * FLOAT_NULL` exactly, so a Float32 narrowing PRESERVES the sentinel and its
 * `isNone` (measured 2026-07-30, probe PR-10).
 *
 * **Correction to NOTE-12 lore:** that note attributed a spurious ~0 chart point
 * to `FLOAT_NULL` "not surviving" the Float32 narrowing. The measurement
 * contradicts that mechanism — `8849b44` was a two-part fix (structural: drop
 * EVID=1 dose rows upstream; plus this residual guard) and the evidence says the
 * STRUCTURAL half did the work. Do not re-derive the narrowing story.
 *
 * Hoisted from three hand-spelled copies in the chart builders so Arm 2 of the
 * sentinel drift-guard is an executable pin rather than a source-grep
 * (nca-dg-native-missing-values F4 Arm 2).
 */
export function toPlatformNull(v: number): number {
  return (v === DGNS.FLOAT_NULL || Number.isNaN(v)) ? DGNS.FLOAT_NULL : v;
}

/**
 * `log10(v)` for a genuinely positive value, else the platform null — the
 * pre-derived semi-log Y channel used by the chart frames.
 *
 * The explicit `!== FLOAT_NULL` leg is separately load-bearing, and for a
 * DIFFERENT reason than in {@link toPlatformNull}: `DG.FLOAT_NULL` is a tiny
 * POSITIVE number, so it PASSES `> 0`, and without the check this would store a
 * REAL `log10(2.68e-34) ≈ −33.6` instead of a gap — a point 33 decades below the
 * data that silently rescales the axis. A raw `NaN` needs no separate leg
 * (`NaN > 0` is false).
 */
export function log10OrNull(v: number): number {
  return (v !== DGNS.FLOAT_NULL && v > 0) ? Math.log10(v) : DGNS.FLOAT_NULL;
}

/** One anomalous numeric cell found by {@link scanNumericSentinels}. */
export interface SentinelScanEntry {
  readonly column: string;
  readonly row: number;
  /** `raw-nan` — an IEEE NaN sitting in a DG column, i.e. a value the platform
   *  reads as PRESENT that no consumer can use. This is the defect the migration
   *  removes and it must FAIL the suite.
   *  `non-finite` — `±Infinity`. REPORT-ONLY: it is a separate, user-deferred
   *  defect (`bioavailability`'s zero-exposure IV reference), and `f64Column`
   *  deliberately passes it through (AC-3). */
  readonly kind: 'raw-nan' | 'non-finite';
  readonly value: number;
}

/**
 * Scan every numeric column of a built frame for values that should not be
 * there. The permanent data-driven drift guard behind the sentinel migration
 * (nca-dg-native-missing-values F4 Arm 1 / AC-12 / AC-12b).
 *
 * **Rule 16 — this is meant to run over REAL generated pipeline output**, not
 * over a mirror of a writer's logic. A writer that forgets `f64Column` produces
 * a frame this scan fails on, naming column + row; no amount of source-level
 * review does that as reliably.
 *
 * Two buckets, deliberately asymmetric:
 *  - `raw-nan` FAILS. A raw NaN is invisible to `isNone`, is counted as present
 *    by group-by / filters / export, and POISONS every statistic on its column
 *    (probe PR-3: min/max/avg/sum/stdev all NaN from one cell).
 *  - `non-finite` is REPORTED but does not fail. Reverting AC-3 left the
 *    detector able to certify a frame as clean while it holds a `+Infinity`, and
 *    the migration makes that WORSE rather than neutral: on a column holding
 *    both, today's NaN dominates every statistic and MASKS the Inf, whereas
 *    afterwards the NaN is excluded and the `+Inf` becomes the sole surviving
 *    poisoner. Note the failure modes differ and must not be conflated — NaN
 *    collapses `min`/`max`/`avg`/`sum`/`stdev`; `+Inf` leaves `min` finite and
 *    sends `avg`/`sum` to `+Inf`.
 *
 * Int columns are skipped: `Int32Array` cannot hold NaN or Infinity (ECMAScript
 * `ToInt32` maps all three to `+0`), so there is nothing for this scan to find
 * there — and that same coercion is why an Int32 null must be written as
 * `DG.INT_NULL` explicitly.
 */
export function scanNumericSentinels(df: DG.DataFrame): SentinelScanEntry[] {
  const out: SentinelScanEntry[] = [];
  for (const col of df.columns as Iterable<DG.Column>) {
    if (col.type !== DGNS.COLUMN_TYPE.FLOAT && col.type !== DGNS.COLUMN_TYPE.QNUM) continue;
    const raw = col.getRawData();
    if (!(raw instanceof Float64Array) && !(raw instanceof Float32Array)) continue;
    for (let i = 0; i < raw.length; i++) {
      const v = raw[i];
      if (Number.isNaN(v)) out.push({column: col.name, row: i, kind: 'raw-nan', value: v});
      else if (!Number.isFinite(v)) out.push({column: col.name, row: i, kind: 'non-finite', value: v});
    }
  }
  return out;
}

/** Render {@link scanNumericSentinels} output as a one-line failure message —
 *  frame + column + row, so a red suite says exactly which writer regressed. */
export function describeSentinelEntries(
  frameName: string, entries: ReadonlyArray<SentinelScanEntry>,
): string {
  return entries
    .map((e) => `${frameName}.${e.column}[${e.row}]=${e.kind === 'raw-nan' ? 'NaN' : e.value}`)
    .join(', ');
}

/**
 * `[min, max]` of the finite values in one column, or `null` when the column is
 * absent / has no qualifying value. The shared axis-range scan.
 *
 * The guard is a TRIPLE: `isNone` → `=== FLOAT_NULL` → `!isFinite`. The second
 * leg is deliberate belt-and-braces for a column type where `isNone` and the raw
 * bits could disagree (e.g. a Float32 narrowing) — do not "simplify" it away.
 *
 * `positiveOnly` additionally drops `v <= 0`, for a LOG axis where a zero or
 * negative datum has no position. Without it a real `0` would reach `Math.log10`
 * and produce `-Infinity`; with the sentinel migration, an unguarded
 * `isFinite(v) && v > 0` would instead admit `FLOAT_NULL` and stretch a semi-log
 * axis ~40 decades. Both failure modes are closed here, in one place.
 *
 * Hoisted from `views/nca-mapping-pane.ts` into `adapter/` so the study-overview
 * thumbnail can reuse it WITHOUT `views/ → views/` import edges (see the module
 * header). `positiveOnly` is additive; the three original callers keep the
 * default.
 */
export function finiteMinMax(
  col: DG.Column | null, positiveOnly = false,
): {min: number; max: number} | null {
  if (col === null) return null;
  let mn = Infinity; let mx = -Infinity;
  for (let i = 0; i < col.length; i++) {
    if (col.isNone(i)) continue;
    const v = col.get(i) as number;
    if (v === DGNS.FLOAT_NULL || !Number.isFinite(v)) continue;
    if (positiveOnly && v <= 0) continue;
    if (v < mn) mn = v; if (v > mx) mx = v;
  }
  return (Number.isFinite(mn) && Number.isFinite(mx)) ? {min: mn, max: mx} : null;
}

/** Name of the BOOL column {@link ensureBoolBlq} derives from a non-bool BLQ flag column. */
export const DERIVED_BLQ_COL = '~BLQ';

/**
 * The frame's BLQ column as a BOOL column: accepts bool, 0/1, or Y / YES / TRUE / BLQ / BLOQ /
 * <LLOQ text (case-insensitive). A non-bool column gets a derived {@link DERIVED_BLQ_COL} column
 * added to `df`; the returned name is the one to read. The ONE reading of a BLQ flag for the
 * pipeline and for both sparse paths (BUG-19 AD-8; moved from `pipeline/pipeline-compute.ts`).
 */
export function ensureBoolBlq(df: DG.DataFrame, name: string): string {
  const col = df.col(name)!;
  if (col.type === DGNS.COLUMN_TYPE.BOOL) return name;
  const truthy = /^(1|y|yes|true|blq|bloq|<lloq)$/i;
  if (df.col(DERIVED_BLQ_COL) !== null) df.columns.remove(DERIVED_BLQ_COL);
  df.columns.add(DGNS.Column.fromList(DGNS.COLUMN_TYPE.BOOL, DERIVED_BLQ_COL,
    Array.from({length: col.length}, (_, i) => !col.isNone(i) && truthy.test(String(col.get(i)).trim()))));
  return DERIVED_BLQ_COL;
}
