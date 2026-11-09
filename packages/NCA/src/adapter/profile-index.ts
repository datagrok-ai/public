/**
 * Profile-index builder — groups observation rows by subject.
 *
 * Scope (UC-01): single occasion, single analyte, single profile per subject.
 * The interface §3 long-term `ProfileIndex` carries `(subjectCode, analyteCode,
 * occasionCode, compoundCode)` plus `dosingRow`; the UC-01 shape is the
 * narrowed `{subjectKey, observationRows}` pair (no dosing DataFrame in this
 * slice — dose is per-row or scalar from the mapping form). Multi-analyte and
 * multi-occasion routing are OQ-2 / future work.
 *
 * Determinism: profiles are returned sorted lexicographically by
 * `subjectKey` so that re-imports of the same CSV produce the same order, and
 * the results DataFrame (Phase 3 results-writer) has a stable row order across
 * runs. Within a profile the row indices preserve their original DataFrame
 * order; per-profile time-sort is the (Phase 3) profile-inputs extractor's
 * responsibility (and only when defensive sort is enabled).
 */

import type * as DG from 'datagrok-api/dg';

import type {ColumnContract} from './column-contract';
import {stringifyKey} from './adapter-utils';

/**
 * UC-01-scoped profile index. One entry per subject. `observationRows` is an
 * `Int32Array` of row indices into the observations DataFrame, in the order
 * encountered during the single-pass scan. Sort by time is performed
 * downstream by the profile-inputs extractor when its `defensiveSort` option
 * is on (UC-01 default).
 */
export interface ProfileIndex {
  readonly subjectKey: string;
  readonly observationRows: Int32Array;
}

/**
 * Build a list of `ProfileIndex` from one observations DataFrame.
 *
 * Single pass O(n_rows) over the subject column; one `Map<string, number[]>`
 * keyed by stringified subject id. Empty subjects (null / undefined cell) are
 * grouped under the empty-string key — validators (Phase 3 validators.ts) are
 * expected to flag missing-subject conditions upstream when present; this
 * builder is intentionally lenient (it never throws).
 *
 * Lexicographic sort uses `Array.sort()` on the bucket keys after the scan;
 * for the UC-01 fixtures (12 theoph + 6 indometh subjects) this is sub-ms.
 */
export function buildProfileIndices(
  df: DG.DataFrame,
  contract: ColumnContract,
): ProfileIndex[] {
  const subjectColName = contract.observations.subjectId.name;
  const subj = df.col(subjectColName);
  if (subj === null) return [];

  const n = df.rowCount;
  const buckets = new Map<string, number[]>();
  for (let i = 0; i < n; i++) {
    const sk = stringifyKey(subj.get(i));
    let rows = buckets.get(sk);
    if (rows === undefined) {
      rows = [];
      buckets.set(sk, rows);
    }
    rows.push(i);
  }

  const keys = Array.from(buckets.keys()).sort();
  const out: ProfileIndex[] = new Array(keys.length);
  for (let k = 0; k < keys.length; k++) {
    const rows = buckets.get(keys[k])!;
    out[k] = {
      subjectKey: keys[k],
      observationRows: Int32Array.from(rows),
    };
  }
  return out;
}

