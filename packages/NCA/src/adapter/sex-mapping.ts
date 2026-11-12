/**
 * Value→sex (M/F/U) mapping for the sex stratifier.
 *
 * **Why this is not a presentation-only relabel.** Stratification is by
 * value-equality: the decorators (`computeBioavailability`, dose-proportionality,
 * TK) and the NONMEM composite-pool cell key (`composite-pool.ts`) key strata on
 * the raw sex-column value. A study coded `1`/`2` has no universal `1=male`
 * convention, and a string column can carry token variants (`M` / `Male` /
 * `MALE`). Collapsing those to canonical M/F/U strata CHANGES which subjects pool
 * together — i.e. it can change F% — so it must be an EXPLICIT, user-controlled
 * map, never a silent guess (peer-review FLAWED finding). This module is the pure
 * core of that map: a heuristic default for unambiguous tokens, and an in-place
 * applier that rewrites the designated sex column to its resolved strata.
 *
 * **Honest uncertainty (CLAUDE.md §4).** Ambiguous numeric codes are NOT guessed
 * — {@link defaultSexAssignment} returns null and the UI leaves them unassigned
 * and blocks Compute until the scientist resolves them from study metadata. A
 * blank cell maps to a VISIBLE `U` stratum, never a silent drop.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

/** Resolved sex stratum. `U` = unknown / undifferentiated (a real stratum, not a
 *  drop). */
export type SexCode = 'M' | 'F' | 'U';

/** Column-name synonyms for the sex column — the mapping pane's Subject Sex
 *  suggestion and the pipeline mapping's detection. Wider than `/^sex$/i` so a
 *  `Gender` / coded `SEXCD` column is found, but anchored (`^…$`) so a substring
 *  like `Unisex` does not match. */
export const SEX_SYNONYMS = /^(sex|gender|sexcd|sex_?code|sexc)$/i;

/**
 * Heuristic default sex assignment for a raw token. Returns:
 *  - `M` / `F` / `U` for UNAMBIGUOUS text tokens (case-folded, trimmed);
 *  - `U` for a blank / missing value (a visible Unknown stratum, never a drop);
 *  - `null` for an AMBIGUOUS numeric code (`1` / `2` / `0`) or any unrecognized
 *    token — the caller must leave it unassigned (no `1=male` convention to
 *    guess from; peer-review LBC-3).
 */
export function defaultSexAssignment(raw: string): SexCode | null {
  const t = raw.trim().toLowerCase();
  if (t === '') return 'U';
  if (t === 'm' || t === 'male') return 'M';
  if (t === 'f' || t === 'female') return 'F';
  if (t === 'u' || t === 'unk' || t === 'unknown' || t === 'undifferentiated') return 'U';
  return null;
}

/**
 * Distinct raw values of a column over ALL rows (a missing / NA cell collapses
 * to the empty-string `''` bucket — the "blank" value). Both arms (IV and EV)
 * are scanned by the caller passing the whole column, so no token is missed and
 * the IV↔EV value-equality join cannot break. Insertion order is preserved for a
 * stable UI.
 */
export function uniqueSexValues(col: DG.Column): string[] {
  const seen = new Set<string>();
  const out: string[] = [];
  const n = col.length;
  for (let i = 0; i < n; i++) {
    const v = col.isNone(i) ? '' : col.getString(i);
    if (!seen.has(v)) {seen.add(v); out.push(v);}
  }
  return out;
}

/**
 * Apply a value→sex map to `df`'s designated sex column IN PLACE, replacing it
 * with a STRING column of the same name holding the resolved M/F/U strata.
 *
 * A raw value NOT present in `map` is left UNCHANGED — so the transform is
 * idempotent across re-Compute (after the first apply the column holds M/F/U,
 * which the raw-keyed map no longer matches, so a second apply is a no-op). The
 * Compute-disable gate guarantees every scanned value IS mapped at apply time,
 * so "left unchanged" only ever fires for already-resolved strata.
 *
 * Replacing (not mutating) the column lets a coded INT column (`1`/`2`) become a
 * STRING M/F column. Row order is preserved, so index-aligned carriage / pooling
 * are unaffected. Pure w.r.t. every other column.
 *
 * @throws if `sexColName` is absent from `df`.
 */
export function applySexMap(
  df: DG.DataFrame, sexColName: string, map: ReadonlyMap<string, SexCode>,
): void {
  const col = df.col(sexColName);
  if (col === null)
    throw new Error(`applySexMap: sex column '${sexColName}' not found.`);
  const n = col.length;
  const out = new Array<string>(n);
  for (let i = 0; i < n; i++) {
    const raw = col.isNone(i) ? '' : col.getString(i);
    out[i] = map.get(raw) ?? raw;
  }
  df.columns.remove(sexColName);
  df.columns.add(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, sexColName, out));
}

/** True when at least one entry resolves a raw value to a DIFFERENT string —
 *  i.e. applying the map would change the column. An all-identity map (e.g.
 *  `M→M` / `F→F` on an already-canonical column) is a no-op, so the caller can
 *  skip the clone + apply (preserves the Slice-A fast path for clean data). */
export function sexMapChangesValues(map: ReadonlyMap<string, SexCode>): boolean {
  for (const [raw, code] of map) if (raw !== code) return true;
  return false;
}
