/**
 * Composite group-key helper (multi-column-group-by F3b) — the single shared
 * join every composite-grouping consumer uses so their keys are byte-identical.
 *
 * **Why one helper (rule 5 / simplicity rationale).** Composite grouping keys a
 * (compound × dose × sex …) tuple. Several call sites must produce the SAME key
 * string for the same row or a silent cross-consumer mismatch appears — the exact
 * defect class this cycle opens with (the drill-in mask fails to match, or a
 * per-group warning is misattributed). The pre-existing `.join('')` `keyStr` in
 * `nca-analysis-open.ts` was a live instance of that hazard (`['IV','bolus']`
 * collides with `['IVbolus','']`). Every consumer now joins with {@link KEY_SEP}.
 *
 * Consumers: the Mean±SD chart composite key, the compound-tier dedup / curve
 * keying, the drill-in tuple reconstruction, `computeFpctByCompound` /
 * `buildCurvesByCompound` key columns, `buildSummaryCache`'s per-group warning
 * join, and `isMixedStratifier`'s composite-cell bucketing.
 */

import type * as DG from 'datagrok-api/dg';

/** U+0001 join delimiter — a control char absent from any real route / compound /
 *  dose / sex value, so `['IV','bolus']` can never collide with `['IVbolus','']`.
 *  It renders invisibly in diffs; it is NOT an empty join. */
export const KEY_SEP = String.fromCharCode(1);

/**
 * The byte-identical composite group key for row `i` over `cols`. Each column's
 * value is stringified (a null cell → `''`, matching the pre-existing per-group
 * key behaviour) and joined by {@link KEY_SEP}. Order-sensitive: the caller must
 * pass the columns in a stable order (the same order every consumer uses).
 */
export function compositeGroupKey(cols: ReadonlyArray<DG.Column>, i: number): string {
  const parts = new Array<string>(cols.length);
  for (let k = 0; k < cols.length; k++) {
    const c = cols[k];
    parts[k] = c.isNone(i) ? '' : c.getString(i);
  }
  return parts.join(KEY_SEP);
}
