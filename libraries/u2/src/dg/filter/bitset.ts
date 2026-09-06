/* DataFrame evaluation on the platform: the core masks run over the frame's own columns
   (`getRawData`, `categories`, `get`), and the result becomes one `DG.BitSet` in a single copy. */
import * as DG from 'datagrok-api/dg';
import {toMask} from '../../core/filter/mask.js';
import type {FilterGroup} from '../../core/filter/model.js';
import type {MaskFrameLike} from '../../core/filter/mask.js';

/** A `DG.DataFrame` as the structural frame `Filters.toMask` evaluates over. */
export function frameLike(df: DG.DataFrame): MaskFrameLike {
  return {rowCount: df.rowCount, column: (name) => df.columns.byName(name)};
}

/** `toMask` over a DataFrame, delivered as a platform BitSet — `df.filter.and(bitset)` applies it. */
export async function toBitSet(df: DG.DataFrame, root: FilterGroup,
  options?: {signal?: AbortSignal, now?: Date}): Promise<DG.BitSet> {
  const mask = await toMask(frameLike(df), root, options);
  return DG.BitSet.fromBytes(mask.bits.buffer as ArrayBuffer, mask.length);
}
