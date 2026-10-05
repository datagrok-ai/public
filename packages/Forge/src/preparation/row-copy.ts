import * as DG from 'datagrok-api/dg';

/** A copy of the rows of [frame] in [mask]. A masked clone of a text column keeps every category of the source,
 * the empty one included: the copy drops those no kept row uses. */
export function rowCopy(frame: DG.DataFrame, mask: DG.BitSet): DG.DataFrame {
  const copy = frame.clone(mask);
  for (const col of copy.columns.toList())
    compactCategories(col);
  return copy;
}

/** {@link rowCopy} for one column. */
export function columnRowCopy(col: DG.Column, mask: DG.BitSet): DG.Column {
  const copy = col.clone(mask);
  compactCategories(copy);
  return copy;
}

/** Drops the categories no row of [col] uses; call it on Forge's own copies only, never on the user's columns. */
export function compactCategories(col: DG.Column): void {
  if (col.type === DG.COLUMN_TYPE.STRING)
    col.compact();
}
