import * as DG from 'datagrok-api/dg';

export function toDataFrame(rows: Record<string, any>[]): DG.DataFrame {
  return DG.DataFrame.fromObjects(rows) ?? DG.DataFrame.create();
}

/**
 * Builds a fixed-length row set aligned to the input order via each result's client-provided id.
 * The server may drop or reorder rows (failed/filtered molecules), so positional joins are unsafe;
 * every row carries the same keys (missing rows are null-filled) to keep a stable column set.
 */
export function scatterRows(
  count: number,
  entries: {index: number; row: Record<string, any>}[],
): Record<string, any>[] {
  const keys = new Set<string>();
  for (const e of entries)
    for (const k in e.row) keys.add(k);
  const empty = (): Record<string, any> => Object.fromEntries([...keys].map((k) => [k, null]));
  const rows = Array.from({length: count}, empty);
  for (const e of entries)
    if (Number.isInteger(e.index) && e.index >= 0 && e.index < count)
      rows[e.index] = e.row;
  return rows;
}
