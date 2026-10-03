import * as DG from 'datagrok-api/dg';

export function defaultFeatures(table: DG.DataFrame, target: DG.Column): DG.Column[] {
  return table.columns.toList().filter((c) => c.name !== target.name && c.matches('numerical') && !isIdLike(c));
}

function isIdLike(col: DG.Column): boolean {
  if (col.type !== DG.TYPE.INT && col.type !== DG.TYPE.BIG_INT)
    return false;
  const stats = col.stats;
  return stats.missingValueCount === 0 && stats.uniqueCount === col.length;
}
