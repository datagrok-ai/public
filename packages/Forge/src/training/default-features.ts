import * as DG from 'datagrok-api/dg';

export function defaultFeatures(table: DG.DataFrame, target: DG.Column): DG.Column[] {
  return table.columns.toList().filter((c) => c.name !== target.name && isReadableNumber(c) && !isIdLike(c));
}

/** Whether the methods read [col] as numbers: numerical, not dates and not very large whole numbers (bigint). */
export function isReadableNumber(col: DG.Column): boolean {
  return col.matches(DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME) && col.type !== DG.COLUMN_TYPE.BIG_INT;
}

/** Why [method] cannot read the bigint column [name]; [ending] says what to do instead. */
export function bigIntProblem(name: string, method: string, ending: string): string {
  return `'${name}' holds very large whole numbers, which ${method} cannot read. ` +
    `Convert the column to a decimal type or ${ending}.`;
}

function isIdLike(col: DG.Column): boolean {
  if (col.type !== DG.TYPE.INT)
    return false;
  const stats = col.stats;
  return stats.missingValueCount === 0 && stats.uniqueCount === col.length;
}
