import * as DG from 'datagrok-api/dg';

/** A string column of [values], a missing value empty. */
export function textColumn(name: string, values: (string | null | undefined)[]): DG.Column {
  return DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, values.map((value) => value ?? ''));
}

export function boolColumn(name: string, values: boolean[]): DG.Column {
  return DG.Column.fromList(DG.COLUMN_TYPE.BOOL, name, values);
}
