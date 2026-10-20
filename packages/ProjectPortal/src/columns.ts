import * as DG from 'datagrok-api/dg';

/** A string column of [values]. `DG.Column.fromStrings` guesses the type from the values: a version
 * "0.3" comes back as the float 0.30000001, a code "007" as the integer 7. Ids, codes, versions and
 * free text use this instead. */
export const textColumn = (name: string, values: string[]): DG.Column =>
  DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, values);
