import * as DG from 'datagrok-api/dg';

/** A non-numerical column of a semantic type in [semTypes] (lower-cased, `featureSemTypes`): a method that declares it
 * takes it as it is, so no preparation step touches it. */
export function isSemTypeFeature(col: DG.Column, semTypes: readonly string[]): boolean {
  return !col.matches(DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME) && !!col.semType &&
    semTypes.includes(col.semType.toLowerCase());
}
