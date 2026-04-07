import type * as _DG from 'datagrok-api/dg';
/**
 * For columns of string type. Checks whether column contains empty values and removes corresponding rows in case user selects to remove.
 *
 */
export declare function removeEmptyStringRows(table: _DG.DataFrame, col: _DG.Column): number[];
export declare function hashDataFrame(table: _DG.DataFrame, names?: string[]): Uint8Array;
export declare const testData: _DG.DataFrame;
//# sourceMappingURL=dataframe-utils.d.ts.map