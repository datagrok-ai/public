import * as DG from 'datagrok-api/dg';
/**
 * For columns of string type. Checks whether column contains empty values and removes corresponding rows in case user selects to remove.
 *
 */
export function removeEmptyStringRows(table: DG.DataFrame, col: DG.Column): void {
    const cats = col.categories;
    const emptyRawInd = cats.map((val, ind) => !val ? ind : null).filter(it => it !== null);
    const rawData = [...col.getRawData()];
    let removedRowsCounter = 0;
    for (let i = 0; i < table.rowCount; i++) {
        if (emptyRawInd.includes(rawData[i])) {
            table.rows.removeAt(i - removedRowsCounter);
            removedRowsCounter += 1;
        }
    }
}