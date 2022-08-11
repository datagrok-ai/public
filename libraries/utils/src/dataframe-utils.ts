import * as DG from 'datagrok-api/dg';
/**
 * For columns of string type. Checks whether column contains empty values and removes corresponding rows in case user selects to remove.
 *
 */
export function removeEmptyStringRows(table: DG.DataFrame, col: DG.Column): number[] {
    const cats = col.categories;
    const emptyRawInd = cats.map((val, ind) => !val ? ind : null).filter(it => it !== null);
    const rawData = [...col.getRawData()];
    const emptyRawsIndexes = [];
    let removedRowsCounter = 0;
    for (let i = 0; i < table.rowCount; i++) {
        if (emptyRawInd.includes(rawData[i])) {
            table.rows.removeAt(i - removedRowsCounter);
            emptyRawsIndexes.push(i);
            removedRowsCounter += 1;
        }
    }
    return emptyRawsIndexes;
}