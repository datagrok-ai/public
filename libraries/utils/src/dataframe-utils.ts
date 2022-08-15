import * as DG from 'datagrok-api/dg';
import * as sha256 from 'fast-sha256';
/**
 * For columns of string type. Checks whether column contains empty values and removes corresponding rows in case user selects to remove.
 *
 */
export function removeEmptyStringRows(table: DG.DataFrame, col: DG.Column): number[] {
  const cats = col.categories;
  const emptyRawInd = cats.map((val, ind) => !val ? ind : null).filter((it) => it !== null);
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

export function hashDataFrame(table: DG.DataFrame, names?: string[]): Uint8Array {
  names ??= table.columns.names();
  const hasher = new sha256.Hash();
  const order = table.getSortedOrder(names);
  const encoder = new TextEncoder();
  for (const name of names) {
    const column = table.columns.byName(name);
    const dataArray = column.getRawData();
    const isString = column.type == DG.TYPE.STRING;
    const cats = column.categories;
    for (let i = 0; i < dataArray.length; i++) {
      if (isString) {
        const data = cats[dataArray[order[i]]];
        hasher.update(encoder.encode(data));
      } else {
        const data = dataArray[order[i]];
        hasher.update(Uint8Array.from([data]));
      }
    }
  }
  return hasher.digest();
}
