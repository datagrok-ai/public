import * as DG from 'datagrok-api/dg';
import * as grok from "datagrok-api/grok";
import {expect} from '@datagrok-libraries/utils/src/test';
import {runKalign} from '../utils/multiple-sequence-alignment';
import { _package} from '../package-test';

export async function loadFileAsText(name: string): Promise<string> {
  return await _package.files.readAsText(name);
}

export async function readDataframe(tableName: string): Promise<DG.DataFrame> {
  const file = await loadFileAsText(tableName);
  const df = DG.DataFrame.fromCsv(file);
  df.name = tableName.replace('.csv', '');
  return df;
}


/**
 * Tests if a table has non zero rows and columns.
 *
 * @param {DG.DataFrame} table Target table.
 */
export function _testTableIsNotEmpty(table: DG.DataFrame): void {
  expect(table.columns.length > 0 && table.rowCount > 0, true);
}


/**
 * Tests if MSA works and returns consistent result.
 *
 * @export
 * @param {DG.Column} col Aligned sequences column.
 */
export async function _testMSAIsCorrect(col: DG.Column): Promise<void> {
  const msaCol = await runKalign(col, true);
  expect(msaCol.toList().every((v, i) => (v == col.get(i) || v == null)), true);

}
