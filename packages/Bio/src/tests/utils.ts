import * as DG from 'datagrok-api/dg';

import {expect} from '@datagrok-libraries/utils/src/test';
import {runKalign} from '../utils/multiple-sequence-alignment';

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
