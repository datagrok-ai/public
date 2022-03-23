import * as DG from 'datagrok-api/dg';
import {runTests} from '@datagrok-libraries/utils/src/test';
import './tests/smiles-tests';

export const _package = new DG.Package();

//name: test
//output: dataframe result
export async function test(): Promise<DG.DataFrame> {
  const data = await runTests();
  return DG.DataFrame.fromObjects(data)!;
}
