import * as DG from 'datagrok-api/dg';

import {runTests} from '@datagrok-libraries/utils/src/test';

import './tests/peptide-space-test';
import './tests/peptides-tests';
import './tests/msa-tests';

export const _package = new DG.Package();

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//output: dataframe result
//top-menu: Tools | Dev | JS API Tests
export async function test(category: string, test: string): Promise<DG.DataFrame> {
  const data = await runTests({category, test});
  return DG.DataFrame.fromObjects(data)!;
}
