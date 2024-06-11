import * as DG from 'datagrok-api/dg';
import {TestContext, runTests, tests} from '@datagrok-libraries/utils/src/test';

import './tests/compute-api/richFunctionViewTests';
import './tests/compute-api/compositionPipelineTests';
import './tests/utils/expectTests';
import './tests/utils/jsonSerializationTests';
import './tests/compute-utils/richFunctionViewTests';
import './tests/compute-utils/compositionPipelineTests';

export const _package = new DG.Package();
export {tests};

//name: test
//input: string category {optional: true}
//input: string test {optional: true}
//input: object testContext {optional: true}
//output: dataframe result
export async function test(category: string, test: string, testContext: TestContext): Promise<DG.DataFrame> {
  const data = await runTests({category, test, testContext});
  return DG.DataFrame.fromObjects(data)!;
}
