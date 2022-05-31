import * as grok from 'datagrok-api/grok';
import {expect, expectArray, expectFloat, expectObject} from '@datagrok-libraries/utils/src/test';


export async function check(cases: {[expression: string]: any}): Promise<void> {
  for (const [expression, expected] of Object.entries(cases)) {
    const result = await grok.functions.eval(expression);
    if (Array.isArray(expected))
      expectArray(result, expected);
    else if (typeof expected === 'object')
      expectObject(result, expected);
    else if (typeof expected === 'number' && !Number.isInteger(expected))
      expectFloat(result, expected);
    else
      expect(result, expected);
  }
}
