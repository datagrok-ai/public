import * as grok from 'datagrok-api/grok';
import dayjs from 'dayjs';
import {expect, expectArray, expectFloat, expectObject} from '@datagrok-libraries/utils/src/test';


export async function check(cases: {[expression: string]: any}): Promise<void> {
  for (const [expression, expected] of Object.entries(cases)) {
    const result = await grok.functions.eval(expression);
    if (Array.isArray(expected))
      expectArray(result, expected);
    else if (typeof expected === 'object' && expected instanceof dayjs)
      expect((<dayjs.Dayjs>expected).isSame(result), true);
    else if (typeof expected === 'object' && expected !== null)
      expectObject(result, expected);
    else if (typeof expected === 'number' && Number.isFinite(expected) && !Number.isInteger(expected))
      expectFloat(result, expected);
    else
      expect(result, expected);
  }
}

export async function checkRandomInt(cases: {[expression: string]: [number, number]}): Promise<void> {
  const iterationCount = 10;
  for (const [expression, [a, b]] of Object.entries(cases)) {
    for (let i = 0; i < iterationCount; i++) {
      const x = await grok.functions.eval(expression);
      expect((a < b) ? (a <= x && x < b) : (b <= x && x < a), true);
    }
  }
}
