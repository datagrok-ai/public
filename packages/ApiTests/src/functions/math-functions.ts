import * as grok from 'datagrok-api/grok';
import { category, expect, test } from '@datagrok-libraries/utils/src/test';


category('Math functions', () => {
  const gfe = grok.functions.eval;

  test('Add', async () => {
    expect(await gfe('Add(2, 10)'), 12);
  });
});
