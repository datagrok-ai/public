import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';


category('Logical functions', () => {
  const gfe = grok.functions.eval;

  test('And', async () => {
    expect(await gfe('And(true, true)'), true);
    expect(await gfe('And(true, false)'), false);
    expect(await gfe('And(false, true)'), false);
    expect(await gfe('And(false, false)'), false);
    expect(await gfe('And(5 == 5, 10 < 20)'), true);
    expect(await gfe('And(2, 3)'), 3);
    expect(await gfe('And(0, 1)'), 0);
    expect(await gfe('And(1, 1)'), 1);
  });

  test('Not', async () => {
    expect(await gfe('Not(true)'), false);
    expect(await gfe('Not(false)'), true);
    expect(await gfe('Not(1)'), false);
    expect(await gfe('Not(0)'), true);
  });

  test('Or', async () => {
    expect(await gfe('Or(true, true)'), true);
    expect(await gfe('Or(true, false)'), true);
    expect(await gfe('Or(false, true)'), true);
    expect(await gfe('Or(false, false)'), false);
    expect(await gfe('Or(5 == 6, 20 < 10)'), false);
    expect(await gfe('Or(2, 3)'), 2);
    expect(await gfe('Or(0, 1)'), 1);
    expect(await gfe('Or(1, 1)'), 1);
  });


  test('Xor', async () => {
    expect(await gfe('Xor(true, true)'), false);
    expect(await gfe('Xor(true, false)'), true);
    expect(await gfe('Xor(false, true)'), true);
    expect(await gfe('Xor(false, false)'), false);
    expect(await gfe('Xor(5 == 6, 20 < 10)'), false);
    expect(await gfe('Xor(5 == 5, 10 < 20)'), false);
    expect(await gfe('Xor(2, 3)'), false);
    expect(await gfe('Xor(2, 2)'), false);
    expect(await gfe('Xor(1, 0)'), true);
    expect(await gfe('Xor(1, 1)'), false);
  });
});
