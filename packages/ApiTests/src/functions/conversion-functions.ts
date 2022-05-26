import * as grok from 'datagrok-api/grok';
import {category, expect, expectFloat, test} from '@datagrok-libraries/utils/src/test';


category('Conversion functions', () => {
  const gfe = grok.functions.eval;

  test('Boolean', async () => {
    expect(await gfe('Boolean(true)'), true);
    expect(await gfe('Boolean("true")'), true);
    expect(await gfe('Boolean("y")'), true);
    expect(await gfe('Boolean(1)'), true);
    expect(await gfe('Boolean(10)'), true);
    expect(await gfe('Boolean(DateParse("20200131T132700"))'), true);

    expect(await gfe('Boolean(false)'), false);
    expect(await gfe('Boolean("false")'), false);
    expect(await gfe('Boolean("n")'), false);
    expect(await gfe('Boolean("abc")'), false);
    expect(await gfe('Boolean("")'), false);
    expect(await gfe('Boolean(null)'), false);
    expect(await gfe('Boolean(0)'), false);
  });

  test('ParseFloat', async () => {
    expect(await gfe('ParseFloat("2025")'), 2025);
    expectFloat(await gfe('ParseFloat("12.78")'), 12.78);
    expectFloat(await gfe('ParseFloat("-012.150")'), -12.15);
  });

  test('ParseInt', async () => {
    expect(await gfe('ParseInt("2025")'), 2025);
    expect(await gfe('ParseInt("-012")'), -12);
    expect(await gfe('ParseInt(" 0101 ")'), 101);
  });

  test('ToString', async () => {
    expect(await gfe('ToString(1)'), '1');
    expect(await gfe('ToString(3.14)'), '3.14');
    expect(await gfe('ToString(true)'), 'true');
  });
});
