import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Widgets: Markdown', () => {
  test('Math block renders via KaTeX or fallback', async () => {
    const el = ui.markdown(String.raw`# Heading

$$\frac{a}{b}$$`);
    const html = el.innerHTML;
    expect(html.includes('katex') || html.includes('d4-math-fallback'), true);
  });

  test('Money amounts are not treated as inline math', async () => {
    const el = ui.markdown('The cost is $100 and $200');
    expect(el.innerHTML.includes('katex'), false);
  });

  test('Plain markdown without math renders normally', async () => {
    const el = ui.markdown('# Plain');
    expect(/<h1[^>]*>Plain<\/h1>/.test(el.innerHTML), true);
  });
}, {owner: 'dkovalyov@datagrok.ai'});
