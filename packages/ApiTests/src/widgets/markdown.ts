import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('Markdown', () => {
  test('plain text unchanged', async () => {
    const el = ui.markdown('# Hello\n\nThis is a paragraph.');
    expect(el instanceof HTMLElement, true);
    expect(el.querySelector('h1') !== null, true);
    expect(el.querySelector('.katex') === null, true);
    expect(el.querySelector('.d4-math-fallback') === null, true);
  });

  test('block math renders KaTeX or fallback', async () => {
    const el = ui.markdown('$$a^2 + b^2 = c^2$$');
    const mathEl = el.querySelector('.katex, .d4-math-fallback');
    expect(mathEl !== null, true);
  });

  test('inline math renders KaTeX or fallback', async () => {
    const el = ui.markdown('Euler: $e^{i\\pi} + 1 = 0$ is elegant.');
    const mathEl = el.querySelector('.katex, .d4-math-fallback');
    expect(mathEl !== null, true);
  });

  test('money amounts are not math', async () => {
    const el = ui.markdown('The cost is $100 and $200');
    expect(el.querySelector('.katex') === null, true);
    expect(el.querySelector('.d4-math-fallback') === null, true);
    expect(el.innerHTML.indexOf('$100') >= 0, true);
    expect(el.innerHTML.indexOf('$200') >= 0, true);
  });

  test('code span stays literal', async () => {
    const el = ui.markdown('Inline code: `$x$` should stay literal.');
    const code = el.querySelector('code');
    expect(code !== null, true);
    expect((code?.textContent ?? '').indexOf('$x$') >= 0, true);
    expect(code?.querySelector('.katex') == null, true);
    expect(code?.querySelector('.d4-math-fallback') == null, true);
  });
}, {owner: 'dkovalyov@datagrok.ai'});
