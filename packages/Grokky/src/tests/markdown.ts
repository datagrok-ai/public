import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {createStyledMarkdown, tryMarkdown} from '../utils';

// The platform's markdown parser throws a RangeError on these: a paragraph ending in `](`, or in an unclosed
// link title. A model answer can end that way, and so does every streamed chunk that stops inside a link.
const UNPARSABLE = [
  'Docs: [Basic tasks → Share](',
  'Steps first.\n\nDocs: [Basic tasks → Share](\n\nMore text',
  'See [docs](https://datagrok.ai "Basic tasks',
];

category('Markdown', () => {
  test('tryMarkdown does not throw', async () => {
    for (const text of UNPARSABLE)
      tryMarkdown(text);
  });

  test('createStyledMarkdown falls back to text', async () => {
    for (const text of UNPARSABLE)
      expect(createStyledMarkdown(text).textContent!.includes('Basic tasks'), true);
  });

  test('createStyledMarkdown renders links', async () => {
    const url = 'https://datagrok.ai/help/datagrok/navigation/basic-tasks';
    const el = createStyledMarkdown(`Docs: [Basic tasks → Share](${url})`);
    expect(el.querySelector('a')?.getAttribute('href'), url);
  });
});
