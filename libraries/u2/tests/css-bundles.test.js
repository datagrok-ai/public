/* all.css against css/: it lists every sheet but the opt-in ones. The EMS bundle's own half of this
   test lives with EMS (`libraries/ems/tests/css-bundles.test.js`), which reads both folders. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readdirSync, readFileSync} from 'node:fs';

const CSS = new URL('../css/', import.meta.url);
const OPT_IN = ['icons.css', 'theme-dark.css'];

const sheets = readdirSync(CSS).filter((f) => f.endsWith('.css') && f !== 'all.css').sort();
const listed = (text, pattern) => [...text.matchAll(pattern)].map((m) => m[1]);
const all = listed(readFileSync(new URL('all.css', CSS), 'utf8'), /@import '\.\/([\w-]+\.css)';/g);

test('all.css lists every sheet but the opt-in ones, tokens first', () => {
  assert.deepEqual([...all].sort(), sheets.filter((s) => !OPT_IN.includes(s)));
  assert.equal(all[0], 'tokens.css');
});
