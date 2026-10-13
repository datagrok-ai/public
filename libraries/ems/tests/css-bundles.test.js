// the EMS sheet bundle (src/dg/styles.ts) against both libraries' css folders

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readdirSync, readFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';

const EMS_CSS = new URL('../css/', import.meta.url);
const U2_CSS = new URL('./', import.meta.resolve('@datagrok-libraries/u2/css/tokens.css'));

const NOT_IN_BUNDLE = ['accordion.css', 'adaptive.css', 'card.css', 'combobox.css', 'designer.css',
  'func-call-history-browser.css', 'func-call-input.css', 'function-input.css', 'functions-browser.css',
  'menu-bar.css', 'message-input.css', 'progress.css', 'property-grid.css', 'range-slider.css', 'spec.css',
  'table.css', 'toolbar.css', 'tour.css', 'icons.css', 'theme-dark.css', 'all.css'];

const sheets = (dir) => readdirSync(fileURLToPath(dir)).filter((f) => f.endsWith('.css')).sort();
const emsSheets = sheets(EMS_CSS);
const u2Sheets = sheets(U2_CSS);

const styles = readFileSync(new URL('../src/dg/styles.ts', import.meta.url), 'utf8');
const taken = [...styles.matchAll(/import '(?:@datagrok-libraries\/u2\/css|\.\.\/\.\.\/css)\/([\w-]+\.css)';/g)]
  .map((m) => m[1]);

test('the bundle takes every sheet it does not leave out on purpose, tokens first', () => {
  const available = [...emsSheets, ...u2Sheets];
  for (const s of taken)
    assert.ok(available.includes(s), `${s} is imported but exists in neither library`);
  assert.deepEqual([...taken].sort(),
    available.filter((s) => !NOT_IN_BUNDLE.includes(s)).sort());
  assert.equal(taken[0], 'tokens.css', 'tokens first: every other sheet reads them');
});

test('every sheet EMS ships is in the bundle — they have no other home', () => {
  assert.deepEqual(emsSheets.filter((s) => !taken.includes(s)), []);
});
