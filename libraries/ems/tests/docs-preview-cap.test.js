// the CRUD recipe states how many rows the import preview checks: it must be the number the wizard checks

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';

const recipe = readFileSync(new URL('../docs/recipes/crud-app.md', import.meta.url), 'utf8').replace(/\s+/g, ' ');
const code = readFileSync(new URL('../src/dg/import.ts', import.meta.url), 'utf8');

test('the CRUD recipe states the import preview cap the code uses', () => {
  const rows = Number(/const PREVIEW_ROWS = (\d+);/.exec(code)?.[1]);
  assert.ok(rows > 0, 'PREVIEW_ROWS not found in src/dg/import.ts');
  const stated = /the first ([\d,]+k?) rows/.exec(recipe)?.[1];
  assert.ok(stated, 'the recipe states no preview cap');
  const value = stated.endsWith('k') ? Number(stated.slice(0, -1)) * 1000 : Number(stated.replace(/,/g, ''));
  assert.equal(value, rows, `the recipe says the preview checks the first ${stated} rows; PREVIEW_ROWS is ${rows}`);
});
