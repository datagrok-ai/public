// the EMS sheet bundle (src/dg/styles.ts): every EMS sheet, and none of u2's — u2 installs its own, once per page

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readdirSync, readFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';

const EMS_CSS = new URL('../css/', import.meta.url);
const emsSheets = readdirSync(fileURLToPath(EMS_CSS)).filter((f) => f.endsWith('.css')).sort();
const styles = readFileSync(new URL('../src/dg/styles.ts', import.meta.url), 'utf8');
const taken = [...styles.matchAll(/import '([^']+\.css)';/g)].map((m) => m[1]);

test('the bundle imports no u2 sheet: the first u2 control installs u2\'s stylesheet', () => {
  assert.deepEqual(taken.filter((s) => !s.startsWith('../../css/')), []);
});

test('every sheet EMS ships is in the bundle — they have no other home', () => {
  assert.deepEqual(taken.map((s) => s.slice('../../css/'.length)).sort(), emsSheets);
});
