import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync, readdirSync} from 'node:fs';
import {U2_KINDS} from '../src/core/kinds.js';

const css = readdirSync(new URL('../css/', import.meta.url)).filter((f) => f.endsWith('.css'))
  .map((f) => readFileSync(new URL(`../css/${f}`, import.meta.url), 'utf8')).join('\n');
const styling = [...new Set(U2_KINDS.flatMap((k) => Object.keys(k.styling ?? {})))].sort();

test('styling: every property a kind offers is one a u2 sheet reads', () => {
  assert.ok(styling.length >= 10, `${styling.length} properties`);
  assert.deepEqual(styling.filter((p) => !css.includes(`var(${p}`)), []);
});

test('styling: the README lists exactly the properties the kinds offer', () => {
  const readme = readFileSync(new URL('../README.md', import.meta.url), 'utf8');
  const section = readme.slice(readme.indexOf('\n## Styling\n'), readme.indexOf('\n## Testing\n'));
  const listed = [...new Set([...section.matchAll(/^\|([^|]*)\|/gm)].flatMap((m) =>
    [...m[1].matchAll(/`(--u2-[\w-]+)`/g)].map((p) => p[1])))].sort();
  assert.deepEqual(listed, styling);
});
