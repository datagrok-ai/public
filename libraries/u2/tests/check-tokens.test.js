import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {checkStates, MIN_DELTA_E} from '../scripts/check-tokens.mjs';

const root = fileURLToPath(new URL('..', import.meta.url));
const read = (file) => readFileSync(new URL(`../${file}`, import.meta.url), 'utf8');
const edited = (file, from, to) => (f) => {
  const text = read(f);
  if (f !== file)
    return text;
  assert.ok(text.includes(from), `${file} no longer has ${from}`);
  return text.replace(from, to);
};

test('lint:tokens passes on the sheets as they are, row states included', () => {
  const r = spawnSync(process.execPath, ['scripts/check-tokens.mjs'], {cwd: root, encoding: 'utf8'});
  assert.equal(r.status, 0, r.stderr);
  assert.match(r.stdout, /row states apart, closest \w+\/\w+ ΔE [\d.]+/);
});

test('the row states read apart standalone, under the dark theme and with the platform palette', () => {
  const {problems, closest} = checkStates();
  assert.deepEqual(problems, []);
  assert.ok(closest.d >= MIN_DELTA_E, JSON.stringify(closest));
});

test('mutation: a hover painted as the section header fails (B-6)', () => {
  const {problems} = checkStates(edited('css/tokens.css',
    '--u2-row-hover-bg: var(--dg-grid-mouse-over-row);', '--u2-row-hover-bg: var(--dg-hover-subtle);'));
  assert.ok(problems.length > 0);
  for (const p of problems)
    assert.match(p, /: the hover row \(--u2-row-hover-bg, rgb\(\d+, \d+, \d+\)\) and the header \(--u2-section-header-bg, /);
  assert.ok(problems.some((p) => p.startsWith('standalone: ') && p.includes('ΔE 0.00')));
});

test('mutation: a header too close to the hover to tell apart fails', () => {
  const {problems} = checkStates(edited('css/tokens.css',
    '--u2-section-header-bg: var(--dg-bg-secondary);', '--u2-section-header-bg: rgb(240, 246, 251);'));
  assert.ok(problems.some((p) => p.startsWith('standalone: the hover row') && p.includes('the header')), problems.join('\n'));
});

test('mutation: the dark theme is checked on its own', () => {
  const {problems} = checkStates(edited('css/theme-dark.css', '--dg-selection-bg: #274D6E;', '--dg-selection-bg: #26282C;'));
  assert.ok(problems.length > 0);
  assert.ok(problems.every((p) => p.startsWith('dark: ') && p.includes('selected')), problems.join('\n'));
});

test('mutation: a state token missing from tokens.css fails', () => {
  const {problems} = checkStates(edited('css/tokens.css', '--u2-row-current-bg: var(--dg-grid-current-row);', ''));
  assert.deepEqual(problems, ['--u2-row-current-bg (the current state) is not declared in css/tokens.css']);
});
