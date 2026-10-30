/* Layout helpers space their children on the token scale, so nobody hand-writes `style.gap`: `divH(children,
   {gap: 'm'})` is a row 8px apart, the `u2-gap-<size>` classes work on any flex or grid element. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import './dom-shim.js';
import {div, divH, divV, span} from '../src/index.js';

const css = (name) => readFileSync(new URL(`../css/${name}`, import.meta.url), 'utf8');

test('divH / divV / div take a gap on the token scale and an extra class', () => {
  const row = divH([span('a'), span('b')], {gap: 'm', cls: 'my-row'});
  assert.deepEqual(row.className.split(' '), ['u2-div', 'u2-div-h', 'u2-gap-m', 'my-row']);
  assert.equal(row.children.length, 2);
  assert.equal(divV([], {gap: 'xl'}).classList.contains('u2-gap-xl'), true);
  assert.equal(div([], {gap: 's'}).classList.contains('u2-gap-s'), true);
  assert.equal(divH([], 'legacy-class').className, 'u2-div u2-div-h legacy-class', 'a class string still works');
  assert.equal(divH([]).style.getPropertyValue('gap'), '', 'never an inline gap');
});

test('a form\'s button row centres what stands beside its buttons (a statusLine "Saved")', () => {
  assert.match(css('form.css'), /\.u2-form-buttons \{[^}]*align-items: center;/);
});

test('the DOM shim: firstElementChild, and a short inspect so a failing assert on elements does not hang', async () => {
  const {inspect} = await import('node:util');
  const host = div([span('a'), span('b')]);
  assert.ok(host.firstElementChild === host.children[0]);
  assert.equal(div([]).firstElementChild, null);
  for (let i = 0; i < 200; i++)
    host.append(div([span(String(i))]));
  document.body.append(host);
  const started = Date.now();
  const text = inspect(host, {depth: 10});
  assert.ok(Date.now() - started < 200, 'inspect returns at once');
  assert.match(text, /^<div class="u2-div">/);
  host.remove();
});

test('every gap class is the token of its size', () => {
  const sheet = css('elements.css');
  for (const size of ['xs', 's', 'm', 'l', 'xl', 'xxl'])
    assert.match(sheet, new RegExp(`\\.u2-gap-${size} \\{\\s*gap: var\\(--dg-space-${size}\\);`));
});
