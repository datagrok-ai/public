/* The two sheet bundles against css/: all.css lists every sheet but the opt-in ones, and the
   domain bundle (src/dg/domain/styles.ts) takes each sheet or leaves it out on purpose — a new
   sheet fails here until it is placed. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readdirSync, readFileSync} from 'node:fs';

const CSS = new URL('../css/', import.meta.url);
const OPT_IN = ['icons.css', 'theme-dark.css'];
const NOT_IN_DOMAIN = ['accordion.css', 'adaptive.css', 'card.css', 'combobox.css', 'designer.css',
  'func-call-history-browser.css', 'func-call-input.css', 'function-input.css', 'functions-browser.css',
  'menu-bar.css', 'message-input.css', 'progress.css', 'property-grid.css', 'range-slider.css', 'spec.css',
  'table.css', 'toolbar.css', 'tour.css', ...OPT_IN];

const sheets = readdirSync(CSS).filter((f) => f.endsWith('.css') && f !== 'all.css').sort();
const listed = (text, pattern) => [...text.matchAll(pattern)].map((m) => m[1]);
const all = listed(readFileSync(new URL('all.css', CSS), 'utf8'), /@import '\.\/([\w-]+\.css)';/g);
const domain = listed(readFileSync(new URL('../src/dg/domain/styles.ts', import.meta.url), 'utf8'),
  /import '\.\.\/\.\.\/\.\.\/css\/([\w-]+\.css)';/g);

test('all.css lists every sheet but the opt-in ones, tokens first', () => {
  assert.deepEqual([...all].sort(), sheets.filter((s) => !OPT_IN.includes(s)));
  assert.equal(all[0], 'tokens.css');
});

test('the domain bundle takes every sheet it does not leave out on purpose, tokens first', () => {
  for (const s of domain)
    assert.ok(sheets.includes(s), `${s} is imported but does not exist`);
  assert.deepEqual([...domain].sort(), sheets.filter((s) => !NOT_IN_DOMAIN.includes(s)));
  assert.equal(domain[0], 'tokens.css');
});
