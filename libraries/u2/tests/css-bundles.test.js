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

/** A sheet's top-level rules, one entry per selector of a selector list; at-rule blocks are skipped. */
function rules(text) {
  const out = [];
  const src = text.replace(/\/\*[\s\S]*?\*\//g, '');
  let i = 0;
  for (let open = src.indexOf('{'); open >= 0; open = src.indexOf('{', i)) {
    const head = src.slice(i, open).trim();
    let close = src.indexOf('}', open);
    if (head.startsWith('@')) {
      for (let depth = 1, j = open + 1; depth > 0 && j < src.length; j++) {
        depth += src[j] === '{' ? 1 : src[j] === '}' ? -1 : 0;
        close = j;
      }
    }
    else {
      const decls = src.slice(open + 1, close).split(';').map((d) => d.split(':'))
        .filter((d) => d.length > 1).map(([prop, ...value]) => [prop.trim(), value.join(':').trim()]);
      for (const selector of head.split(','))
        out.push({selector: selector.replace(/\s+/g, ' ').trim(), decls});
    }
    i = close + 1;
  }
  return out;
}

test('no two sheets in all.css give one selector conflicting declarations; no rule mixes the two chip skins', () => {
  const seen = new Map();
  const conflicts = [];
  for (const sheet of all) {
    for (const {selector, decls} of rules(readFileSync(new URL(sheet, CSS), 'utf8'))) {
      assert.equal(/\.u2-chip\b(?!s)/.test(selector) && selector.includes('.u2-chips-chip'), false,
        `${sheet}: ${selector} styles an entity chip and a toggle chip at once`);
      for (const [prop, value] of decls) {
        const prev = seen.get(`${selector} ${prop}`);
        if (prev !== undefined && prev.sheet !== sheet && prev.value !== value)
          conflicts.push(`${selector} { ${prop} }: ${prev.sheet} says ${prev.value}, ${sheet} says ${value}`);
        seen.set(`${selector} ${prop}`, {sheet, value});
      }
    }
  }
  assert.ok(seen.size > 1000, 'the parser read the sheets');
  assert.deepEqual(conflicts, []);
});

test('a toggle chip is never an entity chip', async () => {
  await import('./dom-shim.js');
  const {ChipsInput} = await import('../src/components/inputs/chips-input.js');
  const input = new ChipsInput({items: ['Sales', 'Chemists'], value: ['Sales']});
  const chips = [...input.root.querySelectorAll('.u2-chips-chip')];
  assert.equal(chips.length, 2);
  assert.equal(chips.some((chip) => chip.classList.contains('u2-chip')), false);
  input.dispose();
});
