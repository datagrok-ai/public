/* The kind catalog (R11) against the code: every data-u2 value u2 stamps has an entry in
   src/core/kinds.ts and every entry names a stamped value, every visual registry tag builds a
   catalogued kind, and manifest.json carries what the generator makes of them. Contract:
   core/docs/features/ui2/CONTRACT.md §4. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync, readdirSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {join} from 'node:path';
import {U2_KINDS} from '../src/core/kinds.js';
import {buildManifest} from '../scripts/gen-manifest.mjs';

const SRC = fileURLToPath(new URL('../src/', import.meta.url));

/** How a data-u2 value is written in the sources. */
const STAMPS = [
  /dataset\.u2\s*=\s*'([^']+)'/g,
  /setAttribute\('data-u2',\s*'([^']+)'\)/g,
  /data-u2="([^"]+)"/g,
  /\bu2: '([^']+)'/g, // the VirtualRows chrome: list, data-table
  /_paneRow\('([^']+)'/g, // the functions browser's pane rows: fb-role, fb-tag
];

/** Bare words that would take the suffix of existing phrases ("Run history", "console output"). */
const GENERIC_WORDS = ['history', 'status', 'value', 'field', 'output', 'query'];

function sources(dir) {
  const files = [];
  for (const entry of readdirSync(dir, {withFileTypes: true})) {
    const path = join(dir, entry.name);
    if (entry.isDirectory())
      files.push(...sources(path));
    else if (entry.name.endsWith('.ts') && !entry.name.endsWith('.d.ts') && path !== join(SRC, 'core', 'kinds.ts'))
      files.push(path);
  }
  return files;
}

function stamped() {
  const values = new Set();
  for (const file of sources(SRC)) {
    const text = readFileSync(file, 'utf8');
    for (const pattern of STAMPS) {
      for (const match of text.matchAll(pattern))
        values.add(match[1]);
    }
  }
  return values;
}

const kinds = new Set(U2_KINDS.map((k) => k.kind));
const manifest = await buildManifest('../src/');

test('kinds: every data-u2 value stamped under src has an entry', () => {
  const missing = [...stamped()].filter((v) => !kinds.has(v)).sort();
  assert.deepEqual(missing, [], `add these to src/core/kinds.ts: ${missing.join(', ')}`);
});

test('kinds: every entry names a data-u2 value stamped under src', () => {
  const values = stamped();
  const stale = [...kinds].filter((k) => !values.has(k)).sort();
  assert.deepEqual(stale, [], `nothing stamps these: ${stale.join(', ')}`);
});

test('kinds: unique, owned parts resolve, nouns are no bare generic words', () => {
  assert.equal(kinds.size, U2_KINDS.length, 'a kind is catalogued once');
  for (const k of U2_KINDS) {
    if (k.partOf !== undefined) {
      const owner = U2_KINDS.find((o) => o.kind === k.partOf);
      assert.ok(owner !== undefined && owner.partOf === undefined, `${k.kind}: no owner kind ${k.partOf}`);
      continue;
    }
    for (const noun of [k.noun, ...k.aliases ?? [], k.items?.noun].filter((n) => n !== undefined))
      assert.ok(!GENERIC_WORDS.includes(noun), `${k.kind}: "${noun}" is a bare generic word`);
    if (k.family === 'input' && k.kind !== 'readonly-field') {
      assert.deepEqual(Object.keys(k.parts).slice(0, 4), ['label', 'editor', 'options', 'error'],
        `${k.kind}: an input carries the shared input parts`);
    }
  }
});

test('kinds: every visual registry tag builds a catalogued kind', () => {
  for (const c of manifest.components) {
    if (c.visual === false)
      continue;
    assert.notEqual(c.kind, undefined, `${c.tag}: its example root carries no data-u2`);
    assert.ok(kinds.has(c.kind), `${c.tag} builds "${c.kind}", which src/core/kinds.ts lacks`);
  }
  const owner = (kind) => manifest.components.find((c) => c.kind === kind)?.tag;
  for (const k of manifest.kinds)
    assert.equal(k.tag, owner(k.kind), `${k.kind}: the tag building it`);
});

test('kinds: a container\'s parts are anchored at its root, so a nested component\'s never answer (P110)', () => {
  for (const kind of ['dialog', 'section', 'card', 'wizard']) {
    for (const [name, selector] of Object.entries(U2_KINDS.find((k) => k.kind === kind).parts))
      assert.ok(selector.startsWith(':scope > '), `${kind}: part "${name}" is not anchored: ${selector}`);
  }
});

test('kinds: a partOf entry is a part of its owner', () => {
  const browser = manifest.kinds.find((k) => k.kind === 'functions-browser');
  assert.equal(browser.parts['clear button'], '[data-u2="fb-clear"]');
  assert.equal(manifest.kinds.find((k) => k.kind === 'func-form').parts['run button'], '[data-u2="ff-run"]');
});

test('kinds: manifest.json carries the generated catalog (`pnpm run manifest` regenerates it)', () => {
  const file = JSON.parse(readFileSync(new URL('../manifest.json', import.meta.url), 'utf8'));
  assert.equal(file.kindsVersion, 1);
  assert.deepEqual(file.kinds, JSON.parse(JSON.stringify(manifest.kinds)));
  const kindsOf = (components) =>
    JSON.parse(JSON.stringify(components.map((c) => ({tag: c.tag, kind: c.kind, states: c.states}))));
  assert.deepEqual(kindsOf(file.components), kindsOf(manifest.components));
});
