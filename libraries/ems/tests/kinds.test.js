/* EMS's kind catalog (R11) against the code: every data-u2 value EMS stamps has an entry in
   src/core/kinds.ts and every entry names a stamped value, no kind collides with u2's, and
   manifest.json carries what the generator makes of them. Contract:
   core/docs/features/ui2/overhaul/CONTRACT.md §4. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync, readdirSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {join} from 'node:path';
import {U2_KINDS} from '@datagrok-libraries/u2/src/core/kinds.js';
import {EMS_KINDS} from '../src/core/kinds.js';
import {buildManifest} from '../scripts/gen-manifest.mjs';

const SRC = fileURLToPath(new URL('../src/', import.meta.url));

const STAMPS = [
  /dataset\.u2\s*=\s*'([^']+)'/g,
  /setAttribute\('data-u2',\s*'([^']+)'\)/g,
  /data-u2="([^"]+)"/g,
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

const kinds = new Set(EMS_KINDS.map((k) => k.kind));
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

test('kinds: unique here and against u2, nouns are no bare generic words', () => {
  assert.equal(kinds.size, EMS_KINDS.length, 'a kind is catalogued once');
  const shared = U2_KINDS.filter((k) => kinds.has(k.kind)).map((k) => k.kind);
  assert.deepEqual(shared, [], 'u2 catalogues these already');
  for (const k of EMS_KINDS) {
    for (const noun of [k.noun, ...k.aliases ?? [], k.items?.noun].filter((n) => n !== undefined))
      assert.ok(!GENERIC_WORDS.includes(noun), `${k.kind}: "${noun}" is a bare generic word`);
  }
});

test('kinds: every visual registry tag builds a catalogued kind', () => {
  for (const c of manifest.components.filter((c) => c.visual !== false))
    assert.ok(kinds.has(c.kind), `${c.tag} builds "${c.kind}", which src/core/kinds.ts lacks`);
});

test('kinds: manifest.json carries the generated catalog (`pnpm run manifest` regenerates it)', () => {
  const file = JSON.parse(readFileSync(new URL('../manifest.json', import.meta.url), 'utf8'));
  assert.equal(file.kindsVersion, 1);
  assert.deepEqual(file.kinds, JSON.parse(JSON.stringify(manifest.kinds)));
  assert.deepEqual(file.components.map((c) => c.kind), manifest.components.map((c) => c.kind));
});
