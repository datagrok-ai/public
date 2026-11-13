// EMS's boundary: src/core is platform-free, src/dg is the only platform layer

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readdirSync, readFileSync, statSync} from 'node:fs';

const SRC = new URL('../src/', import.meta.url);

const U2CORE = 'datagrok-api/u2core';

const U2_ENTRIES = ['@datagrok-libraries/u2', '@datagrok-libraries/u2/dg', '@datagrok-libraries/u2/filter',
  '@datagrok-libraries/u2/dg/filter'];

function sources(dir = SRC, prefix = '') {
  const out = [];
  for (const entry of readdirSync(dir)) {
    const url = new URL(entry, dir);
    if (statSync(url).isDirectory()) {
      out.push(...sources(new URL(`${entry}/`, dir), `${prefix}${entry}/`));
      continue;
    }
    if (entry.endsWith('.ts') && !entry.endsWith('.d.ts'))
      out.push({path: `${prefix}${entry}`, text: readFileSync(url, 'utf8')});
  }
  return out;
}

function specifiers(text) {
  return [...text.matchAll(/(?:from|import)\s*\(?\s*'([^']+)'/g)].map((m) => m[1]);
}

const FILES = sources();

test('the scanner still reads the tree', () => {
  assert.ok(FILES.length > 25, `only ${FILES.length} sources found`);
});

test('src/core runs without a platform: only the signal core comes from datagrok-api', () => {
  for (const file of FILES) {
    if (!file.path.startsWith('core/'))
      continue;
    for (const s of specifiers(file.text)) {
      if (s.startsWith('datagrok-api'))
        assert.equal(s, U2CORE, `core/${file.path} imports ${s}`);
    }
  }
});

test('src/dg is the only layer that imports the platform API', () => {
  for (const file of FILES) {
    for (const s of specifiers(file.text)) {
      if (s.startsWith('datagrok-api') && s !== U2CORE)
        assert.ok(file.path.startsWith('dg/'), `${file.path} imports ${s} outside src/dg`);
    }
  }
});

test('u2 is consumed through its entry points', () => {
  for (const file of FILES) {
    for (const s of specifiers(file.text)) {
      if (s.startsWith('@datagrok-libraries/u2'))
        assert.ok(U2_ENTRIES.includes(s), `${file.path} imports ${s}; use an entry point`);
    }
  }
});
