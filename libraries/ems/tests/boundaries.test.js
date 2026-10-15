// EMS's boundary: src/core is platform-free, src/dg is the only platform layer

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readdirSync, readFileSync, statSync} from 'node:fs';

const SRC = new URL('../src/', import.meta.url);

const U2CORE = 'datagrok-api/u2core';

// deep paths: the u2 dg barrel loads the whole platform layer, and authoring must not
const DEEP_U2 = {
  'dg/authoring/manifest-panel.ts': ['@datagrok-libraries/u2/src/dg/forms/object-form.js'],
  'dg/authoring/binding-dialog.ts': ['@datagrok-libraries/u2/src/dg/inputs/group-input.js'],
  'dg/authoring/binding-edit-dialog.ts': ['@datagrok-libraries/u2/src/dg/inputs/group-input.js'],
};

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

test('u2 is consumed through its entry points, bar the recorded deep paths', () => {
  for (const file of FILES) {
    for (const s of specifiers(file.text)) {
      if (!s.startsWith('@datagrok-libraries/u2/src/'))
        continue;
      const allowed = s === '@datagrok-libraries/u2/src/dg/index.js' ||
        (DEEP_U2[file.path] ?? []).includes(s);
      assert.ok(allowed, `${file.path} imports ${s}; add it to DEEP_U2 with the reason, or use an entry point`);
    }
  }
});

test('every recorded deep path is still used — the list does not rot', () => {
  for (const [file, specs] of Object.entries(DEEP_U2)) {
    const found = FILES.find((f) => f.path === file);
    assert.ok(found, `DEEP_U2 names ${file}, which is gone`);
    for (const s of specs)
      assert.ok(specifiers(found.text).includes(s), `${file} no longer imports ${s}`);
  }
});
