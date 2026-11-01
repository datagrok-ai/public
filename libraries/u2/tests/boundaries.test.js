// u2's boundary: no platform below src/dg except u2core, no EMS imports, no domain server API

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readdirSync, readFileSync, statSync} from 'node:fs';

const SRC = new URL('../src/', import.meta.url);

const PLATFORM_FREE = ['core', 'components', 'spec', 'sources'];

// shared with js-api so one instance serves the page
const U2CORE = 'datagrok-api/u2core';

// pinned: a new importer fails the test
const U2CORE_IMPORTERS = [
  'core/component.ts',
  'core/filter/evaluate.ts',
  'core/filter/kinds.ts',
  'core/filter/model.ts',
  'core/filter/operators.ts',
  'core/property-like.ts',
  'core/scope.ts',
  'core/signals.ts',
  'core/widget-like.ts',
  'sources/df-bindings.ts',
  'sources/rows-like.ts',
];

// guards against the regex silently matching nothing
const MIN_SPECIFIERS = 70;

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

/** Every module specifier a file names: `from '…'`, a bare `import '…'`, and `import('…')`. */
function specifiers(text) {
  return [...text.matchAll(/(?:from|import)\s*\(?\s*'([^']+)'/g)].map((m) => m[1]);
}

const FILES = sources();
// the top-level files (the index, the styles) are platform-free too
const platformFree = (path) => !path.includes('/') || PLATFORM_FREE.includes(path.split('/')[0]);
const platform = FILES.flatMap((f) =>
  specifiers(f.text).filter((s) => s.startsWith('datagrok-api')).map((s) => ({file: f.path, s})));

test('the scanner still reads the tree', () => {
  assert.ok(FILES.length > 150, `only ${FILES.length} sources found`);
  assert.ok(FILES.some((f) => f.path === 'index.ts'), 'the top-level files are read');
  assert.ok(platform.length >= MIN_SPECIFIERS,
    `only ${platform.length} platform specifiers found — the import regex stopped matching`);
});

test('the platform-free layers import the signal core and nothing else of the platform', () => {
  for (const {file, s} of platform) {
    if (platformFree(file))
      assert.equal(s, U2CORE, `${file} imports ${s}; only ${U2CORE} is allowed below src/dg`);
  }
});

test('exactly the pinned platform-free files import the signal core', () => {
  const found = [...new Set(platform.filter((p) => platformFree(p.file)).map((p) => p.file))].sort();
  assert.deepEqual(found, [...U2CORE_IMPORTERS].sort());
});

test('src/dg is the only layer that imports the platform API', () => {
  for (const {file, s} of platform) {
    if (s !== U2CORE)
      assert.ok(file.startsWith('dg/'), `${file} imports ${s} outside src/dg`);
  }
});

test('u2 does not import the EMS library', () => {
  for (const file of FILES) {
    for (const s of specifiers(file.text))
      assert.ok(!s.startsWith('@datagrok-libraries/ems'),
        `${file.path} imports ${s}; EMS depends on u2, never the other way round`);
  }
});

test('u2 does not call the domain (EMS) server API', () => {
  for (const file of FILES) {
    assert.ok(!/\bdapi\.domains\b/.test(file.text.replace(/\/\*[\s\S]*?\*\/|\/\/.*/g, '')),
      `${file.path} uses grok.dapi.domains; the domain server API belongs to @datagrok-libraries/ems`);
  }
});
