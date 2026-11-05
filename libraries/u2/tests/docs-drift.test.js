/* The reference docs name only APIs that exist. A code span `Class.member` whose class u2 exports
   must name a member of it, and a bare call `name()` must be something u2 (or EMS) exports or a
   member of an exported class — `host()` stayed in docs for weeks after it was deleted. The
   reference is this library's README and recipes, and the two shipped contracts in the core repo
   (checked when core/ sits next to public/); the rest of core/docs/features/ui2 is design history. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {existsSync, readdirSync, readFileSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {join} from 'node:path';
import './dom-shim.js';
import * as u2 from '../src/index.js';
import * as core from 'datagrok-api/u2core';
import {U2_RULES} from '../src/core/rules.js';

const root = fileURLToPath(new URL('..', import.meta.url));
const CORE_UI2 = join(root, '..', '..', '..', 'core', 'docs', 'features', 'ui2');
const CORE_DOCS = ['CONTRACT.md', 'AUTOMATION.md'];

/** A doc's text less the rules' wrong lines, which show on purpose what does not exist or must not be written. */
const WRONG = U2_RULES.map((r) => `\`${r.wrong}\``);
const readDoc = (file) => WRONG.reduce((text, wrong) => text.split(wrong).join('``'), readFileSync(file, 'utf8'));

function docs() {
  const out = [join(root, 'README.md'), join(root, 'AGENTS.md'), join(root, 'manifest.json')];
  for (const dir of ['recipes', 'samples'])
    for (const f of readdirSync(join(root, 'docs', dir)))
      if (f.endsWith('.md'))
        out.push(join(root, 'docs', dir, f));
  if (existsSync(CORE_UI2))
    out.push(...CORE_DOCS.map((f) => join(CORE_UI2, f)).filter(existsSync));
  return out;
}

/** Every name declared in u2's sources and EMS's — the platform layer needs the platform to import — and
 * per class, its base and the members its body declares (fields included, which no prototype shows). */
const classes = new Map();
function declared() {
  const names = new Set();
  const walk = (dir) => {
    for (const e of readdirSync(dir, {withFileTypes: true})) {
      const p = join(dir, e.name);
      if (e.isDirectory())
        walk(p);
      else if (e.name.endsWith('.ts') && !e.name.endsWith('.d.ts')) {
        const text = readFileSync(p, 'utf8');
        for (const m of text.matchAll(/export\s+(?:abstract\s+)?(?:class|function|const|let|interface|type|enum)\s+(\w+)/g))
          names.add(m[1]);
        for (const m of text.matchAll(/export\s+(?:type\s+)?\{([^}]+)\}/g))
          for (const part of m[1].split(','))
            names.add(part.trim().split(/\s+as\s+/).pop().trim());
        for (const m of text.matchAll(/^\s+(?:static\s+|async\s+|get\s+|set\s+|protected\s+|private\s+|readonly\s+)*(\w+)\s*[(<:=]/gm))
          names.add(m[1]);
        const heads = [...text.matchAll(/^export\s+(?:abstract\s+)?class\s+(\w+)([^{]*)\{/gm)];
        heads.forEach((m, i) => {
          // the base is the `extends` left once the type parameters are gone (they carry their own)
          let head = m[2];
          while (/<[^<>]*>/.test(head))
            head = head.replace(/<[^<>]*>/g, '');
          const base = /\bextends\s+(\w+)/.exec(head)?.[1];
          const body = text.slice(m.index + m[0].length, heads[i + 1]?.index ?? text.length);
          const members = new Set();
          for (const f of body.matchAll(/^  (?:static |async |get |set |protected |private |readonly |override |abstract )*(\w+)\s*[?!]?\s*[(<:=;]/gm))
            members.add(f[1]);
          classes.set(m[1], {base, members});
        });
      }
    }
  };
  walk(join(root, 'src'));
  walk(join(root, '..', '..', 'js-api', 'src', 'u2core'));
  if (existsSync(join(root, '..', 'ems', 'src')))
    walk(join(root, '..', 'ems', 'src'));
  return names;
}

function membersOf(x) {
  const names = new Set();
  for (let o = x; o && o !== Function.prototype && o !== Object.prototype; o = Object.getPrototypeOf(o))
    for (const k of Object.getOwnPropertyNames(o))
      names.add(k);
  if (typeof x === 'function')
    for (let o = x.prototype; o && o !== Object.prototype; o = Object.getPrototypeOf(o))
      for (const k of Object.getOwnPropertyNames(o))
        names.add(k);
  return names;
}

const exported = {...core, ...u2};
const known = new Set([...Object.keys(exported), ...declared()]);
for (const value of Object.values(exported))
  if (typeof value === 'function' || (typeof value === 'object' && value !== null))
    for (const k of membersOf(value))
      known.add(k);

function declares(cls, name) {
  for (let c = classes.get(cls); c !== undefined; c = classes.get(c.base))
    if (c.members.has(name))
      return true;
  return false;
}

/** Bare calls the docs make to the language, the platform and the test tools — not u2 API. */
const ELSEWHERE = new Set(['fetch', 'setTimeout', 'requestAnimationFrame', 'queueMicrotask', 'require', 'import',
  'expect', 'test', 'describe', 'it', 'await', 'grok', 'ui', 'DG', 'console', 'JSON', 'Promise', 'Math', 'String',
  'Number', 'Array', 'Object', 'Date', 'Error', 'Map', 'Set', 'WeakMap', 'Symbol', 'getComputedStyle',
  'querySelector', 'querySelectorAll', 'addEventListener', 'removeEventListener', 'closest', 'contains', 'append',
  'appendChild', 'remove', 'click', 'focus', 'blur', 'scrollIntoView', 'getBoundingClientRect', 'evaluate',
  'locator', 'waitFor', 'fill', 'press', 'goto', 'screenshot', 'step', 'stepName', 'toDate']);

/** Where a package path lands under `exports`, as a source file: `./src/dg/index.js` → `src/dg/index.ts`. */
function exportSource(subpath) {
  const exports = JSON.parse(readFileSync(join(root, 'package.json'), 'utf8')).exports;
  for (const [key, entry] of Object.entries(exports)) {
    const target = typeof entry === 'string' ? entry : entry.default;
    const star = key.indexOf('*');
    const match = star < 0 ? (key === subpath ? '' : null) :
      subpath.startsWith(key.slice(0, star)) && subpath.endsWith(key.slice(star + 1)) ?
        subpath.slice(star, subpath.length - (key.length - star - 1)) : null;
    if (match !== null)
      return target.replace('*', match).replace(/^\.\/dist\//, './').replace(/\.js$/, '.ts');
  }
  return null;
}

test('every u2 import path the living docs show resolves', () => {
  const broken = [];
  const files = [...docs(), join(root, 'src', 'dg', 'README.md')];
  for (const file of files) {
    const text = readDoc(file);
    for (const [, path] of text.matchAll(/from '(@datagrok-libraries\/u2[^']*)'/g)) {
      const source = exportSource(`.${path.slice('@datagrok-libraries/u2'.length)}`);
      if (source === null || !existsSync(join(root, source)))
        broken.push(`${file.slice(root.length - 1)}: ${path}`);
    }
  }
  assert.deepEqual(broken, []);
});

test('every API the living docs name exists', () => {
  const missing = [];
  for (const file of docs()) {
    const text = readDoc(file);
    for (const [, span] of text.matchAll(/`([^`\n]+)`/g)) {
      const member = /^([A-Z]\w*)\.(\w+)(?:\(.*\))?$/.exec(span);
      if (member !== null && exported[member[1]] !== undefined) {
        if (!membersOf(exported[member[1]]).has(member[2]) && !declares(member[1], member[2]))
          missing.push(`${file.slice(root.length - 1)}: \`${span}\` — ${member[1]} has no ${member[2]}`);
        continue;
      }
      const call = /^([a-z]\w*)\((?:[^()]*)\)$/.exec(span);
      if (call !== null && !known.has(call[1]) && !ELSEWHERE.has(call[1]))
        missing.push(`${file.slice(root.length - 1)}: \`${span}\` — no u2 export or member is called ${call[1]}`);
    }
  }
  assert.deepEqual(missing, []);
});

/** Claims u2's own docs once made that no longer hold, each with what holds instead. */
const STALE_CLAIMS = [
  [/custom-element tags|framework wrappers|fluent `u2\.\*`/, 'u2 has no custom elements, framework wrappers or fluent API'],
  [/enforced by eslint|eslint enforces/i, 'the layer boundary is held by tests/boundaries.test.js'],
  [/^\| `signals-core\/`/m, 'the signal core is the engine\'s copy (js-api/src/u2core/signals-core.ts)'],
  [/\b(?:appView|ownView|viewSignals|showInView|tableViewChrome)\b/, 'gone: a view is owned through ViewHost'],
  [/\bokEnabled\(|\brunInScope\(|\bControl\.build\(/, 'gone: okDisabled(reason), owner.build(fn), component(fn)'],
];

test('u2\'s own docs make none of the stale claims', () => {
  const files = [...docs().filter((f) => f.startsWith(root)), join(root, 'src', 'dg', 'README.md'),
    join(root, 'vendor', 'LICENSES.md')];
  const found = [];
  for (const file of files)
    for (const [claim, now] of STALE_CLAIMS)
      if (claim.test(readDoc(file)))
        found.push(`${file.slice(root.length - 1)}: ${claim} — ${now}`);
  assert.deepEqual(found, []);
});

test('mutation: a stale claim put back is found', () => {
  assert.ok(STALE_CLAIMS.some(([claim]) => claim.test('- **L3 delivery** — fluent `u2.*` API, custom-element tags')));
  assert.ok(STALE_CLAIMS.some(([claim]) => claim.test('A view is owned through `ownView`.')));
});
