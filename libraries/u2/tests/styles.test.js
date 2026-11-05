/* One u2 stylesheet per page, whatever number of package copies load and in whatever order: the
   newest revision's text wins, an equal or older copy changes nothing, and a host that links the
   sheets itself opts out. The generated sheet is css/all.css inlined, and its lock is current. */

import {test} from './owned-test.js';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import './dom-shim.js';
import {Styles} from '../src/styles.js';
import {U2_STYLES} from '../src/styles.g.js';
import {TextInput} from '../src/index.js';

const sheets = () => document.head.querySelectorAll('style[data-u2-styles]');

function fresh() {
  for (const el of sheets())
    el.remove();
  for (const el of document.head.querySelectorAll('meta[name="u2-styles"]'))
    el.remove();
}

test('loading u2 installs exactly one stylesheet, at this copy\'s revision', () => {
  new TextInput({label: 'x'}).dispose();
  assert.equal(sheets().length, 1);
  assert.deepEqual({revision: Styles.installed().revision, hash: Styles.installed().hash},
    {revision: Styles.revision, hash: Styles.hash});
  assert.ok(sheets()[0].textContent.includes('.u2-btn'));
});

test('an older build\'s own copy of the u2 sheet is reported once and left in place', () => {
  fresh();
  const legacy = document.createElement('style');
  legacy.textContent = readFileSync(new URL('../css/inputs.css', import.meta.url), 'utf8');
  const override = document.createElement('style');
  override.textContent = '.my-pkg .u2-input-root{display:block}';
  document.head.appendChild(legacy);
  document.head.appendChild(override);
  const warnings = [];
  const warn = console.warn;
  console.warn = (m) => warnings.push(m);
  try {
    assert.equal(Styles.legacyCopies(document), 1, 'a package that only restyles u2 is no legacy copy');
    Styles.install('.c{}', 9, 'h9');
    Styles.install('.c{}', 9, 'h9');
  }
  finally {
    console.warn = warn;
  }
  assert.equal(warnings.filter((w) => w.includes('before stylesheet revisions')).length, 1);
  assert.ok(legacy.isConnected, 'the old copy is reported, not removed');
  legacy.remove();
  override.remove();
});

test('a newer copy replaces the text in place; an older or identical one leaves it', () => {
  fresh();
  Styles.install('.a{}', 5, 'h5');
  const el = Styles.installed().element;
  Styles.install('.b{}', 7, 'h7');
  assert.equal(sheets().length, 1);
  assert.equal(Styles.installed().element, el, 'the same element keeps its place in the head');
  assert.equal(el.textContent, '.b{}');
  Styles.install('.old{}', 6, 'h6');
  Styles.install('.b{}', 7, 'h7');
  assert.equal(el.textContent, '.b{}');
  assert.equal(Styles.installed().revision, 7);
});

test('two different copies of one revision keep the first and warn once', () => {
  fresh();
  const warnings = [];
  const warn = console.warn;
  console.warn = (...a) => warnings.push(a.join(' '));
  try {
    Styles.install('.first{}', 9, 'aaa');
    Styles.install('.second{}', 9, 'bbb');
    Styles.install('.third{}', 9, 'ccc');
  }
  finally {
    console.warn = warn;
  }
  assert.equal(Styles.installed().element.textContent, '.first{}');
  assert.equal(warnings.length, 1);
  assert.match(warnings[0], /revision 9 differ/);
});

test('a host that links the sheets itself opts out', () => {
  fresh();
  const meta = document.createElement('meta');
  meta.setAttribute('name', 'u2-styles');
  meta.setAttribute('content', 'linked');
  document.head.append(meta);
  assert.equal(Styles.install(), null);
  assert.equal(sheets().length, 0);
  fresh();
  Styles.install();
});

test('the generated sheet is every sheet all.css lists, and css/styles.lock.json is current', () => {
  const all = readFileSync(new URL('../css/all.css', import.meta.url), 'utf8');
  for (const [, sheet] of all.matchAll(/@import '\.\/([\w-]+)\.css';/g)) {
    const first = readFileSync(new URL(`../css/${sheet}.css`, import.meta.url), 'utf8')
      .replace(/\/\*[\s\S]*?\*\//g, '').split(/\r?\n/).map((l) => l.trim()).find((l) => l);
    if (first)
      assert.ok(U2_STYLES.css.includes(first), `${sheet}.css is missing from the generated sheet`);
  }
  const lock = JSON.parse(readFileSync(new URL('../css/styles.lock.json', import.meta.url), 'utf8'));
  assert.deepEqual({revision: lock.revision, hash: lock.hash}, {revision: U2_STYLES.revision, hash: U2_STYLES.hash});
  assert.equal(createHash('sha256').update(U2_STYLES.css).digest('hex').slice(0, 12), U2_STYLES.hash);
});

test('the sheet goes after the platform\'s links and before the packages\' styles, whatever loaded first', () => {
  fresh();
  const platform = document.createElement('link');
  platform.setAttribute('rel', 'stylesheet');
  const packageStyle = document.createElement('style');
  packageStyle.textContent = '.my-package .u2-btn{}';
  document.head.append(platform, packageStyle);
  try {
    const el = Styles.install('.u2-btn{}', 99, 'h99').element;
    const order = [...document.head.children];
    assert.ok(order.indexOf(platform) < order.indexOf(el), 'after the platform');
    assert.ok(order.indexOf(el) < order.indexOf(packageStyle), 'before a package style that loaded earlier');
  }
  finally {
    platform.remove();
    packageStyle.remove();
    fresh();
    Styles.install();
  }
});
