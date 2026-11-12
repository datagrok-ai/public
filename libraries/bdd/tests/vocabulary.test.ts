/* The u2 vocabulary (R11): the generator writes the manifests' kinds deterministically, kindDefsFrom
   turns them into one kind per noun merged with the bdd overrides, and compareKinds proves a kind table
   keeps every key, selector, strategy, part and Dart name of the one taken before the rewrite
   (fixtures/kinds-before.json, tool/snapshot-kinds.mjs). The registered table is the generated one:
   it equals the manifests (drift), names every manifest kind and no other (coverage), and keeps the
   snapshot (compat). */
import assert from 'node:assert/strict';
import {spawnSync} from 'node:child_process';
import {existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {dirname, join} from 'node:path';
import {after, before, test} from 'node:test';
import {fileURLToPath, pathToFileURL} from 'node:url';
import {isDeepStrictEqual} from 'node:util';
import {type Browser, chromium, type Page} from '@playwright/test';

import '../bindings/common/kinds.js';
import {U2_KINDS} from '../bindings/common/u2-kinds.g.js';
import {allKinds, type KindDef, lookupKind, type MatchStrategy} from '../src/registry.js';
import {el} from '../src/runtime/args.js';
import {describeCandidates, effectiveOrder, locate} from '../src/runtime/locate.js';
import {kindDefsFrom, registerU2Kinds, selectorList, type U2KindSpec} from '../src/u2-vocabulary.js';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const FIXTURE = join(ROOT, 'tests/fixtures/manifest-kinds.json');
const BEFORE: Record<string, KindDef & {name: string}> = JSON.parse(readFileSync(join(ROOT, 'tests/fixtures/kinds-before.json'), 'utf8'));
const tmp = mkdtempSync(join(tmpdir(), 'bdd-vocab-'));
after(() => rmSync(tmp, {recursive: true, force: true}));

/** Kinds that only lay out their children: no noun has to name them (CONTRACT §5.3). */
const LAYOUT_ONLY = ['panel', 'div-v', 'div-h'];
const ITEM_KINDS = ['item', 'tree node', 'table row', 'tab', 'wizard step', 'menu item', 'breadcrumb', 'chip'];

/** Selector alternatives of the snapshot dropped on purpose, by kind: the alternative that replaced
 * each, or null when nothing replaces it. */
const REPLACED: Record<string, Record<string, string | null>> = {
  // nothing stamps data-u2="item" (u2, EMS, U2Demo); u2 rows are .u2-list-row / [role=option] with data-u2-item
  'item': {'[data-u2="item"]': null},
  // the Dart root inside a u2 dart-input is the same field, not a second input (R9 ambiguity)
  'input': {'.ui-input-root': '.ui-input-root:not([data-u2="dart-input"] .ui-input-root)'},
};

const generate = (args: string[], cwd = ROOT) =>
  spawnSync(process.execPath, [join(ROOT, 'tool/gen-u2-kinds.mjs'), ...args], {cwd, encoding: 'utf8'});

/** The fixture's kinds as the generated file holds them. */
async function fixtureSpecs(): Promise<U2KindSpec[]> {
  const out = join(tmp, 'fixture.g.ts');
  assert.equal(generate(['--manifest', FIXTURE, '--out', out]).status, 0);
  return (await import(pathToFileURL(out).href)).U2_KINDS;
}

const defsOf = (specs: U2KindSpec[], overrides = {}): Map<string, KindDef> =>
  new Map(kindDefsFrom(specs, overrides).map(({name, def}) => [name, def]));

/** What `after` lost of `before`, one line per loss: a key that no longer resolves; a selector, label
 * or editor alternative, a strategy, a part or one of its alternatives, a Dart name that is gone;
 * today's strategies tried in another order or after a new one; other gestures. Alternatives are
 * compared as top-level comma parts (`REPLACED` lists the dropped ones); strategies in the R9 order,
 * where a new one may only follow today's (a label or title strategy without a label selector and a
 * dart one without names find nothing) — `item` may also come before `name` and `dart`: it decides only
 * the phrases those decide today, and the resolve trace (BDD_RESOLVE_TRACE) reports each of them. */
function compareKinds(before: Record<string, KindDef & {name?: string}>, after: (key: string) => KindDef | undefined): string[] {
  const lost: string[] = [];
  const list = (s?: string): string[] => s === undefined ? [] : selectorList(s);
  const keep = (what: string, b: string[] | undefined, a: string[] | undefined, replaced: Record<string, string | null> = {}) => {
    for (const x of b ?? []) {
      const by = replaced[x];
      if (!(a ?? []).includes(x) && !(by === null || by !== undefined && (a ?? []).includes(by)))
        lost.push(`${what} ${x}`);
    }
  };
  for (const [key, b] of Object.entries(before)) {
    const a = after(key);
    if (a === undefined) {
      lost.push(`"${key}" no longer resolves`);
      continue;
    }
    keep(`"${key}" selector`, list(b.selector), list(a.selector), REPLACED[b.name ?? key]);
    keep(`"${key}" match`, b.match, a.match);
    const today = effectiveOrder(b.match);
    const inert = (s: MatchStrategy) => (s === 'label' || s === 'title') && !a.labelSelector || s === 'dart' && !a.dartNames?.length;
    const lateItem = (s: MatchStrategy, i: number, all: MatchStrategy[]) => s === 'item' && !today.includes('item') &&
      all.slice(i + 1).every((x) => x === 'name' || x === 'dart');
    const order = effectiveOrder(a.match).filter((s, i, all) => today.includes(s) || !inert(s) && !lateItem(s, i, all));
    if (!isDeepStrictEqual(order.slice(0, today.length), today))
      lost.push(`"${key}" match order ${today.join(',')} → ${order.join(',')}`);
    keep(`"${key}" label`, list(b.labelSelector), list(a.labelSelector));
    keep(`"${key}" editor`, list(b.editorSelector), list(a.editorSelector));
    keep(`"${key}" dart name`, b.dartNames, a.dartNames);
    for (const [part, selector] of Object.entries(b.parts ?? {})) {
      if (a.parts?.[part] === undefined)
        lost.push(`"${key}" part "${part}"`);
      else
        keep(`"${key}" part "${part}"`, list(selector), list(a.parts[part]));
    }
    if (!isDeepStrictEqual(b.gestures, a.gestures))
      lost.push(`"${key}" gestures ${JSON.stringify(b.gestures)} → ${JSON.stringify(a.gestures)}`);
  }
  return lost;
}

test('the generator writes the same file from any directory, sorted by kind, and checks it', () => {
  const [a, b] = [join(tmp, 'a.ts'), join(tmp, 'b.ts')];
  assert.equal(generate(['--manifest', FIXTURE, '--out', a]).status, 0);
  assert.equal(generate(['--manifest', FIXTURE, '--out', b], tmp).status, 0);
  const text = readFileSync(a, 'utf8');
  assert.equal(readFileSync(b, 'utf8'), text);
  assert.match(text, /^\/\* generated by tool\/gen-u2-kinds\.mjs from libraries\/bdd\/tests\/fixtures\/manifest-kinds\.json\n/);
  const kinds = [...text.matchAll(/^ {4}"kind": "([^"]+)"/gm)].map((m) => m[1]);
  assert.deepEqual(kinds, [...kinds].sort());
  assert.equal(kinds.length, 14);
  assert.equal(generate(['--manifest', FIXTURE, '--out', a, '--check']).status, 0);
  writeFileSync(a, text.replace('"noun": "dialog"', '"noun": "modal"'));
  const stale = generate(['--manifest', FIXTURE, '--out', a, '--check']);
  assert.equal(stale.status, 1);
  assert.match(stale.stderr, /a\.ts is stale[\s\S]*"noun": "modal"/);
  const bare = join(tmp, 'bare.json');
  writeFileSync(bare, '{"components": []}');
  const none = generate(['--manifest', bare, '--out', join(tmp, 'c.ts')]);
  assert.equal(none.status, 2);
  assert.match(none.stderr, /has no "kinds" yet/);
});

test('the kinds sharing a noun are one kind; a part of another kind is no noun', async () => {
  const specs = await fixtureSpecs();
  assert.deepEqual(specs.map((s) => s.source).filter((s) => s === 'ems'), ['ems']);
  const defs = defsOf(specs);
  assert.deepEqual([...defs.keys()].sort(), ['button', 'checkbox', 'choice input', 'combobox', 'dialog', 'domain list',
    'functions browser', 'list', 'number input', 'text input']);
  const number = defs.get('number input')!;
  assert.equal(number.selector, '[data-u2="bigint-input"], [data-u2="number-input"], [data-u2="qnum-input"]');
  assert.deepEqual(number.aliases, ['numeric input', 'number field']);
  assert.deepEqual(number.match, ['label', 'aria', 'placeholder', 'name']);
  assert.equal(number.labelSelector, '[data-u2-part="label"]');
  assert.equal(number.editorSelector, '[data-u2-part="editor"] input, [data-u2-part="editor"] select, ' +
    '[data-u2-part="editor"] textarea, [data-u2-part="editor"] [contenteditable]');
  assert.equal(number.parts!.error, '[data-u2-part="error"]');
  const dialog = defs.get('dialog')!;
  assert.equal(dialog.labelSelector, '.u2-dialog-title-text');
  assert.equal(dialog.editorSelector, undefined);
  assert.deepEqual(defs.get('button')!.match, ['text', 'aria', 'title', 'name']);
  const [button] = specs.filter((s) => s.kind === 'button');
  const dropdown = defsOf([{...button, kind: 'actions-menu', noun: 'dropdown button', nameFrom: ['aria', 'title']},
    {...button, kind: 'dropdown-button', noun: 'dropdown button'}]).get('dropdown button')!;
  assert.deepEqual(dropdown.match, ['text', 'aria', 'title', 'name'], 'the kind the noun names leads');
});

test('an override adds selectors, names and parts, and replaces the strategies and the label', async () => {
  const specs = await fixtureSpecs();
  const defs = defsOf(specs, {
    'text input': {selector: '.ui-input-text', aliases: ['textbox', 'input box'], match: ['name', 'label', 'placeholder', 'dart'],
      labelSelector: '.ui-input-label', dartNames: ['input-host-{q}'], parts: {editor: '.ui-input-editor'}},
    'dialog': {dartNames: ['dialog-{q}']},
    'list': {u2Selector: '[data-u2="list"]:not([role="tree"])', selector: '[role="listbox"]'},
  });
  assert.equal(defs.get('list')!.selector, '[data-u2="list"]:not([role="tree"]), [role="listbox"]');
  const text = defs.get('text input')!;
  assert.equal(text.selector, '[data-u2="text-input"], .ui-input-text');
  assert.deepEqual(text.aliases, ['text field', 'textbox', 'input box']);
  assert.deepEqual(text.match, ['name', 'label', 'placeholder', 'dart']);
  assert.equal(text.labelSelector, '.ui-input-label');
  assert.deepEqual(text.dartNames, ['input-host-{q}']);
  assert.equal(text.parts!.editor, '.ui-input-editor');
  assert.equal(text.parts!.label, '[data-u2-part="label"]');
  assert.deepEqual(defs.get('dialog')!.match, ['title', 'aria', 'name', 'dart']);
  assert.throws(() => defsOf(specs, {'modal': {}}), /no manifest kind has the noun "modal"/);
});

test('a part of another kind is a part of its owner, by its noun or a legacy name', async () => {
  const specs = await fixtureSpecs();
  assert.deepEqual(defsOf(specs).get('functions browser')!.parts,
    {'list': '[data-u2="fb-list"]', 'clear button': '[data-u2="fb-clear"]', 'search': '[data-u2="fb-search"]'});
  const legacy = defsOf(specs, {'functions browser': {partNames: {'fb-search': 'search box'}}}).get('functions browser')!;
  assert.equal(legacy.parts!['search box'], '[data-u2="fb-search"]');
  assert.equal(legacy.parts!.search, undefined);
});

test('an item noun of a collection is an item kind, a bare generic word is no alias', async () => {
  const [button] = (await fixtureSpecs()).filter((s) => s.kind === 'button');
  const specs: U2KindSpec[] = [
    {...button, kind: 'chips-input', noun: 'chips input', aliases: ['history', 'chip list'], items: {role: 'button', noun: 'chip', gestures: ['toggle']}},
    {...button, kind: 'chip', noun: 'chip'},
  ];
  const defs = defsOf(specs);
  assert.equal(defs.get('chip')!.items, true);
  assert.equal(defs.get('chips input')!.items, undefined);
  assert.deepEqual(defs.get('chips input')!.aliases, ['chip list']);
  assert.equal(defsOf(specs, {chip: {items: false}}).get('chip')!.items, undefined);
});

test('registering the u2 kinds makes their nouns resolve', async () => {
  assert.equal(lookupKind('domain catalog'), undefined);
  registerU2Kinds((await fixtureSpecs()).filter((s) => s.source === 'ems').map((s) => ({...s, noun: 'domain catalog'})), {});
  assert.equal(lookupKind('domain catalog')!.selector, '[data-u2="domain-list"]');
});

test('the vocabulary is the one the manifests generate', (t) => {
  if (!existsSync(join(ROOT, '../u2/manifest.json'))) {
    t.skip('no libraries/u2/manifest.json in this checkout');
    return;
  }
  const check = generate(['--check']);
  assert.equal(check.status, 0, check.stderr);
});

test('every manifest kind has a noun, and no selector names a kind the manifests lack', () => {
  const kinds = allKinds();
  const named = (selectors: (string | undefined)[]) => new Set(selectors.flatMap((s) => [...(s ?? '').matchAll(/\[data-u2="([^"]+)"\]/g)])
    .map((m) => m[1]));
  const roots = named(kinds.map((k) => k.selector));
  assert.deepEqual(U2_KINDS.filter((s) => s.partOf === undefined && !LAYOUT_ONLY.includes(s.kind) && !roots.has(s.kind)).map((s) => s.kind), []);
  const manifest = new Set(U2_KINDS.map((s) => s.kind));
  const all = named(kinds.flatMap((k) => [k.selector, k.labelSelector, k.editorSelector, ...Object.values(k.parts ?? {})]));
  assert.deepEqual([...all].filter((k) => !manifest.has(k)), []);
});

test('a u2 noun and its aliases resolve to the generated kind, never to a hand-written one', () => {
  for (const spec of U2_KINDS.filter((s) => s.partOf === undefined)) {
    assert.equal(lookupKind(spec.noun)?.name, spec.noun);
    for (const alias of spec.aliases ?? [])
      assert.equal(lookupKind(alias)?.name ?? spec.noun, spec.noun, `"${alias}" of ${spec.kind}`);
  }
});

test('items of a collection are item kinds, matched by key after what the user sees; a control is no item', () => {
  for (const name of ITEM_KINDS) {
    const k = lookupKind(name)!;
    assert.equal(k.items, true, name);
    const order = effectiveOrder(k.match);
    assert.ok(order.slice(order.indexOf('item') + 1).every((s) => s === 'name' || s === 'dart'), `${name}: ${order.join(',')}`);
  }
  for (const name of ['button', 'checkbox', 'tag', 'input'])
    assert.equal(lookupKind(name)!.items, undefined, name);
});

test('the kind table keeps everything the snapshot has', () => {
  assert.ok(Object.keys(BEFORE).length > 150);
  assert.deepEqual(compareKinds(BEFORE, lookupKind), []);
});

test('the compat check reports what a rewrite lost', () => {
  const changed = (key: string): KindDef | undefined => {
    const k = lookupKind(key);
    if (key === 'dropdown' || k === undefined)
      return undefined;
    switch (k.name) {
      case 'text input': return {...k, selector: k.selector.replace(', .ui-input-text', ''), match: ['label', 'aria', 'placeholder', 'name', 'dart']};
      case 'button': return {...k, match: ['text', 'name'], dartNames: []};
      case 'icon button': return {...k, match: ['aria', 'title', 'name', 'dart']};
      case 'notification': return {...k, match: ['text', 'name']};
      case 'column input': return {...k, match: ['dart', 'label', 'name'], gestures: undefined};
      case 'dialog': return {...k, parts: {...k.parts, footer: '.u2-dialog-footer'}};
      case 'tab': return {...k, match: ['name', 'item', 'label', 'text', 'aria', 'dart']};
      case 'item': return {...k, selector: k.selector.replace('.u2-list-row, ', '')};
      default: return k;
    }
  };
  const lost = compareKinds(BEFORE, changed);
  assert.ok(lost.some((l) => l.startsWith('"tab" match order ')), 'an item key tried before the label');
  assert.ok(lost.includes('"item" selector .u2-list-row'));
  assert.ok(!lost.some((l) => l.includes('[data-u2="item"]')), 'a retired alternative');
  assert.ok(lost.includes('"dropdown" no longer resolves'));
  assert.ok(lost.includes('"text input" selector .ui-input-text'));
  assert.ok(lost.includes('"textbox" selector .ui-input-text'));
  assert.ok(lost.includes('"button" match aria'));
  assert.ok(lost.includes('"button" dart name button-{q}'));
  assert.ok(lost.includes('"column input" match order label,accessible,name,dart → dart,label,accessible,name'));
  assert.ok(lost.includes('"column input" gestures {"open":"mousedown"} → undefined'));
  assert.ok(lost.includes('"dialog" part "footer" .d4-dialog-footer'));
  assert.ok(lost.some((l) => l.startsWith('"text input" match order ')), 'aria moved before the placeholder');
  assert.ok(!lost.some((l) => l.startsWith('"icon button"') || l.startsWith('"notification"')),
    'a title strategy without a label selector, and a strategy after today\'s, change nothing');
});

let missing = '';
let browser: Browser | undefined;
let page: Page | undefined;

// no __name shim on this page: what runs in it must not depend on one
before(async () => {
  try {
    browser = await chromium.launch();
  }
  catch (e) {
    missing = `Chromium does not start here: ${String(e).slice(0, 80)}`;
    return;
  }
  page = await browser.newPage();
  await page.goto('about:blank');
});

after(async () => {
  await browser?.close();
});

const BRIDGED = '<div data-u2="dart-input"><div class="ui-input-root ui-input-text" name="input-host-Weight">' +
  '<label class="ui-input-label">Weight</label><input class="ui-input-editor" value="70"></div></div>' +
  '<div class="ui-input-root ui-input-text" name="input-host-Height"><label class="ui-input-label">Height</label>' +
  '<input class="ui-input-editor" value="180"></div>' +
  '<div data-u2="chips-input"><button class="u2-chips-chip" data-u2-item="k1"><span data-u2-part="item-label">Kinase</span></button></div>';

test('a u2 input bridging a Dart one is one input; a chip of a chips input goes by its key', async (t) => {
  if (page === undefined) {
    t.skip(missing || 'no page');
    return;
  }
  await page.setContent(BRIDGED);
  const bothRoots = page.locator('[data-u2$="-input"], .ui-input-root').filter({has: page.locator('.ui-input-label', {hasText: 'Weight'})});
  assert.equal(await bothRoots.count(), 2, 'the field the selector used to find twice');
  const weight = await locate(page, el('Weight input'));
  assert.equal(await weight.count(), 1);
  assert.equal(await weight.getAttribute('data-u2'), 'dart-input');
  assert.equal(await (await locate(page, el('Weight text input'))).getAttribute('name'), 'input-host-Weight');
  assert.equal(await (await locate(page, el('Height input'))).count(), 1);
  assert.equal(await (await locate(page, el('"k1" chip'))).textContent(), 'Kinase');
  assert.equal(await (await locate(page, el('Kinase chip'))).count(), 1);
  assert.equal(await describeCandidates(bothRoots), '"Weight" ×2');
});
