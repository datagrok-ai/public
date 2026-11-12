/* R9 addressing against static pages: what the user sees resolves first, data-u2-name is the
   fallback, Dart names keep their place, an ambiguity names its candidates, the trace mode reports a
   phrase whose target moved, and u2 items resolve by key and are revealed for a gesture. Runs in the
   library's Chromium; no stand. */
import assert from 'node:assert/strict';
import {mkdtempSync, readFileSync, rmSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {after, before, test} from 'node:test';
import {type Browser, chromium, type Page} from '@playwright/test';

import '../bindings/common/kinds.js';
import '../bindings/platform/elements.js';
import {kind, lookupKind, type MatchStrategy} from '../src/registry.js';
import {el} from '../src/runtime/args.js';
import {expectState} from '../src/runtime/assertions.js';
import {failure, isWaitFailure} from '../src/runtime/failure.js';
import {click} from '../src/runtime/gestures.js';
import {describeCandidates, effectiveOrder, explain, locate} from '../src/runtime/locate.js';

let missing = '';
let browser: Browser | undefined;
let page: Page | undefined;

const scenario = (name: string, fn: () => Promise<void>): void => void test(name, async (t) => {
  if (page === undefined) {
    t.skip(missing || 'no page');
    return;
  }
  await fn();
});

before(async () => {
  try {
    browser = await chromium.launch();
  }
  catch (e) {
    missing = `Chromium does not start here: ${String(e).slice(0, 80)}`;
    return;
  }
  page = await browser.newPage();
  // tsx names the functions it compiles with a helper the page does not have
  await page.addInitScript(() => { (window as any).__name = (fn: unknown) => fn; });
  await page.goto('about:blank');
});

after(async () => {
  await browser?.close();
});

const INPUT = (label: string, attrs = '', value = '') => `<div data-u2="text-input"${attrs}>` +
  `<label data-u2-part="label">${label}</label><span data-u2-part="editor"><input value="${value}"></span></div>`;
const LABEL_OR_NAME = INPUT('Organization', '', 'A') + INPUT('Company', ' data-u2-name="organization"', 'B') +
  '<div data-u2="number-input" data-u2-name="dose"><span data-u2-part="editor"><input value="5"></span></div>';
const valueOf = async (phrase: string): Promise<string> => (await locate(page!, el(phrase))).locator('input').inputValue();
const count = async (phrase: string): Promise<number> => (await locate(page!, el(phrase))).count();

test('the R9 order: what the user sees first, data-u2-name after it, dart where it was', () => {
  const order = (match: MatchStrategy[]) => effectiveOrder(match).join(',');
  assert.equal(order(['name', 'label', 'placeholder', 'dart']), 'label,accessible,placeholder,name,dart');
  assert.equal(order(['name', 'text', 'aria', 'dart']), 'text,aria,accessible,name,dart');
  assert.equal(order(['dart', 'name']), 'dart,name');
  assert.equal(order(['name', 'dart']), 'name,dart');
  assert.equal(order(['name', 'title', 'dart']), 'title,name,dart');
});

scenario('a label the user reads wins over a data-u2-name that spells the same word', async () => {
  await page!.setContent(LABEL_OR_NAME);
  assert.equal(await count('Organization input'), 1);
  assert.equal(await valueOf('Organization input'), 'A');
});

scenario('an input with no visible name still resolves by its data-u2-name', async () => {
  await page!.setContent(LABEL_OR_NAME);
  assert.equal(await valueOf('dose input'), '5');
});

scenario('a label the browser associates from outside the root (aria-labelledby, label for) resolves the input', async () => {
  await page!.setContent(`<span id="weight-label">Weight, kg</span>
    <div data-u2="number-input"><span data-u2-part="editor"><input aria-labelledby="weight-label" value="70"></span></div>
    <label for="height-editor">Height</label>
    <div data-u2="number-input"><span data-u2-part="editor"><input id="height-editor" value="180"></span></div>`);
  assert.equal(await valueOf('"Weight, kg" input'), '70');
  assert.equal(await valueOf('Height number input'), '180');
});

scenario('a placeholder still names an input that shows no label', async () => {
  await page!.setContent('<div data-u2="typeahead"><span data-u2-part="editor"><input placeholder="Search users…"></span></div>');
  assert.equal(await count('"Search users…" typeahead'), 1);
});

scenario('a Dart input resolves by its label, and by its input-host name, as before', async () => {
  await page!.setContent(`<div class="ui-input-root ui-input-text" name="input-host-Caption"><label class="ui-input-label">Caption</label>
    <input class="ui-input-editor" value="x"></div>
    <div class="ui-input-root ui-input-int" name="input-host-Count"><input class="ui-input-editor" value="3"></div>`);
  assert.equal(await valueOf('Caption input'), 'x');
  assert.equal(await valueOf('Count input'), '3');
  assert.equal(await valueOf('Count number input'), '3');
});

scenario('two visible matches fail the gesture, and the report names each by its label and its section', async () => {
  await page!.setContent(`<div data-u2="dialog"><div class="u2-dialog-title"><span class="u2-dialog-title-text">Visit</span></div>
    <div data-u2="section"><div class="u2-section-header"><span class="u2-section-title">Patient</span></div>${INPUT('Name')}</div>
    <div data-u2="section" data-u2-name="doctor"><div data-u2-part="title">Doctor</div>${INPUT('Name', ' data-u2-name="doctorName"')}</div>
  </div>`);
  const e = await click(page!, el('Name text input')).then(() => undefined, (x: unknown) => x);
  assert.ok(e, 'the click must not pick one of two inputs');
  assert.ok(isWaitFailure(e));
  const f = failure('features/visit.feature:9', 'When user clicks on Name text input', e, await explain(page!));
  const report = f.message.split('\n');
  assert.match(report[3], /^strict mode violation: .* resolved to 2 elements/);
  assert.match(f.message, /visible text inputs: "Name" in section "Patient" \| "Name" in section "Doctor" \[name=doctorName\]/);
});

scenario('candidates go by the label the u2 status reports for a control root, and identical ones are counted', async () => {
  await page!.setContent(`<div data-u2="list" data-fake="Results" data-u2-name="recent"><div class="u2-list-row">first row</div></div>
    <div data-u2="list" data-fake="Results">x</div><div data-u2="list" data-fake="Results">y</div>`);
  await page!.evaluate(() => {
    (window as any).DG = {U2: {
      control: (e: Element) => e.closest('[data-fake]') ? {root: e.closest('[data-fake]')} : undefined,
      state: (e: Element) => ({label: e.getAttribute('data-fake')}),
    }};
  });
  assert.equal(await describeCandidates(page!.locator('[data-u2="list"]')), '"Results" [name=recent] | "Results" ×2');
  assert.equal(await describeCandidates(page!.locator('.u2-list-row')), '"first row"', 'a row is not named by the list it is in');
  await page!.evaluate(() => { delete (window as any).DG; });
});

scenario('the trace records a phrase the R9 order sends to another element, and nothing for one it does not move', async () => {
  const dir = mkdtempSync(join(tmpdir(), 'bdd-trace-'));
  const file = join(dir, 'trace.jsonl');
  process.env.BDD_RESOLVE_TRACE = file;
  try {
    await page!.setContent(LABEL_OR_NAME);
    await locate(page!, el('Organization input'));
    await locate(page!, el('dose input'));
    const lines = readFileSync(file, 'utf8').trim().split('\n').map((l) => JSON.parse(l));
    assert.equal(lines.length, 1);
    assert.deepEqual({...lines[0], test: undefined}, {test: undefined, phrase: 'Organization input', kind: 'input',
      qualifier: 'organization', legacy: {strategy: 'name', count: 1}, chosen: {strategy: 'label', count: 1}, same: false});
  }
  finally {
    delete process.env.BDD_RESOLVE_TRACE;
    rmSync(dir, {recursive: true, force: true});
  }
});

scenario('an item a gesture names in its collection is revealed by u2 when not rendered; an absence check reveals nothing', async () => {
  await page!.setContent('<div data-u2="list" data-u2-name="cities"><div class="u2-list-row" role="option" data-u2-item="oslo">Oslo</div></div>');
  await page!.evaluate(() => {
    const w = window as any;
    w.__revealed = [];
    w.DG = {U2: {revealItem: async (list: Element, q: string) => {
      w.__revealed.push(q);
      const row = document.createElement('div');
      row.className = 'u2-list-row';
      row.setAttribute('role', 'option');
      row.setAttribute('data-u2-item', q.toLowerCase());
      row.textContent = q;
      row.addEventListener('click', () => row.setAttribute('data-clicked', 'true'));
      list.appendChild(row);
      return row;
    }}};
  });
  await click(page!, el('Bergen item in cities list'));
  assert.equal(await page!.locator('[data-u2-item="bergen"]').getAttribute('data-clicked'), 'true');
  await click(page!, el('Oslo item in cities list'));
  await expectState(page!, el('Tromsø item in cities list'), 'absent');
  assert.deepEqual(await page!.evaluate(() => (window as any).__revealed), ['Bergen'], 'asked once, as the feature spells it');
  await page!.evaluate(() => {
    (window as any).DG.U2.revealItem = async () => { throw new Error('u2: 2 items match "narvik": narvik-1, narvik-2'); };
  });
  await assert.rejects(click(page!, el('Narvik item in cities list')), /u2: 2 items match "narvik": narvik-1, narvik-2/);
  await page!.evaluate(() => { delete (window as any).DG; });
});

scenario('an item of a u2 collection resolves by its key where its text says something else', async () => {
  kind('item', {...lookupKind('item')!, match: ['text', 'item', 'label', 'aria', 'name'], items: true});
  await page!.setContent(`<div data-u2="list" data-u2-name="cities">
    <div class="u2-list-row" role="option" data-u2-item="kyiv">Київ</div><div class="u2-list-row" role="option" data-u2-item="oslo">Oslo</div></div>`);
  assert.equal(await count('"kyiv" item in cities list'), 1);
  assert.equal(await (await locate(page!, el('"kyiv" item in cities list'))).textContent(), 'Київ');
});
