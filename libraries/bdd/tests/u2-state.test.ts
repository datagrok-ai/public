/* State and value claims through the u2 status against a static page with a stand-in `DG.U2`: the
   status decides for a control root that reports the key, the DOM for everything else, and the
   u2-only phrases fail in one sentence. Runs in the library's Chromium; no stand. */
import assert from 'node:assert/strict';
import {after, before, test} from 'node:test';
import {type Browser, chromium, type Page} from '@playwright/test';

import '../bindings/common/kinds.js';
import {itemShouldBe, shouldBeOnStep, shouldHaveError, shouldHaveItemsInTotal, shouldHaveNoError, shouldHaveNothingSelected,
  shouldHaveSelection, shouldListItem, shouldNotListItem, shouldReportValue, shouldShowTab} from '../bindings/common/u2-steps.js';
import {el} from '../src/runtime/args.js';
import {expectCount, expectState, expectText, expectValue, readValue} from '../src/runtime/assertions.js';
import {reasonOf} from '../src/runtime/failure.js';
import {locate} from '../src/runtime/locate.js';
import {whileExpectedToFail} from '../src/runtime/patience.js';
import {readItems, readState, summary} from '../src/runtime/u2-state.js';

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

/** The page, with a `DG.U2` that answers for every element under a `data-fake-state` root. */
async function show(html: string): Promise<void> {
  await page!.setContent(html);
  await page!.evaluate(() => {
    const root = (e: Element) => e.closest('[data-fake-state]');
    const state = (e: Element) => JSON.parse(root(e)!.getAttribute('data-fake-state')!);
    (window as any).DG = {U2: {
      control: (e: Element) => root(e) ? {root: root(e)} : undefined,
      state: (e: Element) => root(e) ? state(e) : null,
      items: (e: Element) => root(e) ? state(e).items ?? null : null,
    }};
  });
}

function fake(state: object): string {
  return JSON.stringify({enabled: true, visible: true, focused: false, ...state}).replace(/'/g, '&#39;');
}

function input(label: string, state: object, editor = '<input>', attrs = ''): string {
  return `<div data-u2="text-input" data-fake-state='${fake({kind: 'text-input', label, ...state})}' ${attrs}>` +
    `<label data-u2-part="label">${label}</label><span data-u2-part="editor">${editor}</span><span data-u2-part="error"></span></div>`;
}

/** The reasons failing claims report, failing fast and together: the claims poll, and a stand-in page
 * never changes. */
async function reasonsFor(...bodies: (() => Promise<unknown>)[]): Promise<string[]> {
  return whileExpectedToFail(() => Promise.all(bodies.map(async (body) => {
    try {
      await body();
    }
    catch (e) {
      return reasonOf(e);
    }
    return assert.fail('the claim held');
  })));
}

/** What the claims printed to the console while `body` ran. */
async function warnings(body: () => Promise<void>): Promise<string[]> {
  const out: string[] = [];
  const warn = console.warn;
  console.warn = (...args: unknown[]) => { out.push(args.join(' ')); };
  try {
    await body();
  }
  finally {
    console.warn = warn;
  }
  return out;
}

scenario('the status decides a state the DOM says otherwise, and the disagreement prints once per phrase', async () => {
  await show(input('Code', {enabled: false}, '<input value="A1">'));
  const printed = await warnings(async () => {
    await expectState(page!, el('Code input'), 'disabled');
    await expectState(page!, el('Code input'), 'disabled');
    const [reason] = await reasonsFor(() => expectState(page!, el('Code input'), 'enabled'));
    assert.match(reason, /^enabled expected/);
  });
  assert.deepEqual(printed, ['bdd: u2 contract: Code input: the status says disabled, the DOM says enabled (aria-disabled)']);
});

scenario('a key the status does not report, a Dart input and a part are read from the DOM', async () => {
  await show(input('Dose', {}, '<input aria-invalid="true">') +
    '<div class="ui-input-root ui-input-text" name="input-host-Caption"><label class="ui-input-label">Caption</label>' +
    '<input class="ui-input-editor d4-invalid" value="x" readonly required></div>' +
    input('Lot', {enabled: false}));
  const printed = await warnings(async () => {
    await expectState(page!, el('Dose input'), 'invalid');
    await expectState(page!, el('Caption input'), 'invalid');
    await expectState(page!, el('Caption input'), 'read-only');
    await expectState(page!, el('Caption input'), 'required');
    // the part is not the control's root: the status of the control around it does not answer for it
    await expectState(page!, el('error of Lot input'), 'enabled');
  });
  assert.deepEqual(printed, []);
  assert.equal((await readState(page!.locator('[name="input-host-Caption"]')))?.state, null);
});

scenario('the value is the text the status reports, and a failure names the control', async () => {
  await show(`<div data-u2="choice-input" data-fake-state='${fake({kind: 'choice-input', label: 'Stage', value: 'disc', text: 'Discovery', valid: true})}'>
    <label data-u2-part="label">Stage</label><span data-u2-part="editor"><select><option>Discovery</option><option>Phase II</option></select></span></div>
    ${input('Batch', {value: [1, 2]}, '<input value="1, 2">')}${input('Lab', {text: 'North'}, '<input value="South">')}`);
  await expectValue(page!, el('Stage input'), 'Discovery');
  await expectValue(page!, el('Batch input'), '1, 2');
  assert.deepEqual(await reasonsFor(() => expectValue(page!, el('Stage input'), 'Phase II'), () => expectValue(page!, el('Stage input'), 'Discovery', true)), [
    'the value of Stage input is "Discovery", expected "Phase II" (u2: choice-input "Stage", enabled, valid)',
    'the value of Stage input is "Discovery", expected anything but "Discovery" (u2: choice-input "Stage", enabled, valid)']);
  const printed = await warnings(async () => assert.equal(await readValue(page!, el('Lab input')), 'North'));
  assert.deepEqual(printed, ['bdd: u2 contract: Lab input: the status says "North", the DOM says "South" (the editor\'s value)']);
  assert.equal(summary((await readState(await locate(page!, el('Stage input'))))!.state!),
    'choice-input "Stage", enabled, valid, value "Discovery"');
});

scenario('read-only, editable, busy and required, from the status and from the DOM', async () => {
  await show(`${input('Owner', {readOnly: true, required: false, busy: false})}
    ${input('Plain', {}, '<input readonly aria-required="true">')}
    <div aria-busy="true">${input('Loading', {})}</div>`);
  const printed = await warnings(async () => {
    await expectState(page!, el('Owner input'), 'read-only');
    const [reason] = await reasonsFor(() => expectState(page!, el('Owner input'), 'editable'));
    assert.match(reason, /^editable expected/);
  });
  assert.deepEqual(printed, ['bdd: u2 contract: Owner input: the status says read-only, the DOM says editable (aria-readonly)']);
  await expectState(page!, el('Owner input'), 'required', true);
  await expectState(page!, el('Owner input'), 'busy', true);
  await expectState(page!, el('Plain input'), 'read-only');
  await expectState(page!, el('Plain input'), 'required');
  await expectState(page!, el('Plain input'), 'busy', true);
  await expectState(page!, el('Loading input'), 'busy');
  await expectState(page!, el('Owner input'), 'editable', true);
});

const STAGES = fake({kind: 'list', label: 'Stages', selection: ['d', 'p2'], items: {count: 3, truncated: false, list: [
  {key: 'd', label: 'Discovery', index: 0, selected: true, rendered: true},
  {key: 'p1', label: 'Phase I', index: 1, selected: false, disabled: true, rendered: true},
  {key: 'p2', label: 'Phase II', index: 2, selected: true, rendered: false},
]}});

scenario('all the items, one item\'s state and the selection come from the status', async () => {
  await show(`<div data-u2="list" data-u2-name="stages" data-fake-state='${STAGES}'>
    <div class="u2-list-row" role="option" data-u2-item="d">Discovery</div><div class="u2-list-row" role="option" data-u2-item="p1">Phase I</div></div>`);
  const stages = el('stages list');
  await shouldHaveItemsInTotal(page!, stages, 3);
  await expectCount(page!, stages, 2);
  await shouldListItem(page!, stages, 'phase ii');
  await shouldListItem(page!, stages, 'p2');
  await shouldNotListItem(page!, stages, 'Phase III');
  await itemShouldBe(page!, 'Discovery', stages, 'selected');
  await itemShouldBe(page!, 'Phase I', stages, 'disabled');
  await itemShouldBe(page!, 'Discovery', stages, 'enabled');
  await itemShouldBe(page!, 'p2', stages, 'valid');
  await shouldHaveSelection(page!, stages, 'Discovery, p2');
  assert.deepEqual(await reasonsFor(
    () => shouldListItem(page!, stages, 'Phase III'),
    () => shouldNotListItem(page!, stages, 'Discovery'),
    () => itemShouldBe(page!, 'Phase I', stages, 'selected'),
    () => itemShouldBe(page!, 'Discovery', stages, 'checked'),
    () => shouldHaveSelection(page!, stages, 'Discovery'),
    () => shouldHaveNothingSelected(page!, stages),
    () => shouldHaveItemsInTotal(page!, stages, 2),
  ), [
    'stages list does not list "Phase III"; it lists: "Discovery", "Phase I", "Phase II"',
    'stages list lists "Discovery" (key d)',
    'item "Phase I" of stages list is not selected; it reports: selected false, disabled true, rendered true',
    'item "Discovery" of stages list reports no checked state; it reports: selected true, rendered true',
    'the selection of stages list is "Discovery, Phase II", expected "Discovery"',
    'stages list has the selection "Discovery, Phase II"',
    'stages list has 3 items in total, expected 2',
  ]);
  assert.equal((await readItems(page!.locator('[data-u2="list"]'), true))?.count, 3);
});

scenario('the current tab and step come from the status, by key or by label', async () => {
  const strip = (kind: string, name: string, current: string) => `<div data-u2="${kind}" data-u2-name="${name}" data-fake-state='${fake({kind, current,
    items: {count: 2, truncated: false, list: [{key: 'data', label: 'Data', index: 0, rendered: true}, {key: 'ana', label: 'Analysis', index: 1, rendered: true}]}})}'>x</div>`;
  await show(strip('tabs', 'results', 'ana') + strip('wizard', 'setup', 'data'));
  await shouldShowTab(page!, el('results tabs'), 'Analysis');
  await shouldShowTab(page!, el('results tabs'), 'ana');
  await shouldBeOnStep(page!, el('setup wizard'), 'data');
  assert.deepEqual(await reasonsFor(
    () => shouldShowTab(page!, el('results tabs'), 'Data'),
    () => shouldBeOnStep(page!, el('setup wizard'), 'Analysis'),
    () => shouldHaveNothingSelected(page!, el('results tabs')),
  ), [
    'results tabs shows the "Analysis" tab, expected "Data"',
    'setup wizard is on the "Data" step, expected "Analysis"',
    'results tabs (u2: tabs, enabled) reports no selection',
  ]);
});

scenario('the error and the raw value: the status, the error part as the fallback, and nothing for a Dart input', async () => {
  await show(`${input('Code', {valid: false, error: 'At least 3 characters', value: 12, text: '12.00'})}
    ${input('Ref', {}, '<input>').replace('<span data-u2-part="error"></span>', '<span data-u2-part="error">Required</span>')}
    ${input('Memo', {valid: true, error: null, value: {a: [1, 'x']}})}
    <div class="ui-input-root ui-input-text" name="input-host-Caption"><label class="ui-input-label">Caption</label><input class="ui-input-editor"></div>`);
  await shouldHaveError(page!, el('Code input'), 'At least 3 characters');
  await shouldHaveError(page!, el('Ref input'), 'Required');
  await shouldHaveNoError(page!, el('Memo input'));
  await shouldReportValue(page!, el('Code input'), '12');
  await shouldReportValue(page!, el('Memo input'), '{"a": [1, "x"]}');
  assert.deepEqual(await reasonsFor(
    () => shouldHaveNoError(page!, el('Code input')),
    () => shouldHaveError(page!, el('Memo input'), 'x'),
    () => shouldReportValue(page!, el('Code input'), '13'),
    () => shouldReportValue(page!, el('Caption input'), 'x'),
    () => shouldHaveNoError(page!, el('Caption input')),
  ), [
    'Code input has the error "At least 3 characters"',
    'the error of Memo input is none, expected "x"',
    'Code input reports the value 12, expected 13',
    'Caption input is not a u2 control; it reports no value',
    'Caption input is not a u2 control; it reports no error',
  ]);
});

scenario('a hidden first match answers no value, state or text read while a visible one exists; with none visible it still does', async () => {
  await show(`<div hidden>${input('Query', {value: 'old', valid: false})}
    <div class="ui-input-root" name="input-host-Seed"><input class="ui-input-editor" value="1"></div><div data-u2-name="rows">2 of 6</div></div>
    ${input('Query', {value: 'city = "Basel"', valid: true})}
    <div class="ui-input-root" name="input-host-Seed"><input class="ui-input-editor" value="2"></div><div data-u2-name="rows">6 of 6</div>
    <div data-u2-name="away" hidden>Kept</div><div data-u2-name="away" hidden>Other</div>`);
  await expectValue(page!, el('Query input'), 'city = "Basel"');
  assert.equal((await readState(await locate(page!, el('Query input'))))?.state?.valid, true);
  assert.equal(await readValue(page!, el('Seed input')), '2');
  await expectText(page!, el('rows element'), '6 of 6', {exact: true});
  await expectText(page!, el('rows element'), '2 of 6', {negate: true});
  await expectText(page!, el('away element'), 'Kept');
  const [value, text] = await reasonsFor(() => expectValue(page!, el('Query input'), 'old'), () => expectText(page!, el('rows element'), '2 of 6'));
  assert.equal(value, 'the value of Query input is "city = "Basel"", expected "old" (u2: text-input "Query", enabled, valid)');
  assert.match(text, /^the matches of rows element with the text "2 of 6"/);
});

scenario('two visible matches: a u2 phrase names them, a state claim keeps its quantifier', async () => {
  await show(`<div data-u2="section" data-u2-name="patient"><div class="u2-section-title">Patient</div>${input('Name', {enabled: false, value: 'a'})}</div>
    <div data-u2="section" data-u2-name="doctor"><div class="u2-section-title">Doctor</div>${input('Name', {enabled: false, value: 'b'})}</div>`);
  assert.deepEqual(await reasonsFor(() => shouldReportValue(page!, el('Name input'), 'a')),
    ['"Name input" matches 2 visible inputs: "Name" in section "Patient" | "Name" in section "Doctor"']);
  await shouldReportValue(page!, el('first Name input'), 'a');
  const printed = await warnings(async () => {
    await expectState(page!, el('Name input'), 'disabled');
    await page!.evaluate(() => {
      const second = document.querySelectorAll('[data-fake-state]')[1];
      second.setAttribute('data-fake-state', second.getAttribute('data-fake-state')!.replace('"enabled":false', '"enabled":true'));
    });
    const [reason] = await reasonsFor(() => expectState(page!, el('Name input'), 'disabled'));
    assert.match(reason, /^disabled expected/);
  });
  assert.deepEqual(printed, ['bdd: u2 contract: Name input: the status says disabled, the DOM says enabled (aria-disabled)']);
});
