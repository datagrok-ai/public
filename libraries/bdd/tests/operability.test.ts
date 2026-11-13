/* The operability audit, the balloon floor and their two claims over a page running the engine's own probe
   (js-api/dist/src/u2core, served to the library's Chromium): each guard flags the real bug it exists for — the
   old read-only rail rule, a disabled button with no reason, an "Error adding filter" balloon no step claimed —
   and passes its fixed twin; the harness probes a dialog after the step that opened it. No stand. */
import assert from 'node:assert/strict';
import {existsSync, mkdirSync, mkdtempSync, rmSync, writeFileSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {after, before, test} from 'node:test';
import {type Browser, chromium, type Page} from '@playwright/test';

import '../bindings/common/kinds.js';
import {everyControlOperable, shouldBeDisabledBecause} from '../bindings/common/u2-steps.js';
import {errorOrWarningBalloonMatching} from '../bindings/tiers/viewers/steps.js';
import {el} from '../src/runtime/args.js';
import {auditScenario, balloonFindings} from '../src/runtime/audit.js';
import {reasonOf} from '../src/runtime/failure.js';
import {feature} from '../src/runtime/harness.js';
import {readOperability} from '../src/runtime/operability.js';
import {whileExpectedToFail} from '../src/runtime/patience.js';

for (const key of ['BDD_AUDIT', 'BDD_AUDIT_LAYOUT', 'BDD_AUDIT_ERRORS', 'BDD_AUDIT_BALLOONS', 'BDD_AUDIT_OPERABILITY', 'BDD_AUDIT_LEAKS'])
  delete process.env[key];

const tmp = mkdtempSync(join(tmpdir(), 'bdd-operability-'));
function project(name: string, audits: object): string {
  const dir = join(tmp, name);
  mkdirSync(dir, {recursive: true});
  writeFileSync(join(dir, 'bdd.config.json'), JSON.stringify({audits}));
  return dir;
}
const REPORT = project('report', {layout: 'off', errors: 'report', balloons: 'report', operability: 'report'});
const FAIL = project('fail', {layout: 'off', operability: 'fail'});
const PLAIN = project('plain', {});
process.env.BDD_ROOT = REPORT;

const ENGINE = join(import.meta.dirname, '..', '..', '..', 'js-api', 'dist', 'src', 'u2core');
const ORIGIN = 'http://bdd.test';

let missing = '';
let browser: Browser | undefined;
let page: Page | undefined;

const scenario = (name: string, fn: () => Promise<void>): void => void test(name, async (t) => {
  if (page === undefined) {
    t.skip(missing || 'no page');
    return;
  }
  process.env.BDD_ROOT = REPORT;
  await fn();
});

/** The engine's u2core modules served from the js-api build, and loaded as `DG.U2`. */
async function withEngine(p: Page): Promise<void> {
  await p.route(`${ORIGIN}/**`, (route) => {
    const path = new URL(route.request().url()).pathname;
    return path.startsWith('/u2core/') ? route.fulfill({path: join(ENGINE, path.slice('/u2core/'.length)), contentType: 'text/javascript'}) :
      route.fulfill({contentType: 'text/html', body: '<!doctype html><html><head></head><body></body></html>'});
  });
  await p.goto(`${ORIGIN}/`);
  await p.evaluate(async (url) => { (window as any).DG = {U2: await import(url)}; }, `${ORIGIN}/u2core/index.js`);
}

before(async () => {
  if (!existsSync(join(ENGINE, 'inspect.js'))) {
    missing = `js-api is not built: no ${ENGINE}`;
    return;
  }
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
  await withEngine(page);
});

after(async () => {
  await browser?.close();
  rmSync(tmp, {recursive: true, force: true});
});

/** A stand-in shell whose current view holds `view`, `overlays` beside it in the body, every `[data-u2]` element an
 * engine control, and a balloon record (`__bdd.balloons`) the way the harness's in-page runtime keeps it. */
async function standIn(p: Page, view: string, overlays = ''): Promise<void> {
  await p.evaluate(([v, o]) => {
    const w = window as any;
    w.__name = (fn: unknown) => fn;
    delete w.__bddOperability;
    document.body.innerHTML = `<div id="view">${v}</div>${o}`;
    const root = document.getElementById('view');
    w.grok = {shell: {v: {type: 'datagrok', name: 'Home', root}, views: [], windows: {},
      closeAll: () => { root!.innerHTML = ''; }}};
    const balloons: unknown[] = [];
    w.__bdd = {balloons, takeBalloons: () => balloons.splice(0, balloons.length), putBalloons: (b: unknown[]) => { balloons.unshift(...b); }};
    const scope = w.DG.U2.Scope.root('bdd test');
    for (const e of document.querySelectorAll('[data-u2]'))
      new w.DG.U2.Control(e, scope);
  }, [view, overlays]);
}

const STYLE = `<style>
  body { margin: 0; font: 14px sans-serif; }
  [data-u2="dialog"] { position: fixed; left: 100px; top: 100px; width: 420px; padding: 8px; background: white; z-index: 10; }
  .u2-input-root { display: flex; align-items: center; height: 28px; }
  .u2-input-label { width: 100px; }
  .u2-input-box { display: flex; align-items: center; }
  .u2-input-editor { width: 200px; }
  .u2-input-options { display: flex; }
  .u2-input-options > button { width: 20px; height: 20px; }
</style>`;

// the old read-only rail rule of u2's inputs.css: every option on a read-only input's rail inert
const OLD_RAIL = `<style>.u2-input-readonly .u2-input-options,
.u2-input-readonly .u2-input-clear,
.u2-input-readonly input[type='range'] {
  pointer-events: none;
}</style>`;
const RAIL = `<style>.u2-input-readonly .u2-input-options > .u2-input-edit-option,
.u2-input-readonly .u2-input-clear,
.u2-input-readonly input[type='range'] {
  pointer-events: none;
}</style>`;

/** The Submit dialog's Name row as u2 builds it: a read-only text input with the caller's pencil and download. */
const SUBMIT = '<div data-u2="dialog" role="dialog" aria-label="Submit"><div class="u2-dialog-title"><span class="u2-dialog-title-text">Submit</span></div>' +
  '<div data-u2="text-input" class="u2-input-root u2-input-readonly">' +
  '<label data-u2-part="label" class="u2-input-label">Name</label><div class="u2-input-box">' +
  '<input data-u2-part="editor" class="u2-input-editor" readonly value="KHL-3"><div data-u2-part="options" class="u2-input-options">' +
  '<button aria-label="Edit"></button><button aria-label="Download"></button></div></div></div></div>';

const operability = async (p: Page) => (await auditScenario(p, () => []))
  .findings.filter((f) => f.audit === 'operability').map((f) => `${f.rule}: ${f.detail}`);

/** The reasons failing claims report, failing fast: the claims poll, and a stand-in page never changes. */
async function reasonsFor(...bodies: (() => Promise<unknown>)[]): Promise<string[]> {
  return whileExpectedToFail(() => Promise.all(bodies.map(async (body) => {
    try {
      await body();
    }
    catch (e) {
      return reasonOf(e);
    }
    return 'passed';
  })));
}

scenario('the old read-only rail rule leaves the caller\'s pencil and download without a pointer: the audit and the ' +
  'claim name both; under today\'s rule they are operable', async () => {
  const p = page!;
  await standIn(p, '<div data-u2="button">Home</div>', STYLE + OLD_RAIL + SUBMIT);
  assert.deepEqual(await operability(p), [
    'no-pointer: button "Edit" takes no pointer: pointer-events: none on div.u2-input-options — in dialog "Submit"',
    'no-pointer: button "Download" takes no pointer: pointer-events: none on div.u2-input-options — in dialog "Submit"']);
  const [reason] = await reasonsFor(() => everyControlOperable(p, el('"Name" text input')));
  assert.match(reason, /^2 control\(s\) in "Name" text input cannot be operated as offered: no-pointer: button "Edit" takes no pointer: pointer-events: none on div\.u2-input-options; no-pointer: button "Download"/);

  await standIn(p, '<div data-u2="button">Home</div>', STYLE + RAIL + SUBMIT);
  assert.deepEqual(await operability(p), []);
  await everyControlOperable(p, el('"Name" text input'));
  await everyControlOperable(p, el('"Submit" dialog'));
});

scenario('a covered control is named with what covers it; a balloon over it, or a layer opened after its own, is no finding, ' +
  'and an older layer over a newer one is', async () => {
  const p = page!;
  const run = '<div data-u2="form" aria-label="Run form" style="position: relative; height: 80px">' +
    '<button style="position: absolute; left: 10px; top: 10px; width: 80px; height: 24px">Run</button>';
  await standIn(p, `${run}<div class="glass" style="position: absolute; inset: 0"></div></div>`);
  assert.deepEqual(await operability(p), ['covered: button "Run" is covered by div.glass — in view "Home"']);
  const [reason] = await reasonsFor(() => everyControlOperable(p, el('"Run form" form')));
  assert.equal(reason, '1 control(s) in "Run form" form cannot be operated as offered: covered: button "Run" is covered by div.glass');

  await standIn(p, `${run}</div>`, '<div class="d4-balloon" style="position: fixed; left: 0; top: 0; width: 300px; height: 60px"></div>');
  assert.deepEqual(await operability(p), [], 'a balloon goes away on its own');

  const dialog = '<div data-u2="dialog" aria-label="Pick" style="position: fixed; left: 0; top: 0; width: 300px; height: 100px; ' +
    'z-index: 10; background: white"><button class="pick-choose" style="width: 80px; height: 24px">Choose</button></div>';
  const popup = (z: number) => `<div data-u2-part="popup" data-u2-owner="Stage" style="position: fixed; left: 0; top: 0; ` +
    `width: 200px; height: 60px; z-index: ${z}; background: white"><button style="width: 80px; height: 24px">Discovery</button></div>`;
  await standIn(p, '', dialog);
  assert.equal((await readOperability(p, 'new')).probed, 1, 'the dialog, seen after the step that opened it');
  await p.evaluate((html) => document.body.insertAdjacentHTML('beforeend', html), popup(20));
  assert.deepEqual(await operability(p), [], 'the popup the dialog opened lies over the dialog\'s button');
  await standIn(p, '', dialog);
  await readOperability(p, 'new');
  await p.evaluate((html) => document.body.insertAdjacentHTML('beforeend', html), popup(5));
  assert.deepEqual(await operability(p), ['covered: button "Discovery" is covered by button.pick-choose — in popup of "Stage"']);
});

scenario('a disabled button that says no reason is a finding and fails its claim; one with a reason passes it', async () => {
  const p = page!;
  await standIn(p, '', '<div data-u2="dialog" aria-label="Settings"><button data-u2="button" aria-disabled="true">Save</button>' +
    '<button data-u2="button" aria-disabled="true" aria-description="Name the campaign first">Submit</button>' +
    '<button data-u2="button">Cancel</button></div><button disabled title="Pick a column first">OK</button>');
  assert.deepEqual(await operability(p), ['disabled-without-reason: button "Save" is disabled and says no reason — in dialog "Settings"']);
  await shouldBeDisabledBecause(p, el('Submit button'), 'Name the campaign first');
  await shouldBeDisabledBecause(p, el('OK button'), 'Pick a column first');
  assert.deepEqual(await reasonsFor(
    () => shouldBeDisabledBecause(p, el('Save button'), 'Nothing changed'),
    () => shouldBeDisabledBecause(p, el('Submit button'), 'Pick a template'),
    () => shouldBeDisabledBecause(p, el('Cancel button'), 'Busy'),
  ), [
    'Save button is disabled and says no reason, expected "Nothing changed"',
    'Submit button is disabled because "Name the campaign first", expected "Pick a template"',
    'Cancel button is enabled',
  ]);
});

scenario('a busy control is judged as soon as its aria-busy clears; one still busy after 2 s is left out', async () => {
  const p = page!;
  const busy = '<div data-u2="form" aria-label="Results" aria-busy="true"><button data-u2="button" aria-disabled="true">Export</button></div>';
  await standIn(p, busy);
  const started = Date.now();
  await p.evaluate(() => setTimeout(() => document.querySelector('[aria-label="Results"]')!.removeAttribute('aria-busy'), 300));
  assert.deepEqual(await operability(p), ['disabled-without-reason: button "Export" is disabled and says no reason — in view "Home"']);
  const took = Date.now() - started;
  assert.ok(took >= 300 && took < 1500, `judged once it cleared, not at the deadline: ${took} ms`);

  await standIn(p, busy);
  const again = Date.now();
  assert.deepEqual(await operability(p), [], 'still busy: nothing in it is judged');
  assert.ok(Date.now() - again >= 1900, 'it waited the 2 s');
});

const FILTER_BALLOON = {type: 'warning', message: 'Error adding filter: <br><pre>{\n  "type": "histogram",\n  "column": "HeavyAtomCount"\n}</pre>'};

scenario('an "Error adding filter" balloon no step claimed is a balloons finding; claimed by a balloon step, it is not, ' +
  'and neither is an info balloon', async () => {
  const p = page!;
  const shown = (balloons: object[]) => p.evaluate((b) => { (window as any).__bdd.balloons.push(...b); }, balloons);
  const unclaimed = async () => (await auditScenario(p, () => [])).findings.filter((f) => f.audit === 'balloons')
    .map((f) => `${f.rule}: ${f.detail}`);
  await standIn(p, '<div data-u2="button">Campaign</div>');
  await shown([FILTER_BALLOON, {type: 'info', message: 'Campaign saved successfully.'}]);
  assert.deepEqual(await unclaimed(), ['balloon: warning balloon "Error adding filter: { "type": "histogram", "column": "HeavyAtomCount" }" ' +
    'that no step claimed']);
  await shown([FILTER_BALLOON]);
  await errorOrWarningBalloonMatching(p, '^Error adding filter');
  assert.deepEqual(await unclaimed(), []);
});

test('an allowed balloon is matched by its text without markup', () => {
  const balloons = [FILTER_BALLOON, {type: 'error', message: 'Campaign KHL-9 could not be opened'}];
  assert.deepEqual(balloonFindings(balloons, [{message: '^Error adding filter: \\{ "type"', reason: 'r'}]).map((f) => f.detail),
    ['error balloon "Campaign KHL-9 could not be opened" that no step claimed']);
});

/** The harness hooks over the stand-in shell: a fake `test` that records the hooks and the test's info. */
function fakeTest(info: {title: string; status: string; expectedStatus: string; tags: string[]; annotations: {type: string; description?: string}[]}) {
  const hooks: {afterEach?: () => Promise<void>} = {};
  const api = {afterEach: (fn: () => Promise<void>) => { hooks.afterEach = fn; }, afterAll: () => undefined, info: () => ({...info, timeout: 30000}),
    setTimeout: () => undefined, step: async (_: string, body: () => Promise<void>) => body()};
  return {api: api as unknown as Parameters<typeof feature>[0], hooks};
}

scenario('the harness probes a dialog after the step that opened it, so a dialog the scenario closed again is audited; ' +
  'in fail mode that fails the scenario', async () => {
  const info = {title: 'Submit a campaign', status: 'passed', expectedStatus: 'passed', tags: [], annotations: [] as {type: string; description?: string}[]};
  const {api, hooks} = fakeTest(info);
  const session = feature(api);
  const p = await session.page(browser!);
  await p.addInitScript(() => { (window as any).__name = (fn: unknown) => fn; });
  await withEngine(p);
  const open = () => standIn(p, '<div data-u2="button">Campaign</div>', STYLE + OLD_RAIL + SUBMIT);
  const close = () => p.evaluate(() => document.querySelector('[data-u2="dialog"]')!.remove());
  const finding = 'no-pointer: button "Edit" takes no pointer: pointer-events: none on div.u2-input-options — in dialog "Submit"';
  for (const root of [REPORT, FAIL]) {
    process.env.BDD_ROOT = root;
    info.annotations.length = 0;
    await standIn(p, '<div data-u2="button">Campaign</div>');
    await session.page(browser!);
    await session.step(5, 'When user opens the Submit dialog', open);
    await session.step(6, 'And user closes the Submit dialog', close);
    const warn = console.warn;
    console.warn = () => undefined;
    try {
      if (root === REPORT)
        await hooks.afterEach!();
      else
        await assert.rejects(hooks.afterEach!(), (e: Error) => e.message.includes('audit: 2 finding(s) after "Submit a campaign"') &&
          e.message.includes(finding));
    }
    finally {
      console.warn = warn;
    }
    assert.ok(info.annotations.some((a) => a.type === 'bdd-audits' && a.description!.split(', ').includes('operability')));
    assert.ok(info.annotations.some((a) => a.type === 'bdd-audit' && a.description === finding), JSON.stringify(info.annotations));
  }
});

scenario('with no audits in its config a scenario fails on an error balloon no step claimed, and on a control off with ' +
  'no reason; the platform\'s grid ribbon icons are accepted', async () => {
  process.env.BDD_ROOT = PLAIN;
  const info = {title: 'Open a campaign', status: 'passed', expectedStatus: 'passed', tags: [], annotations: [] as {type: string; description?: string}[]};
  const {api, hooks} = fakeTest(info);
  const session = feature(api);
  const p = await session.page(browser!);
  await p.addInitScript(() => { (window as any).__name = (fn: unknown) => fn; });
  await withEngine(p);
  const cases: [string, () => Promise<unknown>, string][] = [
    ['an error balloon', () => p.evaluate(() => (window as any).__bdd.balloons.push({type: 'error', message: 'Campaign KHL-9 could not be opened'})),
      'balloon: error balloon "Campaign KHL-9 could not be opened" that no step claimed'],
    ['a reasonless off control', () => standIn(p, '<div data-u2="form" aria-label="Settings"><button data-u2="button" aria-disabled="true">Save</button></div>'),
      'disabled-without-reason: button "Save" is disabled and says no reason — in view "Home"'],
  ];
  const warn = console.warn;
  console.warn = () => undefined;
  try {
    for (const [what, show, finding] of cases) {
      info.annotations.length = 0;
      await standIn(p, '<div data-u2="button">Campaign</div>');
      await session.page(browser!);
      await session.step(5, `When user meets ${what}`, show);
      await assert.rejects(hooks.afterEach!(), (e: Error) => e.message.includes('audit: 1 finding(s) after "Open a campaign"') &&
        e.message.includes(finding), what);
    }
  }
  finally {
    console.warn = warn;
  }

  const icon = (cls: string, name: string) => `<i class="grok-icon svg-icon ${cls}" role="button" tabindex="0" aria-disabled="true" ` +
    `aria-label="${name}"></i>`;
  await standIn(p, `<div data-u2="toolbar">${icon('svg-remove-selected-rows', 'Remove selected rows')}` +
    `${icon('svg-remove-selected-columns', 'Remove selected columns')}${icon('svg-add-new-row', 'Add new row')}</div>`);
  assert.deepEqual(await operability(p), ['disabled-without-reason: i "Add new row" is disabled and says no reason — in view "Home"']);
});
