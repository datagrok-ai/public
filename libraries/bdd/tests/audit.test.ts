/* The automatic audits: modes and allowlists, which scenarios are audited, the findings and how they are
   settled, and — in the library's Chromium, over a stand-in shell and a fake `DG.U2.debug` — the leak
   baseline and check, the scenario audit and the harness wiring (audits before the reset, leaks after it,
   a journey's scenarios one by one). No stand. */
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';
import {mkdirSync, mkdtempSync, rmSync, writeFileSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {after, before, test} from 'node:test';
import {type Browser, chromium, type Page} from '@playwright/test';

import {auditModes, auditScenario, BUILTIN_ALLOW, errorFindings, leakCheck, leakFindings, markLeakBaseline, resolveAllow,
  resolveModes, settleAudit, shouldAudit, strictU2, trackLeaks} from '../src/runtime/audit.js';
import {feature, journey} from '../src/runtime/harness.js';

for (const key of ['BDD_AUDIT', 'BDD_AUDIT_LAYOUT', 'BDD_AUDIT_ERRORS', 'BDD_AUDIT_BALLOONS', 'BDD_AUDIT_OPERABILITY', 'BDD_AUDIT_LEAKS'])
  delete process.env[key];

const tmp = mkdtempSync(join(tmpdir(), 'bdd-audit-'));
function project(name: string, audits: object): string {
  const dir = join(tmp, name);
  mkdirSync(dir, {recursive: true});
  writeFileSync(join(dir, 'bdd.config.json'), JSON.stringify({audits}));
  return dir;
}
const REPORT = project('report', {layout: 'report', errors: 'report', balloons: 'report', operability: 'report', leaks: 'report',
  allow: {leaks: [{label: '^tooltip$', reason: 'the tooltip service lives with the page'}]}});
const FAIL = project('fail', {layout: 'fail', errors: 'report', leaks: 'fail'});
process.env.BDD_ROOT = REPORT;

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
  rmSync(tmp, {recursive: true, force: true});
});

/** A stand-in shell (`grok.shell` with a Home view) and a fake `DG.U2.debug` over `records`, each live until
 * its `until` (performance.now() ms) if it has one. `open` adds a record, the way a view's control would. */
async function standIn(p: Page, html = '', debug: 'tracking' | 'off' | 'missing' = 'tracking'): Promise<void> {
  await p.setContent(`<div id="view">${html}</div>`);
  await p.evaluate((mode) => {
    const w = window as any;
    w.__name = (fn: unknown) => fn;
    const root = document.getElementById('view');
    w.grok = {shell: {v: {type: 'datagrok', root}, views: [], windows: {},
      closeAll: () => { root!.innerHTML = ''; w.__closed = (w.__closed ?? 0) + 1; }}};
    const records: any[] = [];
    w.__records = records;
    if (mode === 'missing')
      return;
    w.DG = {U2: {debug: {
      enabled: mode === 'tracking',
      swept: 0,
      snapshot: () => ({ids: records.map((r) => r.id)}),
      sweep: () => { w.DG.U2.debug.swept++; return []; },
      diff: (since: {ids: number[]}) => {
        const live = records.filter((r) => !since.ids.includes(r.id) && (r.until === undefined || performance.now() < r.until));
        return {scopes: live.filter((r) => r.kind === 'scope'), timers: live.filter((r) => r.kind === 'timer'),
          listeners: live.filter((r) => r.kind === 'listener'), total: live.filter((r) => r.label !== 'singleton').length};
      },
    }}};
  }, debug);
}

const open = (p: Page, record: object, liveMs?: number): Promise<void> => p.evaluate(([r, ms]) => {
  const records = (window as any).__records;
  records.push({id: records.length + 100, ...r, ...(ms === undefined ? {} : {until: performance.now() + ms})});
}, [record, liveMs] as [object, number | undefined]);

const STACK = 'Error\n    at new TextInput (https://host/api/packages/published/files/U2Demo/x/dist/package.js:2:345)\n' +
  '    at textInput (https://host/dist/package.js:2:900)\n    at basicInputs (https://host/dist/package.js:9:12)\n    at four (x.js:1:1)';

test('modes: per-audit env over BDD_AUDIT over the project config over the defaults', () => {
  assert.deepEqual(resolveModes(undefined, {}), {layout: 'fail', errors: 'fail', balloons: 'fail', operability: 'fail', leaks: 'off'});
  assert.deepEqual(resolveModes({leaks: 'report', layout: 'report', errors: 'off'}, {}),
    {layout: 'report', errors: 'off', balloons: 'fail', operability: 'fail', leaks: 'report'});
  assert.deepEqual(resolveModes({leaks: 'report', layout: 'fail'}, {BDD_AUDIT: 'off'}),
    {layout: 'off', errors: 'off', balloons: 'off', operability: 'off', leaks: 'off'});
  assert.deepEqual(resolveModes({layout: 'fail'}, {BDD_AUDIT: 'off', BDD_AUDIT_LEAKS: 'fail', BDD_AUDIT_OPERABILITY: 'report', BDD_AUDIT_BALLOONS: 'fail'}),
    {layout: 'off', errors: 'off', balloons: 'fail', operability: 'report', leaks: 'fail'});
  assert.throws(() => resolveModes(undefined, {BDD_AUDIT_LAYOUT: 'strict'}), /BDD_AUDIT_LAYOUT: "strict" is not an audit mode \(fail, report, off\)/);
  assert.throws(() => resolveModes({errors: 'warn' as never}, {}), /bdd.config.json audits.errors: "warn"/);
});

test('modes are read from the project once per worker; a project without audits gets the defaults', () => {
  assert.deepEqual(auditModes(REPORT), {layout: 'report', errors: 'report', balloons: 'report', operability: 'report', leaks: 'report'});
  assert.deepEqual(auditModes(FAIL), {layout: 'fail', errors: 'report', balloons: 'fail', operability: 'fail', leaks: 'fail'});
  const plain = join(tmp, 'plain');
  mkdirSync(plain);
  assert.deepEqual(auditModes(plain), {layout: 'fail', errors: 'fail', balloons: 'fail', operability: 'fail', leaks: 'off'});
  writeFileSync(join(plain, 'bdd.config.json'), JSON.stringify({audits: {layout: 'off'}}));
  assert.equal(auditModes(plain).layout, 'fail');
});

test('every allowed finding says why; its patterns are regular expressions; the library\'s come first', () => {
  assert.throws(() => resolveAllow({allow: {errors: [{message: 'x', reason: ' '}]}}), /audits.allow.errors .*says why \(reason\)/);
  assert.throws(() => resolveAllow({allow: {leaks: [{label: '(', reason: 'r'}]}}), /label is not a regular expression/);
  const entry = {rule: 'overlap' as const, selector: '.x', reason: 'a badge'};
  BUILTIN_ALLOW.layout.push({rule: 'zero-size', selector: '.y', reason: 'library'});
  try {
    assert.deepEqual(resolveAllow({allow: {layout: [entry]}}).layout.map((e) => e.selector), ['.y', '.x']);
  }
  finally {
    BUILTIN_ALLOW.layout.pop();
  }
});

test('a scenario is audited when it passed as expected and carries no @known-failure', () => {
  const info = (status: string, expectedStatus = 'passed', tags: string[] = []) => ({status, expectedStatus, tags}) as never;
  assert.equal(shouldAudit(info('passed')), true);
  assert.equal(shouldAudit(info('passed', 'passed', ['@demo'])), true);
  assert.equal(shouldAudit(info('failed')), false);
  assert.equal(shouldAudit(info('timedOut')), false);
  assert.equal(shouldAudit(info('skipped', 'skipped')), false);
  assert.equal(shouldAudit(info('passed', 'failed')), false);
  assert.equal(shouldAudit(info('passed', 'passed', ['@known-failure'])), false);
  assert.equal(shouldAudit(info('passed', 'passed', ['@journey', '@known-failure:12'])), false);
});

test('a console error is a finding unless an entry matches its text and, when it names one, where it came from', () => {
  const logged = [{message: 'TypeError: x is null\n    at f (a.js:1:1)', url: 'https://host/dist/package.js'},
    {message: 'ResizeObserver loop limit exceeded', url: ''}, {message: 'ResizeObserver loop limit exceeded', url: 'https://other/x.js'}];
  const found = errorFindings(logged, [{message: '^ResizeObserver loop', url: '^$', reason: 'browser noise'}]);
  assert.deepEqual(found.map((f) => f.detail), ['TypeError: x is null (https://host/dist/package.js)',
    'ResizeObserver loop limit exceeded (https://other/x.js)']);
  assert.ok(found.every((f) => f.audit === 'errors' && f.rule === 'console-error'));
});

test('a leaked record is one sentence with its first three frames; an entry matching its label or stack drops it', () => {
  const read = {singletons: 1,
    scopes: [{id: 812, label: 'text-input', name: 'org', stack: STACK}, {id: 9, label: 'tooltip'}],
    timers: [{kind: 'interval', ms: 500, stack: STACK}],
    listeners: [{target: 'window', type: 'keydown', capture: true}]};
  assert.deepEqual(leakFindings(read, []).map((f) => `${f.rule}: ${f.detail}`), [
    'leak: scope #812 text-input (org) still live after the view closed — created at new TextInput (package.js:2:345) < ' +
      'textInput (package.js:2:900) < basicInputs (package.js:9:12)',
    'leak: scope #9 tooltip still live after the view closed',
    'leak: interval 500 ms timer still pending after the view closed — created at new TextInput (package.js:2:345) < ' +
      'textInput (package.js:2:900) < basicInputs (package.js:9:12)',
    'leak: window keydown capture listener still attached after the view closed',
  ]);
  const allowed = leakFindings(read, [{label: '^tooltip$', reason: 'page-wide'}, {stack: 'basicInputs', label: 'timer', reason: 'demo'},
    {label: 'window keydown', reason: 'overlay'}]);
  assert.deepEqual(allowed.map((f) => f.detail.split(' still')[0]), ['scope #812 text-input (org)']);
});

test('report mode annotates and prints; fail mode answers the error that fails the scenario; each audit by its own mode', () => {
  const warned: string[] = [];
  const warn = console.warn;
  console.warn = (m: string) => void warned.push(m);
  try {
    const outcome = {ran: ['layout', 'errors', 'leaks'] as never[], findings: [
      {audit: 'layout' as const, rule: 'overlap' as const, detail: '"Name" text input overlaps "Notes" text input by 150×8 px at (0, 20)'},
      {audit: 'errors' as const, rule: 'console-error' as const, detail: 'TypeError: x'},
      {audit: 'leaks' as const, rule: 'leak' as const, detail: 'scope #3 dialog still live after the view closed'}]};
    process.env.BDD_ROOT = REPORT;
    const reported = {annotations: [] as {type: string; description?: string}[]};
    assert.equal(settleAudit(outcome, 'Feature › Scenario', reported), undefined);
    assert.deepEqual(reported.annotations, [{type: 'bdd-audits', description: 'layout, errors, leaks'},
      {type: 'bdd-audit', description: 'overlap: "Name" text input overlaps "Notes" text input by 150×8 px at (0, 20)'},
      {type: 'bdd-audit', description: 'console-error: TypeError: x'},
      {type: 'bdd-audit', description: 'leak: scope #3 dialog still live after the view closed'}]);
    assert.equal(warned.length, 1);
    assert.match(warned[0], /^bdd: audit \(report\) Feature › Scenario: 3 finding\(s\)\n {2}overlap: "Name"/);

    process.env.BDD_ROOT = FAIL;
    warned.length = 0;
    const failed = {annotations: [] as {type: string; description?: string}[]};
    const error = settleAudit(outcome, 'Journey › Step two', failed, 'Step two');
    assert.equal(error?.message, 'audit: 2 finding(s) after "Journey › Step two"\n' +
      '  overlap: "Name" text input overlaps "Notes" text input by 150×8 px at (0, 20)\n' +
      '  leak: scope #3 dialog still live after the view closed');
    assert.deepEqual(warned, ['bdd: audit (report) Journey › Step two: 1 finding(s)\n  console-error: TypeError: x']);
    assert.equal(failed.annotations[1].description, 'overlap: "Name" text input overlaps "Notes" text input by 150×8 px at (0, 20) (scenario "Step two")');
    assert.equal(settleAudit({ran: [], findings: []}, 'quiet', failed), undefined);
    assert.equal(failed.annotations.length, 4);
  }
  finally {
    console.warn = warn;
    process.env.BDD_ROOT = REPORT;
  }
});

scenario('tracking is switched on before the page\'s own scripts run, and only where leaks are audited', async () => {
  const fresh = await browser!.newPage();
  const quiet = await browser!.newPage();
  try {
    await trackLeaks(fresh);
    await trackLeaks(fresh);
    await fresh.goto('data:text/html,<script>document.title = String(globalThis.__U2_DEBUG__)</script>');
    assert.equal(await fresh.title(), 'true');
    process.env.BDD_ROOT = join(tmp, 'plain');
    await trackLeaks(quiet);
    await quiet.goto('data:text/html,<script>document.title = String(globalThis.__U2_DEBUG__)</script>');
    assert.equal(await quiet.title(), 'undefined');
  }
  finally {
    await fresh.close();
    await quiet.close();
  }
});

scenario('u2 throws on a misuse in every page, set before the page\'s own scripts run; U2_STRICT names another mode', async () => {
  const fresh = await browser!.newPage();
  const named = await browser!.newPage();
  const was = process.env.U2_STRICT;
  const probe = 'data:text/html,<script>document.title = String(globalThis.__U2_STRICT__)</script>';
  try {
    delete process.env.U2_STRICT;
    await strictU2(fresh);
    await strictU2(fresh);
    await fresh.goto(probe);
    assert.equal(await fresh.title(), 'throw');
    process.env.U2_STRICT = 'report';
    await strictU2(named);
    await named.goto(probe);
    assert.equal(await named.title(), 'report');
  }
  finally {
    if (was === undefined)
      delete process.env.U2_STRICT;
    else
      process.env.U2_STRICT = was;
    await fresh.close();
    await named.close();
  }
});

scenario('a record still live after the view closed is a leak; one that clears within 2 s, or a singleton, is not', async () => {
  await standIn(page!);
  await markLeakBaseline(page!);
  await open(page!, {kind: 'scope', label: 'text-input', name: 'org', stack: STACK});
  await open(page!, {kind: 'timer', label: undefined, ms: 300}, 300);
  await open(page!, {kind: 'listener', target: 'document', type: 'keydown', capture: false, label: 'singleton'});
  const started = Date.now();
  const outcome = await leakCheck(page!);
  assert.deepEqual(outcome.ran, ['leaks']);
  assert.deepEqual(outcome.findings.map((f) => f.detail.split(' — ')[0]), ['scope #100 text-input (org) still live after the view closed']);
  assert.ok(Date.now() - started >= 1900, 'the check waits up to 2 s while anything is live');
  assert.ok(await page!.evaluate(() => (window as any).DG.U2.debug.swept) > 1, 'the orphan sweep runs before every read');

  await markLeakBaseline(page!);
  await open(page!, {kind: 'timer', label: undefined, kind2: 'debounce'}, 300);
  const fast = Date.now();
  assert.deepEqual(await leakCheck(page!), {ran: ['leaks'], findings: []});
  assert.ok(Date.now() - fast < 1500, 'nothing live: the check ends as soon as the debounce cleared');
  assert.deepEqual(await leakCheck(page!), {ran: [], findings: []}, 'a baseline is used once');
});

scenario('an allowed leak is dropped; no tracking or a new document means no check', async () => {
  await standIn(page!);
  await markLeakBaseline(page!);
  await open(page!, {kind: 'scope', label: 'tooltip'});
  assert.deepEqual(await leakCheck(page!), {ran: ['leaks'], findings: []});

  await markLeakBaseline(page!);
  await page!.goto('about:blank');
  await standIn(page!);
  assert.deepEqual(await leakCheck(page!), {ran: [], findings: []}, 'a navigation drops the baseline with its document');

  const warned: string[] = [];
  const warn = console.warn;
  console.warn = (m: string) => void warned.push(m);
  try {
    await standIn(page!, '', 'off');
    await markLeakBaseline(page!);
    await markLeakBaseline(page!);
    assert.deepEqual(await leakCheck(page!), {ran: [], findings: []});
    assert.equal(warned.length, 1);
    assert.match(warned[0], /DG.U2.debug is not tracking .*: no leak check on this worker/);
  }
  finally {
    console.warn = warn;
  }
});

const OVERLAPPING = '<style>[data-u2] { display: block; width: 200px; height: 28px; }</style>' +
  '<div data-u2="text-input" aria-label="Name">Name</div><div data-u2="text-input" aria-label="Notes" style="margin-top: -20px">Notes</div>';

scenario('a u2 view is audited for layout and console errors; a view without u2 is not, and its errors stay', async () => {
  await standIn(page!, OVERLAPPING);
  let taken = 0;
  const logged = () => {
    taken++;
    return [{message: 'TypeError: boom', url: 'https://host/dist/package.js'}];
  };
  const outcome = await auditScenario(page!, logged);
  assert.deepEqual(outcome.ran, ['layout', 'errors', 'balloons']);
  assert.deepEqual(outcome.findings.map((f) => f.rule), ['overlap', 'console-error']);
  assert.match(outcome.findings[0].detail, /"Name" text input overlaps "Notes" text input/);
  assert.equal(taken, 1);

  await standIn(page!, '<div class="d4-grid">a Dart grid</div>');
  assert.deepEqual(await auditScenario(page!, logged), {ran: [], findings: []});
  assert.equal(taken, 1, 'errors a platform page logs are left to its own floor');
});

/** The harness hooks over a stand-in shell: a fake `test` that records the hooks and the test's info. */
function fakeTest(info: {title: string; status: string; expectedStatus: string; tags: string[]; annotations: {type: string; description?: string}[]}) {
  const hooks: {afterEach?: () => Promise<void>} = {};
  const api = {afterEach: (fn: () => Promise<void>) => { hooks.afterEach = fn; }, afterAll: () => undefined, info: () => ({...info, timeout: 30000}),
    setTimeout: () => undefined, step: async (_: string, body: () => Promise<void>) => body()};
  return {api: api as unknown as Parameters<typeof feature>[0], hooks};
}

scenario('the harness audits the final state before the shell reset and the leaks after it', async () => {
  const info = {title: 'Basic inputs', status: 'passed', expectedStatus: 'passed', tags: [], annotations: [] as {type: string; description?: string}[]};
  const {api, hooks} = fakeTest(info);
  const session = feature(api);
  const p = await session.page(browser!);
  await standIn(p, OVERLAPPING);
  await session.page(browser!);
  // a control of the view that closeAll does not dispose, and one that it does
  await open(p, {kind: 'scope', label: 'text-input', name: 'org', stack: STACK});
  await p.evaluate(() => {
    const w = window as any;
    w.__records.push({id: 999, kind: 'scope', label: 'form'});
    const close = w.grok.shell.closeAll;
    w.grok.shell.closeAll = () => { close(); w.__records.pop(); };
  });
  const warn = console.warn;
  console.warn = () => undefined;
  try {
    await hooks.afterEach!();
  }
  finally {
    console.warn = warn;
  }
  assert.equal(await p.evaluate(() => (window as any).__closed), 1);
  assert.deepEqual(info.annotations.map((a) => a.description?.split(' — ')[0]), ['layout, errors, balloons, leaks',
    info.annotations[1].description, 'leak: scope #100 text-input (org) still live after the view closed']);
  assert.match(info.annotations[1].description!, /^overlap: "Name" text input overlaps "Notes"/);

  info.status = 'failed';
  info.annotations.length = 0;
  await standIn(p, OVERLAPPING);
  await session.page(browser!);
  await open(p, {kind: 'scope', label: 'text-input'});
  await hooks.afterEach!();
  assert.deepEqual(info.annotations, [], 'a failed scenario is not audited');
});

scenario('a journey audits each scenario as it ends and fails the one with findings in fail mode; the leaks once, at its end', async () => {
  process.env.BDD_ROOT = FAIL;
  const info = {title: 'Journey', status: 'passed', expectedStatus: 'passed', tags: ['@journey'], annotations: [] as {type: string; description?: string}[]};
  const {api, hooks} = fakeTest(info);
  const session = feature(api);
  const p = await session.page(browser!);
  await standIn(p, OVERLAPPING);
  await session.page(browser!);
  const run = journey(api, 2, p);
  await run.scenario('Overlapping', async () => undefined);
  await p.evaluate(() => { document.getElementById('view')!.innerHTML = '<div data-u2="button">OK</div>'; });
  await run.scenario('Clean', async () => undefined);
  assert.throws(() => run.finish(), (e: Error) => /Overlapping/.test(e.message) && /audit: 1 finding\(s\) after "Journey › Overlapping"/.test(e.message) &&
    !/Clean/.test(e.message.split('\n')[0]));
  assert.deepEqual(info.annotations.map((a) => a.description?.replace(/: ".*/, '')), ['layout, errors, balloons (scenario "Overlapping")', 'overlap',
    'layout, errors, balloons (scenario "Clean")']);
  info.annotations.length = 0;
  await open(p, {kind: 'listener', target: 'window', type: 'resize', capture: false});
  await assert.rejects(hooks.afterEach!(), /audit: 1 finding\(s\) after "Journey"\n {2}leak: window resize listener still attached/);
  assert.deepEqual(info.annotations.map((a) => a.type), ['bdd-audits', 'bdd-audit'], 'no second layout audit at the journey\'s end');
});

test('audit-summary counts findings by rule, feature and scenario, and reads the memory log', () => {
  const dir = join(tmp, 'summary');
  mkdirSync(dir);
  const result = (annotations: object[]) => ({status: 'passed', annotations});
  const spec = (file: string, title: string, annotations: object[]) => ({file, title, tests: [{status: 'expected', annotations: [], results: [result(annotations)]}]});
  writeFileSync(join(dir, 'report.json'), JSON.stringify({stats: {expected: 3, unexpected: 0, skipped: 0, flaky: 0, duration: 4000}, suites: [
    {file: 'demo/a.test.ts', specs: [spec('demo/a.test.ts', 'One', [{type: 'bdd-audits', description: 'layout, errors, leaks'},
      {type: 'bdd-audit', description: 'overlap: "A" overlaps "B"'}, {type: 'bdd-audit', description: 'leak: scope #4 dialog still live after the view closed'},
      {type: 'bdd-audit', description: 'leak: window keydown listener still attached after the view closed'}]),
    spec('demo/a.test.ts', 'Two', [{type: 'bdd-audits', description: 'layout, errors, balloons, operability'},
      {type: 'bdd-audit', description: 'no-pointer: button "Edit" takes no pointer: pointer-events: none on div.u2-input-options — in dialog "Submit"'},
      {type: 'bdd-audit', description: 'balloon: warning balloon "Error adding filter" that no step claimed'}])]},
    {file: 'demo/b.test.ts', specs: [spec('demo/b.test.ts', 'J', [{type: 'bdd-audits', description: 'layout, errors (scenario "S1")'},
      {type: 'bdd-audit', description: 'leak: scope #9 dialog still live after the view closed (scenario "S1")'}])]},
  ]}));
  const line = (worker: string, onPage: number, heap: number) => JSON.stringify({worker, onPage, page: {heapUsedMb: heap},
    afterGc: {heapUsedMb: heap - 50}, processes: {renderer: {mb: heap * 4}}, processesAfterGc: {renderer: {mb: heap * 3}}});
  writeFileSync(join(dir, 'memory.jsonl'), [line('0', 1, 150), line('0', 2, 170), line('0', 1, 140), line('1', 1, 160)].join('\n') + '\n');
  const out = execFileSync(process.execPath, [join(import.meta.dirname, '..', 'tool', 'audit-summary.mjs'), join(dir, 'report.json'),
    join(dir, 'memory.jsonl')], {encoding: 'utf8'});
  assert.match(out, /audited 3 test\(s\): layout 3, console errors 3, balloons 1, operability 1, leaks 1; findings 6/);
  assert.match(out, /\n {2}leak {15}3\n {2}balloon {12}1\n {2}no-pointer {9}1\n {2}overlap {12}1\n/);
  assert.match(out, /demo\/a +2 +\| 1 +\| 1 +\| 1 +\| 1, 0, 1/);
  assert.match(out, /demo\/a › Two: 0 \| 1 \| 1 \| 0, 0, 0/);
  assert.match(out, /demo\/b › J › S1: 0 \| 0 \| 0 \| 1, 0, 0/);
  assert.match(out, /2× +leak: scope #n dialog still live after the view closed/);
  assert.match(out, /worker 0: 3 features on 2 page\(s\); JS heap 150 → 140 MB, after GC 100 → 90 MB; renderer 600 → 560 MB \(after GC 450 → 420 MB\)/);
});
