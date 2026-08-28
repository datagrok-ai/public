/* The W4 wave (plan-w4.md WO-6): `visible:`/`enabled:` expressions, `validator:` (regex and
   expression), named `validators:` and categoryGroups on the bench's Expressions & validation
   section — one shared `fceW4Vehicle` call under the Dart InputForm and u2 funcForm at once.
   Selectors stay scoped to the convergence page (the inputs bench shares class names) and to
   `[data-row=…]` cells; every wait polls DOM state. The section's status lines are prefixed
   `w4 sync = ` / `w4 inputs = ` so the func-form/async/tables lanes' positional reads of the
   plain lines stay unmoved. Console discipline: the section's `slow` row rejects the Dart form
   build (an async validator — divergence #21) at page open in EVERY run of this fixture; the
   rejection is caught into the cell and probed console-clean (WO-5 driver, 0 errors). */
import {ok, openApp, shot} from '../local.mjs';
import {APP} from '../lib.mjs';

const PAGE = '.u2demo-page:has([data-row="batchId"])';
const WHOLE = `${PAGE} .u2demo-w4-whole`;

const u2Cell = (n) => `${PAGE} .u2-input-root[data-row="${n}"]`;
const dartCell = (n) => `${PAGE} .ui-input-root[data-row="${n}"]`;

const vis = (page, sel) => page.evaluate((sel) => {
  const el = document.querySelector(sel);
  return el == null ? {present: false} : {present: true, hiddenAttr: el.hidden === true,
    displayNone: getComputedStyle(el).display === 'none'};
}, sel);

const w4Status = (page, prefix) => page.evaluate(({sel, prefix}) =>
  [...document.querySelector(sel)?.querySelectorAll('.u2demo-status') ?? []]
    .map((s) => s.textContent).find((t) => t.startsWith(prefix)) ?? '', {sel: PAGE, prefix});

const counters = async (page) => {
  const m = (await w4Status(page, 'w4 sync = ')).match(/u2 → Dart (\d+) · Dart → u2 (\d+)/);
  return m === null ? {toDart: -1, toU2: -1} : {toDart: Number(m[1]), toU2: Number(m[2])};
};

const errText = (page, cell) => page.evaluate((cell) =>
  document.querySelector(`${cell} .u2-input-error`)?.textContent ?? '', cell);

const settle = (page) => page.evaluate(() => new Promise((resolve) => setTimeout(resolve, 250)));

const flipType = async (page, value, waitSel, shown) => {
  await page.locator(`${u2Cell('type')} select`).selectOption(value);
  await page.waitForFunction(({sel, shown}) =>
    (getComputedStyle(document.querySelector(sel)).display === 'none') !== shown,
  {sel: waitSel, shown}, {timeout: 8000});
};

export async function fixture(page) {
  await page.evaluate(() => localStorage.setItem('u2demo.tab', 'func-convergence'));
  await openApp(page, APP.package, 'u2DemoApp');
  await page.waitForSelector(`${PAGE} [data-row="notes"]`, {timeout: 30000});
  await page.waitForSelector(`${WHOLE} [data-u2="func-form"]`, {timeout: 30000});
  await settle(page);
}

async function checkInitialState(page) {
  const battU2 = await vis(page, u2Cell('batteryCapacity'));
  const battDart = await vis(page, dartCell('batteryCapacity'));
  const cylU2 = await vis(page, u2Cell('cylinders'));
  const cylDart = await vis(page, dartCell('cylinders'));
  const extDisabled = await page.evaluate(({u2, dart}) => ({
    u2: document.querySelector(`${u2} input`)?.disabled === true,
    dart: document.querySelector(`${dart} input`)?.disabled === true,
  }), {u2: u2Cell('tankExtension'), dart: dartCell('tankExtension')});
  const inputs = await w4Status(page, 'w4 inputs = ');
  await shot(page, 'func-form-expressions-1-rest');
  ok('func-form-expressions/1a/the-electric-only-field-opens-hidden-on-both-sides-under-the-ice-preset',
    battU2.present && battU2.hiddenAttr && battU2.displayNone && battDart.displayNone &&
    !cylU2.displayNone && !cylDart.displayNone &&
    extDisabled.u2 && extDisabled.dart &&
    inputs.includes("type: 'ICE'") && inputs.includes('tankVolume: 40'),
    `battery u2=${JSON.stringify(battU2)} dart=${JSON.stringify(battDart)} ` +
    `ext disabled=${JSON.stringify(extDisabled)}`);
}

async function checkTypeFlip(page) {
  const before = await counters(page);
  await flipType(page, 'Electric', u2Cell('cylinders'), false);
  await page.waitForFunction(({u2, dart, battU2, battDart}) =>
    getComputedStyle(document.querySelector(u2)).display === 'none' &&
    getComputedStyle(document.querySelector(dart)).display === 'none' &&
    getComputedStyle(document.querySelector(battU2)).display !== 'none' &&
    getComputedStyle(document.querySelector(battDart)).display !== 'none',
  {u2: u2Cell('tankVolume'), dart: dartCell('tankVolume'),
    battU2: u2Cell('batteryCapacity'), battDart: dartCell('batteryCapacity')}, {timeout: 8000});
  await settle(page);
  const after = await counters(page);
  await shot(page, 'func-form-expressions-2-electric');
  ok('func-form-expressions/2a/the-type-flip-hides-engine-and-shows-battery-on-both-columns-one-tick-no-echo',
    after.toDart === before.toDart + 1 && after.toU2 === before.toU2,
    `counters ${JSON.stringify(before)} → ${JSON.stringify(after)}`);
}

async function checkHeaderAutoHide(page) {
  // still Electric from /2 — the whole-form pair: u2's emptied Engine header auto-hides
  // ([hidden] must be VISUAL: computed display none, the inventory-12 regression class);
  // the Dart form's expression path bypasses its own machinery, so its header stays (#19)
  const headers = await page.evaluate((sel) => {
    const whole = document.querySelector(sel);
    const u2Engine = [...whole?.querySelectorAll('[data-u2="func-form"] .u2-form-category') ?? []]
      .find((h) => h.textContent === 'Engine');
    const dartEngine = [...whole?.firstElementChild?.querySelectorAll('h1, h2, h3') ?? []]
      .find((h) => h.textContent.trim() === 'Engine');
    return {
      u2HiddenAttr: u2Engine?.hidden === true,
      u2DisplayNone: u2Engine != null && getComputedStyle(u2Engine).display === 'none',
      dartStays: dartEngine != null && getComputedStyle(dartEngine).display !== 'none',
    };
  }, WHOLE);
  await shot(page, 'func-form-expressions-3-headers');
  ok('func-form-expressions/3a/the-u2-engine-header-hides-attr-and-computed-display-while-darts-stays',
    headers.u2HiddenAttr && headers.u2DisplayNone && headers.dartStays,
    JSON.stringify(headers));
  await flipType(page, 'ICE', u2Cell('cylinders'), true);
  await page.waitForFunction((sel) => {
    const h = [...document.querySelectorAll(`${sel} [data-u2="func-form"] .u2-form-category`)]
      .find((x) => x.textContent === 'Engine');
    return h != null && h.hidden === false && getComputedStyle(h).display !== 'none';
  }, WHOLE, {timeout: 8000});
  ok('func-form-expressions/3b/the-header-returns-when-the-flip-shows-the-engine-fields-again', true,
    'Engine header visible again under ICE');
}

async function checkEnabledFlip(page) {
  await page.locator(`${u2Cell('tankVolume')} input`).fill('60');
  await page.waitForFunction(({u2, dart}) =>
    document.querySelector(`${u2} input`)?.disabled === false &&
    document.querySelector(`${dart} input`)?.disabled === false,
  {u2: u2Cell('tankExtension'), dart: dartCell('tankExtension')}, {timeout: 8000});
  await page.locator(`${u2Cell('tankVolume')} input`).fill('40');
  await page.waitForFunction(({u2, dart}) =>
    document.querySelector(`${u2} input`)?.disabled === true &&
    document.querySelector(`${dart} input`)?.disabled === true,
  {u2: u2Cell('tankExtension'), dart: dartCell('tankExtension')}, {timeout: 8000});
  ok('func-form-expressions/4a/tankvolume-60-enables-tankextension-on-both-sides-and-40-disables-it-back', true,
    'enabled: tankVolume > 50 — both directions, both columns');
  // the persona's U2 gap: the title must sit on the ROOT (the box is not the checkbox's hover
  // target), pinned on the whole-form pair, which carries no [data-row] anchors
  const whole = await page.evaluate((sel) => {
    const root = [...document.querySelectorAll(`${sel} [data-u2="func-form"] .u2-input-root`)]
      .find((r) => r.querySelector('.u2-input-label')?.textContent === 'Tank Extension');
    return root == null ? {present: false} : {present: true,
      disabled: root.classList.contains('u2-input-disabled'), title: root.getAttribute('title')};
  }, WHOLE);
  ok('func-form-expressions/4b/the-whole-form-disabled-tankextension-root-carries-the-enabled-when-title',
    whole.present && whole.disabled && whole.title === 'Enabled when: tankVolume > 50',
    JSON.stringify(whole));
}

async function checkExpressionValidator(page) {
  await page.locator(`${u2Cell('minAge')} input`).fill('15');
  await page.waitForFunction(({u2, dart}) =>
    (document.querySelector(`${u2} .u2-input-error`)?.textContent ?? '') !== '' &&
    document.querySelector(`${dart} .d4-invalid`) != null,
  {u2: u2Cell('minAge'), dart: dartCell('minAge')}, {timeout: 8000});
  const msg = await errText(page, u2Cell('minAge'));
  await shot(page, 'func-form-expressions-5-invalid');
  await page.locator(`${u2Cell('minAge')} input`).fill('25');
  await page.waitForFunction(({u2, dart}) =>
    (document.querySelector(`${u2} .u2-input-error`)?.textContent ?? '') === '' &&
    document.querySelector(`${dart} .d4-invalid`) == null,
  {u2: u2Cell('minAge'), dart: dartCell('minAge')}, {timeout: 8000});
  ok('func-form-expressions/5a/the-expression-validator-fails-with-the-expression-text-and-clears-on-both-sides',
    msg === 'minAge > 18', `u2 message="${msg}" (the documented false→expression-text parity)`);
}

async function checkNamedValidator(page) {
  // the verdict lands async (evalParamValidators round trip) — poll for it
  await page.locator(`${u2Cell('code')} input`).fill('abc');
  await page.waitForFunction(({u2, dart}) =>
    (document.querySelector(`${u2} .u2-input-error`)?.textContent ?? '') !== '' &&
    document.querySelector(`${dart} .d4-invalid`) != null,
  {u2: u2Cell('code'), dart: dartCell('code')}, {timeout: 8000});
  const msg = await errText(page, u2Cell('code'));
  await page.locator(`${u2Cell('code')} input`).fill('X12');
  await page.waitForFunction(({u2, dart}) =>
    (document.querySelector(`${u2} .u2-input-error`)?.textContent ?? '') === '' &&
    document.querySelector(`${dart} .d4-invalid`) == null,
  {u2: u2Cell('code'), dart: dartCell('code')}, {timeout: 8000});
  ok('func-form-expressions/6a/the-named-validator-verdict-lands-on-u2-and-dart-and-clears-on-a-pass',
    msg === 'Code must start with X', `u2 message="${msg}"`);
}

async function checkUnlistedCategoryRow(page) {
  // the W1 matrix carries its own `notes` param — scope to the W4 grid (the one grid with
  // the vehicle's batteryCapacity row)
  const row = await page.evaluate((page) => {
    const grid = [...document.querySelectorAll(`${page} .u2demo-ab`)]
      .find((g) => g.querySelector('[data-row="batteryCapacity"]') != null);
    return {
      u2Input: grid?.querySelector('.u2-input-root[data-row="notes"] input, ' +
        '.u2-input-root[data-row="notes"] textarea') != null,
      placeholder: [...grid?.querySelectorAll('[data-row="notes"]') ?? []]
        .some((c) => c.textContent.includes('not rendered by the Dart form')),
      dartInput: grid?.querySelector('.ui-input-root[data-row="notes"]') != null,
    };
  }, PAGE);
  await shot(page, 'func-form-expressions-7-notes');
  ok('func-form-expressions/7a/the-unlisted-category-row-renders-on-u2-with-a-placeholder-dart-cell',
    row.u2Input && row.placeholder && !row.dartInput, JSON.stringify(row));
}

async function checkRunGateOnHide(page) {
  const gate = `${PAGE} .u2demo-w4-gate button`;
  const state = (sel) => {
    const b = document.querySelector(sel);
    return {disabled: b?.disabled, title: b?.title};
  };
  const blocked = await page.evaluate(state, gate);
  await flipType(page, 'Electric', u2Cell('cylinders'), false);
  await page.waitForFunction((sel) =>
    document.querySelector(sel)?.disabled === false, gate, {timeout: 8000});
  const unblocked = await page.evaluate(state, gate);
  await shot(page, 'func-form-expressions-8-gate');
  await flipType(page, 'ICE', u2Cell('cylinders'), true);
  await page.waitForFunction((sel) =>
    document.querySelector(sel)?.disabled === true, gate, {timeout: 8000});
  const reblocked = await page.evaluate(state, gate);
  ok('func-form-expressions/8a/hiding-the-required-field-unblocks-run-and-the-tooltip-names-only-visible-fields',
    blocked.disabled === true && blocked.title === 'Fill Cylinders to run' &&
    unblocked.disabled === false && unblocked.title === '' &&
    reblocked.disabled === true && reblocked.title === 'Fill Cylinders to run',
    `blocked=${JSON.stringify(blocked)} unblocked=${JSON.stringify(unblocked)} ` +
    `reblocked=${JSON.stringify(reblocked)}`);
}

// still ICE from /8: every Engine-section input root of the whole-form pair
const setEngineHidden = (page, hidden) => page.evaluate(({sel, hidden}) => {
  const rows = document.querySelector(`${sel} [data-u2="func-form"] .u2-form-rows`);
  const kids = [...(rows?.children ?? [])];
  const start = kids.findIndex((k) => k.classList.contains('u2-form-category') &&
    k.textContent === 'Engine');
  if (start < 0)
    return 0;
  let n = 0;
  for (let i = start + 1; i < kids.length && !kids[i].classList.contains('u2-form-category'); i++)
    if (kids[i].classList.contains('u2-input-root')) {
      kids[i].hidden = hidden;
      n++;
    }
  return n;
}, {sel: WHOLE, hidden});

async function checkDirectHide(page) {
  // a consumer hiding inputs directly (`root.hidden = true`, no expression involved): the form's
  // MutationObserver must still collapse the emptied headers — the dom-shim cannot cover this
  const count = await setEngineHidden(page, true);
  await page.waitForFunction((sel) => {
    const headers = [...document.querySelectorAll(`${sel} [data-u2="func-form"] .u2-form-category`)];
    const engine = headers.find((h) => h.textContent === 'Engine');
    const power = headers.find((h) => h.textContent === 'Power');
    return engine != null && engine.hidden === true &&
      getComputedStyle(engine).display === 'none' && power?.hidden === true;
  }, WHOLE, {timeout: 8000});
  await shot(page, 'func-form-expressions-9-direct-hide');
  ok('func-form-expressions/9a/direct-root-hidden-writes-on-every-engine-field-collapse-the-engine-and-power-headers',
    count === 3, `hid ${count} Engine roots`);
  await setEngineHidden(page, false);
  await page.waitForFunction((sel) => {
    const headers = [...document.querySelectorAll(`${sel} [data-u2="func-form"] .u2-form-category`)];
    const engine = headers.find((h) => h.textContent === 'Engine');
    const power = headers.find((h) => h.textContent === 'Power');
    return engine != null && engine.hidden === false &&
      getComputedStyle(engine).display !== 'none' && power?.hidden === false;
  }, WHOLE, {timeout: 8000});
  ok('func-form-expressions/9b/restoring-the-roots-brings-the-headers-back', true,
    'Engine and Power headers returned');
}

export const checks = [
  {id: 'func-form-expressions/1 the Electric-only field opens hidden on both sides', run: checkInitialState},
  {id: 'func-form-expressions/2 the type flip moves visibility on both columns, one tick', run: checkTypeFlip},
  {id: 'func-form-expressions/3 u2 header auto-hide is visual and reversible; Dart\'s header stays', run: checkHeaderAutoHide},
  {id: 'func-form-expressions/4 the enabled: expression flips tankExtension both ways', run: checkEnabledFlip},
  {id: 'func-form-expressions/5 the expression validator verdict appears and clears', run: checkExpressionValidator},
  {id: 'func-form-expressions/6 the named-validator verdict lands async and clears', run: checkNamedValidator},
  {id: 'func-form-expressions/7 the unlisted-category row is u2-only (#18)', run: checkUnlistedCategoryRow},
  {id: 'func-form-expressions/8 the Run gate flips on expression hide (#20)', run: checkRunGateOnHide},
  {id: 'func-form-expressions/9 direct root.hidden writes collapse and restore the headers', run: checkDirectHide},
];
