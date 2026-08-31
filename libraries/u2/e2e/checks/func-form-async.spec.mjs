/* W2 of the u2 func-call editor: the bench's third section, "Dynamic sources"
   (packages/U2Demo/src/func-convergence.ts) — every source is a function over ONE shared
   FuncCall: `city` re-asks `fceW2Cities(country)` after a country edit, `model` propagates its
   picked row into `mpg`/`cyl`, `site` suggests from the typed text, `cohort` opens with its
   computed default written into the call (R6), `seed`'s default command is broken on purpose,
   and a stale `city` is pruned INTO the call (divergence #9). Section 3 shares the page and the
   `u2demo-ab` grid classes with section 1, so its status lines are the SECOND `sync = ` /
   `inputs = ` lines and every field selector scopes through a section-3-only row name. All
   waits poll DOM state — the dependent re-eval is debounced at 200 ms, never slept for. */
import {ok, openApp, shot} from '../local.mjs';
import {APP} from '../lib.mjs';

/** The bench page, told apart from the inputs-convergence page by the gate row only it has. */
const PAGE = '.u2demo-page:has([data-row="batchId"])';

const u2Cell = (name) => `${PAGE} .u2-input-root[data-row="${name}"]`;
const dartCell = (name) => `${PAGE} .ui-input-root[data-row="${name}"]`;

const statuses = (page) => page.evaluate((sel) =>
  [...document.querySelector(sel)?.querySelectorAll('.u2demo-status') ?? []].map((s) => s.textContent), PAGE);

/** Section 3's status lines are the second of each prefix on the page. */
const dynStatus = async (page, prefix) => (await statuses(page)).filter((t) => t.startsWith(prefix))[1] ?? '';

const counters = async (page) => {
  const m = (await dynStatus(page, 'sync = ')).match(/u2 → Dart (\d+) · Dart → u2 (\d+)/);
  return m === null ? {toDart: -1, toU2: -1} : {toDart: Number(m[1]), toU2: Number(m[2])};
};

const value = (page, sel) => page.evaluate((sel) => document.querySelector(sel)?.value ?? null, sel);

const options = (page, sel) => page.evaluate((sel) =>
  [...document.querySelector(sel)?.options ?? []].map((o) => o.value), sel);

const waitOption = (page, sel, item) => page.waitForFunction(({sel, item}) =>
  [...document.querySelector(sel)?.options ?? []].some((o) => o.value === item),
{sel, item}, {timeout: 10000});

const waitValue = (page, sel, expected) => page.waitForFunction(({sel, expected}) =>
  document.querySelector(sel)?.value === expected, {sel, expected}, {timeout: 8000});

const waitDynInputs = (page, fragment) => page.waitForFunction(({sel, fragment}) =>
  [...document.querySelector(sel)?.querySelectorAll('.u2demo-status') ?? []]
    .filter((s) => s.textContent.startsWith('inputs = '))[1]?.textContent.includes(fragment),
{sel: PAGE, fragment}, {timeout: 10000});

/** One macrotask past the last observed state — what an echo storm would need to show itself. */
const settle = (page) => page.evaluate(() => new Promise((resolve) => setTimeout(resolve, 120)));

export async function fixture(page) {
  await page.evaluate(() => localStorage.setItem('u2demo.tab', 'func-convergence'));
  await openApp(page, APP.package, 'u2DemoApp');
  await page.waitForSelector(`${PAGE} [data-row="regions"]`, {timeout: 30000});
}

async function checkSettledBoot(page) {
  await waitOption(page, `${u2Cell('city')} select`, 'FR-1');
  await waitOption(page, `${dartCell('city')} select`, 'FR-1');
  // the computed default is an async round-trip written into the call at open (R6)
  await waitDynInputs(page, 'cohort: 4');
  const u2City = await options(page, `${u2Cell('city')} select`);
  const dartCity = await options(page, `${dartCell('city')} select`);
  const sync = await dynStatus(page, 'sync = ');
  const inputs = await dynStatus(page, 'inputs = ');
  await shot(page, 'func-form-async-1-rest');
  ok('func-form-async/1a/the-preset-country-populates-both-city-selects-and-the-computed-default-is-in-the-call',
    u2City.includes('FR-1') && u2City.includes('FR-2') && dartCity.includes('FR-1') && dartCity.includes('FR-2') &&
    inputs.includes('cohort: 4') && sync === 'sync = u2 → Dart 0 · Dart → u2 0',
    `u2 city=${JSON.stringify(u2City)} dart city=${JSON.stringify(dartCity)} sync="${sync}" ` +
    `inputs: ${inputs.slice(0, 200)}`);
}

async function checkDependentRefresh(page) {
  const before = await counters(page);
  await page.locator(`${u2Cell('country')} select`).selectOption('DE');
  await waitOption(page, `${u2Cell('city')} select`, 'DE-1');
  await waitOption(page, `${dartCell('city')} select`, 'DE-2');
  const after = await counters(page);
  await settle(page);
  const later = await counters(page);
  const u2City = await options(page, `${u2Cell('city')} select`);
  await shot(page, 'func-form-async-2-country-edit');
  ok('func-form-async/2a/a-u2-country-edit-re-populates-both-city-selects-one-counter-tick-no-echo',
    u2City.includes('DE-1') && u2City.includes('DE-2') &&
    after.toDart === before.toDart + 1 && after.toU2 === before.toU2 &&
    later.toDart === after.toDart && later.toU2 === after.toU2,
    `u2 city=${JSON.stringify(u2City)} counters ${JSON.stringify(before)} → ${JSON.stringify(after)} → ` +
    `settled ${JSON.stringify(later)}`);
}

async function checkSuggestions(page) {
  await page.locator(`${u2Cell('site')} input`).click();
  await page.locator(`${u2Cell('site')} input`).fill('par');
  // the popup is overlay-hosted outside PAGE; nothing else opens a typeahead here
  await page.waitForFunction(() =>
    [...document.querySelectorAll('.u2-typeahead-popup .u2-typeahead-option')]
      .some((el) => el.textContent.includes('par-1')), null, {timeout: 10000});
  const rows = await page.evaluate(() =>
    [...document.querySelectorAll('.u2-typeahead-popup .u2-typeahead-option')].map((el) => el.textContent));
  await page.locator('.u2-typeahead-popup .u2-typeahead-option').first().click();
  await waitValue(page, `${dartCell('site')} input.ui-input-editor`, 'par-1');
  await waitDynInputs(page, "site: 'par-1'");
  ok('func-form-async/3a/the-typed-text-drives-the-popup-and-a-pick-lands-in-the-dart-input-and-the-call',
    rows.length > 0 && rows.every((t) => t.includes('par')) && rows[0].includes('par-1') &&
    await value(page, `${u2Cell('site')} input`) === 'par-1',
    `popup rows=${JSON.stringify(rows)} u2 site="${await value(page, `${u2Cell('site')} input`)}"`);
}

async function checkPropagate(page) {
  await waitOption(page, `${u2Cell('model')} select`, 'Mazda RX4');
  await page.locator(`${u2Cell('model')} select`).selectOption('Mazda RX4');
  await waitValue(page, `${u2Cell('mpg')} input`, '21');
  await waitValue(page, `${dartCell('mpg')} input.ui-input-editor`, '21');
  await waitDynInputs(page, 'mpg: 21');
  const inputs = await dynStatus(page, 'inputs = ');
  ok('func-form-async/4a/a-u2-model-pick-fills-mpg-and-cyl-on-both-sides-and-in-the-call',
    await value(page, `${u2Cell('cyl')} input`) === '6' &&
    await value(page, `${dartCell('cyl')} input.ui-input-editor`) === '6' &&
    inputs.includes('mpg: 21') && inputs.includes('cyl: 6'),
    `u2 cyl="${await value(page, `${u2Cell('cyl')} input`)}" ` +
    `dart cyl="${await value(page, `${dartCell('cyl')} input.ui-input-editor`)}" inputs: ${inputs.slice(0, 250)}`);
}

async function checkBrokenDefault(page) {
  const seed = await page.evaluate((sel) => {
    const err = document.querySelector(`${sel} .u2-param-source-error`);
    return err === null ? null :
      {text: err.textContent, retry: err.querySelector('.u2-param-source-retry')?.textContent ?? null};
  }, u2Cell('seed'));
  ok('func-form-async/5a/the-broken-default-keeps-its-error-on-the-seed-field-with-a-retry',
    seed !== null && seed.text.includes('Couldn\'t compute the default:') &&
    seed.text.includes('Variable "nosuchvar" not found') && seed.retry === 'Retry',
    JSON.stringify(seed).slice(0, 300));
}

async function checkPruneToCall(page) {
  await page.locator(`${u2Cell('city')} select`).selectOption('DE-1');
  await waitDynInputs(page, "city: 'DE-1'");
  await page.locator(`${u2Cell('country')} select`).selectOption('US');
  await waitOption(page, `${u2Cell('city')} select`, 'US-1');
  await waitOption(page, `${dartCell('city')} select`, 'US-1');
  await waitDynInputs(page, 'city: null');
  const u2City = await value(page, `${u2Cell('city')} select`);
  const dartCity = await page.evaluate((sel) => {
    const select = document.querySelector(sel);
    return select === null ? null : {value: select.value, index: select.selectedIndex};
  }, `${dartCell('city')} select`);
  ok('func-form-async/6a/the-stale-city-is-pruned-into-the-call-and-both-selects-go-blank',
    u2City === '' && dartCity !== null && (dartCity.value === '' || dartCity.index === -1),
    `u2 select="${u2City}" dart=${JSON.stringify(dartCity)}`);
}

export const checks = [
  {id: 'func-form-async/1 boot: the dynamic section settles — city loaded both sides, default written', run: checkSettledBoot},
  {id: 'func-form-async/2 dependent refresh: one u2 country edit, both city selects, one tick', run: checkDependentRefresh},
  {id: 'func-form-async/3 suggestions: the typed text drives the popup and a pick lands in the call', run: checkSuggestions},
  {id: 'func-form-async/4 propagate: a model pick fills mpg and cyl on both sides', run: checkPropagate},
  {id: 'func-form-async/5 broken default: the seed row keeps the error and a Retry', run: checkBrokenDefault},
  {id: 'func-form-async/6 prune-to-call: a picked city goes null when country changes', run: checkPruneToCall},
];
