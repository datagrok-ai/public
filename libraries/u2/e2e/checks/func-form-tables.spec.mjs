/* W3 of the u2 func-call editor: the bench's fourth section, "Tables & columns"
   (packages/U2Demo/src/func-convergence.ts) — one shared call over a `DG.Script.create` fixture:
   `df` is a table param both sides auto-fill at open, `num`/`cat` are filtered column params
   auto-picked by name similarity (real writes, in `inputs` before any edit), `metrics` is a
   numerical column_list that is never auto-picked. The call holds OBJECTS: the bench's `fmt`
   prints a real DG.Column as `Column 'x'` and a raw string as `'x'`, so the `inputs` line is the
   object-not-string discriminator too — the page exposes no handle to the call itself. Section 4
   shares the page and the `u2demo-ab` classes with sections 1 and 3, so its status lines are the
   LAST `sync = ` / `inputs = ` pair and every field selector scopes through a section-4-only row
   name. All waits poll DOM state — table retargets and re-defaults are never slept for.

   Deliberately absent: a SAME-frame u2 single-column pick (select `height` while `num` holds
   `age`). Any external same-frame write to a single-column param makes the shipped Dart form's
   ColumnComboBox refresh re-enter its own stream ("Bad state: Cannot fire new event. Controller
   is already firing an event") — reproduced with a raw `DG.InputForm.forFuncCall` +
   `setParamValue('num', column)` and no u2 code at all, so it is a platform defect, not ours.
   The value lands and the combo follows; only the console error is the defect. Until it is
   fixed (or ruled into the NOISE list), the object-write and combo-follow assertions ride the
   cross-frame re-defaults below, which rebuild the Dart combo instead of setting it. */
import {ok, openApp, shot} from '../local.mjs';
import {APP} from '../lib.mjs';

/** The bench page, told apart from the inputs-convergence page by the gate row only it has. */
const PAGE = '.u2demo-page:has([data-row="batchId"])';

const u2Cell = (name) => `${PAGE} .u2-input-root[data-row="${name}"]`;
const dartCell = (name) => `${PAGE} .ui-input-root[data-row="${name}"]`;

/** Section 4's status lines are the last of each prefix on the page. */
const lastStatus = (page, prefix) => page.evaluate(({sel, prefix}) =>
  [...document.querySelector(sel)?.querySelectorAll('.u2demo-status') ?? []]
    .map((s) => s.textContent).filter((t) => t.startsWith(prefix)).pop() ?? '',
{sel: PAGE, prefix});

const counters = async (page) => {
  const m = (await lastStatus(page, 'sync = ')).match(/u2 → Dart (\d+) · Dart → u2 (\d+)/);
  return m === null ? {toDart: -1, toU2: -1} : {toDart: Number(m[1]), toU2: Number(m[2])};
};

const waitW3Inputs = (page, fragment) => page.waitForFunction(({sel, fragment}) => {
  const lines = [...document.querySelector(sel)?.querySelectorAll('.u2demo-status') ?? []]
    .map((s) => s.textContent).filter((t) => t.startsWith('inputs = '));
  return (lines[lines.length - 1] ?? '').includes(fragment);
}, {sel: PAGE, fragment}, {timeout: 10000});

const u2Select = (page, name) => page.evaluate((sel) =>
  document.querySelector(`${sel} select`)?.value ?? null, u2Cell(name));
const u2Options = (page, name) => page.evaluate((sel) =>
  [...document.querySelector(`${sel} select`)?.options ?? []].map((o) => o.value), u2Cell(name));
const dartText = (page, name) => page.evaluate((sel) =>
  document.querySelector(sel)?.textContent ?? '', dartCell(name));
const dartTableText = (page, name) => page.evaluate((sel) => {
  const select = document.querySelector(`${sel} select`);
  return select == null ? null : (select.selectedOptions[0]?.textContent ?? select.value);
}, dartCell(name));

/** One debounce window past the last observed state — what an echo storm would need to show
 * itself (the state itself is always waited for above, never slept for). */
const settle = (page) => page.evaluate(() => new Promise((resolve) => setTimeout(resolve, 250)));

export async function fixture(page) {
  await page.evaluate(() => localStorage.setItem('u2demo.tab', 'func-convergence'));
  await openApp(page, APP.package, 'u2DemoApp');
  // the section-4 rows are appended after both forms of every section built
  await page.waitForSelector(`${PAGE} [data-row="metrics"]`, {timeout: 30000});
}

async function checkBoot(page) {
  await waitW3Inputs(page, "df: DataFrame '");
  const u2Df = await u2Select(page, 'df');
  const dartDf = await dartTableText(page, 'df');
  const u2Items = await u2Options(page, 'df');
  const dartItems = await page.evaluate((sel) =>
    [...document.querySelector(`${sel} select`)?.options ?? []].map((o) => o.textContent), dartCell('df'));
  const inputs = await lastStatus(page, 'inputs = ');
  await shot(page, 'func-form-tables-1-rest');
  ok('func-form-tables/1a/both-table-selects-offer-the-fixtures-and-auto-filled-the-same-table-into-the-call',
    u2Df != null && u2Df !== '' && dartDf === u2Df &&
    ['fceW3Demog', 'fceW3Alt'].every((n) => u2Items.includes(n) && dartItems.includes(n)) &&
    inputs.includes(`df: DataFrame '${u2Df}'`),
    `u2=${JSON.stringify(u2Df)} dart=${JSON.stringify(dartDf)} u2 items=${JSON.stringify(u2Items)} ` +
    `inputs: ${inputs.slice(0, 200)}`);
  // `fmt` prints objects as `Column 'x'` and strings as `'x'` — this is the objects-in-the-call
  // pin; column lists are never auto-picked, so the call honestly holds null at rest
  ok('func-form-tables/1b/num-and-cat-hold-column-objects-in-the-call-metrics-never-auto-picked',
    /num: Column '\w+'/.test(inputs) && /cat: Column '\w+'/.test(inputs) &&
    inputs.includes('metrics: null'),
    inputs.slice(0, 260));
}

async function checkTableRetarget(page) {
  if (await u2Select(page, 'df') === 'fceW3Alt') {
    await page.locator(`${u2Cell('df')} select`).selectOption('fceW3Demog');
    await waitW3Inputs(page, "df: DataFrame 'fceW3Demog'");
    await waitW3Inputs(page, "cat: Column 'sex'");
    await settle(page);
  }
  const before = await counters(page);
  await page.locator(`${u2Cell('df')} select`).selectOption('fceW3Alt');
  await waitW3Inputs(page, "df: DataFrame 'fceW3Alt'");
  // Alt's one numerical is `weight`, one categorical `label` — algorithm-proof re-defaults;
  // `Column '…'` in the line is the object-not-string pin for a u2-side column write
  await waitW3Inputs(page, "num: Column 'weight'");
  await waitW3Inputs(page, "cat: Column 'label'");
  await settle(page);
  const after = await counters(page);
  await settle(page);
  const later = await counters(page);
  const dartNum = await dartText(page, 'num');
  const dartCat = await dartText(page, 'cat');
  await shot(page, 'func-form-tables-2-alt');
  ok('func-form-tables/2a/a-u2-table-pick-re-defaults-both-columns-as-objects-into-the-call-and-the-dart-combos-follow',
    dartNum.includes('weight') && dartCat.includes('label'),
    `dart num cell="${dartNum.slice(0, 60)}" cat cell="${dartCat.slice(0, 60)}"`);
  ok('func-form-tables/2b/counters-flat-a-settle-later-no-echo-storm',
    later.toDart === after.toDart && later.toU2 === after.toU2,
    `${JSON.stringify(before)} → ${JSON.stringify(after)} → settled ${JSON.stringify(later)} ` +
    '(one pick: df + re-defaulted num/cat + cleared metrics all write through u2 inputs)');
}

async function checkNoMatch(page) {
  await page.evaluate(() => {
    if (!grok.shell.tables.some((t) => t.name === 'fceW3Nums')) {
      const df = DG.DataFrame.fromCsv('x,y\n1,2\n3,4');
      df.name = 'fceW3Nums';
      grok.shell.addTable(df);
    }
  });
  await page.waitForFunction((sel) =>
    [...document.querySelector(`${sel} select`)?.options ?? []].some((o) => o.value === 'fceW3Nums'),
  u2Cell('df'), {timeout: 8000});
  await page.locator(`${u2Cell('df')} select`).selectOption('fceW3Nums');
  await waitW3Inputs(page, "df: DataFrame 'fceW3Nums'");
  await waitW3Inputs(page, 'cat: null');
  const inputs = await lastStatus(page, 'inputs = ');
  const cat = await page.evaluate((sel) => ({
    combo: document.querySelector(`${sel} .u2-columns`) != null,
    summary: document.querySelector(`${sel} .u2-columns-summary`)?.textContent ?? null,
  }), u2Cell('cat'));
  await shot(page, 'func-form-tables-3-nomatch');
  ok('func-form-tables/3a/no-categorical-column-cat-re-defaults-to-null-in-the-call-and-the-u2-combo-empties',
    inputs.includes('cat: null') && /num: Column '(x|y)'/.test(inputs) &&
    cat.combo && cat.summary === '',
    `inputs: ${inputs.slice(0, 200)} cat=${JSON.stringify(cat)}`);
}

/** The multi-column popup is the platform ColumnsGrid in checkbox mode (dialog semantics): its own
 * All/None labels drive the checks DOM-wise — the checkboxes themselves are canvas cells. */
async function checkColumnsPick(page) {
  await page.evaluate((sel) =>
    document.querySelector(sel)?.scrollIntoView({block: 'center'}), u2Cell('metrics'));
  await page.locator(`${u2Cell('metrics')} .u2-columns`).click();
  await page.waitForSelector('.u2-columns-popup .d4-column-grid', {timeout: 5000});
  await page.locator('.u2-columns-popup label', {hasText: /^All$/}).click();
  await page.locator('.u2-columns-popup .u2-columns-buttons button', {hasText: 'OK'}).first().click();
  await waitW3Inputs(page, 'metrics: [x, y]');
  const inputs = await lastStatus(page, 'inputs = ');
  const u2Metrics = await page.evaluate((sel) =>
    document.querySelector(`${sel} .u2-columns-summary`)?.textContent ?? '', u2Cell('metrics'));
  const dartMetrics = await dartText(page, 'metrics');
  await shot(page, 'func-form-tables-4-columns');
  ok('func-form-tables/4a/two-columns-checked-through-the-platform-grid-read-back-as-their-names-and-both-sides-show-2',
    inputs.includes('metrics: [x, y]') && u2Metrics.includes('(2)') && dartMetrics.includes('(2)'),
    `inputs: ${inputs.slice(0, 200)} u2 summary="${u2Metrics}" dart cell="${dartMetrics.slice(0, 60)}"`);

  await page.locator(`${u2Cell('metrics')} .u2-columns`).click();
  await page.waitForSelector('.u2-columns-popup .d4-column-grid', {timeout: 5000});
  await page.locator('.u2-columns-popup label', {hasText: /^None$/}).click();
  await page.locator('.u2-columns-popup .u2-columns-buttons button', {hasText: 'CANCEL'}).first().click();
  await settle(page);
  const after = await lastStatus(page, 'inputs = ');
  const closed = await page.evaluate(() => document.querySelector('.u2-columns-popup') == null);
  ok('func-form-tables/4b/cancel-after-unchecking-everything-discards-the-toggles-zero-writes',
    closed && after.includes('metrics: [x, y]'), after.slice(0, 200));
}

async function checkCloseTable(page) {
  await page.evaluate(() =>
    grok.shell.closeTable(grok.shell.tables.find((t) => t.name === 'fceW3Nums')));
  await waitW3Inputs(page, 'df: null');
  await waitW3Inputs(page, 'num: null');
  await waitW3Inputs(page, 'metrics: []');
  await settle(page);
  const after = await counters(page);
  await settle(page);
  const later = await counters(page);
  const u2Df = await u2Select(page, 'df');
  const dartDf = await dartTableText(page, 'df');
  await shot(page, 'func-form-tables-5-closed');
  ok('func-form-tables/5a/the-closed-table-is-pruned-to-null-in-the-call-with-every-dependent-and-the-u2-select-clears',
    (u2Df === '' || u2Df === null) && later.toDart === after.toDart && later.toU2 === after.toU2,
    `u2 df=${JSON.stringify(u2Df)} dart df select displays=${JSON.stringify(dartDf)} ` +
    `(divergence #11: the Dart select does not display what the call held) counters ` +
    `${JSON.stringify(after)} → settled ${JSON.stringify(later)}`);
}

export const checks = [
  {id: 'func-form-tables/1 boot: fixtures offered, one table auto-filled both sides, columns auto-picked', run: checkBoot},
  {id: 'func-form-tables/2 a u2 table pick re-targets and re-defaults both dependent columns', run: checkTableRetarget},
  {id: 'func-form-tables/3 no-match re-default: cat goes null in the call, the u2 select empties', run: checkNoMatch},
  {id: 'func-form-tables/4 a two-column popup pick lands as a ColumnList, (2) on both sides', run: checkColumnsPick},
  {id: 'func-form-tables/5 closing the value table prunes it and every dependent into the call', run: checkCloseTable},
];
