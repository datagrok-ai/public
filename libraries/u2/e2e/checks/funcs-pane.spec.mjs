/* The U2Demo "Forms" pane (packages/U2Demo/src/funcs.ts): the platform's own searchable
   DG.FunctionsWidget on the left drives an old-vs-new A/B on the right — the Dart
   DG.InputForm.forFuncCall and the u2 funcForm over ONE shared call, so edits sync live.
   A showcase function (fceShowcase) is preselected at boot. Selection goes through the REAL
   Dart widget rows (the same click a user makes), so this also covers the FunctionsWidget
   interop event round-trip. */
import {ok, openApp, shot} from '../local.mjs';
import {APP} from '../lib.mjs';

const PAGE = '.u2demo-funcs';
const LIST = `${PAGE} .u2demo-funcs-list`;
const DETAIL = `${PAGE} .u2demo-funcs-detail`;
const AB = `${DETAIL} .u2demo-funcs-ab`;
const DART_COL = `${AB} > :first-child`;
const U2_COL = `${AB} > :last-child`;

const clickRow = (page, name) => page.evaluate(({sel, name}) => {
  const row = [...document.querySelectorAll(`${sel} .grok-actions-browser-table span`)]
    .find((s) => s.textContent.trim() === name && s.offsetParent != null);
  if (row == null)
    return false;
  row.click();
  return true;
}, {sel: LIST, name});

const detailText = (page) => page.evaluate((sel) =>
  document.querySelector(sel)?.innerText ?? '', DETAIL);

const waitDetail = (page, fragment) => page.waitForFunction(({sel, fragment}) =>
  (document.querySelector(sel)?.innerText ?? '').includes(fragment), {sel: DETAIL, fragment}, {timeout: 8000});

const visibleNames = (page) => page.evaluate((sel) =>
  [...document.querySelectorAll(`${sel} .grok-actions-browser-table span`)]
    .filter((s) => s.offsetParent != null && s.textContent.trim() !== '')
    .map((s) => s.textContent.trim()), LIST);

/** Both columns of the current selection have settled: the Dart form landed (or failed inline)
 * and the u2 form rendered. */
const waitColumns = (page) => page.waitForFunction(({dart, u2}) =>
  document.querySelector(`${u2} .u2-input-root`) != null &&
  (document.querySelector(`${dart} .ui-input-root`) != null ||
   document.querySelector(`${dart} .u2demo-error`) != null),
{dart: DART_COL, u2: U2_COL}, {timeout: 15000});

export async function fixture(page) {
  await page.evaluate(() => localStorage.setItem('u2demo.tab', 'funcs'));
  await openApp(page, APP.package, 'u2DemoApp');
  await page.waitForSelector(`${LIST} .grok-actions-browser-table span`, {timeout: 30000});
}

async function checkBoot(page) {
  await waitColumns(page);
  // the computed default `2 + 2` is evaluated async and written into the call (R6)
  await page.waitForFunction((u2) =>
    [...document.querySelectorAll(`${u2} .u2-input-root input`)].some((i) => i.value === '4'),
  U2_COL, {timeout: 8000});
  const shape = await page.evaluate(({list, detail, dart, u2}) => {
    const el = document.querySelector(list);
    return {
      rows: document.querySelectorAll(`${list} .grok-actions-browser-table span`).length,
      search: document.querySelector(`${list} input`) !== null,
      listBounded: el.clientHeight > 200 && el.clientHeight <= window.innerHeight,
      listNarrow: el.getBoundingClientRect().width <= 260,
      text: document.querySelector(detail)?.innerText ?? '',
      dartInputs: document.querySelectorAll(`${dart} .ui-input-root`).length,
      u2Inputs: document.querySelectorAll(`${u2} .u2-input-root`).length,
    };
  }, {list: LIST, detail: DETAIL, dart: DART_COL, u2: U2_COL});
  await shot(page, 'funcs-pane-1-boot');
  ok('funcs-pane/1a/boot-preselects-the-showcase-with-both-form-columns-and-the-computed-default',
    shape.rows > 100 && shape.search && shape.listBounded && shape.listNarrow &&
    shape.text.includes('fceShowcase') && shape.dartInputs >= 10 && shape.u2Inputs >= 10,
    `rows=${shape.rows} search=${shape.search} bounded=${shape.listBounded} ` +
    `narrow=${shape.listNarrow} dart=${shape.dartInputs} u2=${shape.u2Inputs}`);
}

async function checkSync(page) {
  await page.locator(`${U2_COL} .u2-input-root:has(.u2-input-label:text-is("Name")) input`)
    .fill('Showcase A');
  await page.waitForFunction((dart) => {
    const root = [...document.querySelectorAll(`${dart} .ui-input-root`)]
      .find((r) => r.querySelector('.ui-input-label')?.textContent.trim() === 'Name');
    return root?.querySelector('input')?.value === 'Showcase A';
  }, DART_COL, {timeout: 8000});
  await shot(page, 'funcs-pane-2-synced');
  ok('funcs-pane/2a/a-u2-side-edit-reaches-the-dart-side-through-the-shared-call', true,
    'u2 Name → "Showcase A" → Dart Name follows');
}

async function checkSearchAndSelect(page) {
  await page.locator(`${LIST} input`).first().fill('rand bet');
  await page.waitForFunction((sel) => {
    const names = [...document.querySelectorAll(`${sel} .grok-actions-browser-table span`)]
      .filter((s) => s.offsetParent != null);
    return names.some((s) => s.textContent.trim() === 'Rand Between') &&
      !names.some((s) => s.textContent.trim() === 'Abs');
  }, LIST, {timeout: 8000});
  const names = await visibleNames(page);
  const clicked = await clickRow(page, 'Rand Between');
  await waitDetail(page, 'core:RandBetween(int a, int b): int');
  await waitColumns(page);
  const counts = await page.evaluate(({dart, u2}) => ({
    dart: document.querySelectorAll(`${dart} .ui-input-root input`).length,
    u2: document.querySelectorAll(`${u2} .u2-input-root input`).length,
  }), {dart: DART_COL, u2: U2_COL});
  const text = await detailText(page);
  await shot(page, 'funcs-pane-3-selected');
  ok('funcs-pane/3a/typing-filters-the-list-and-a-click-renders-both-forms',
    clicked && counts.dart === 2 && counts.u2 === 2 && text.includes('Rand Between') &&
    text.includes('random integer') && text.includes('result ='),
    `filtered=${JSON.stringify(names.slice(0, 4))} clicked=${clicked} ` +
    `dart=${counts.dart} u2=${counts.u2}`);
}

async function checkRun(page) {
  const gated = await page.evaluate((sel) => {
    const run = [...document.querySelectorAll(`${sel} button`)]
      .find((b) => b.textContent.trim().toLowerCase() === 'run');
    return run?.disabled === true && run.title !== '';
  }, DETAIL);
  await page.locator(`${U2_COL} .u2-input-root input`).first().fill('5');
  await page.locator(`${U2_COL} .u2-input-root input`).nth(1).fill('6');
  await page.evaluate((sel) => {
    [...document.querySelectorAll(`${sel} button`)]
      .find((b) => b.textContent.trim().toLowerCase() === 'run').click();
  }, DETAIL);
  await waitDetail(page, 'result =\n5');
  ok('funcs-pane/4a/run-is-gated-while-the-required-fields-are-empty-then-executes-with-the-form-values',
    gated && (await detailText(page)).includes('result =\n5'),
    `gated=${gated} ${(await detailText(page)).slice(-60)}`);
}

async function checkRebindAndTableParams(page) {
  await page.locator(`${LIST} input`).first().fill('join tables');
  await page.waitForFunction((sel) =>
    [...document.querySelectorAll(`${sel} .grok-actions-browser-table span`)]
      .some((s) => s.offsetParent != null && s.textContent.trim() === 'Join Tables'), LIST, {timeout: 8000});
  await clickRow(page, 'Join Tables');
  await waitDetail(page, 'core:JoinTables');
  await waitColumns(page);
  const dartState = await page.evaluate((dart) => {
    const col = document.querySelector(dart);
    return col?.querySelector('.ui-input-root') != null ? 'form' :
      col?.querySelector('.u2demo-error') != null ? 'error' : 'none';
  }, DART_COL);
  const text = await detailText(page);
  // dataframe/column_list params route to real fields since W3 — nothing left unsupported here
  ok('funcs-pane/5a/reselection-rebinds-both-columns-and-the-w3-table-params-render-as-fields',
    text.includes('Join Tables') && !text.includes('Rand Between') &&
    !text.includes('not supported by this form yet') && text.includes('Table1') &&
    text.includes('Join Type') && dartState !== 'none',
    `dart=${dartState} ${text.slice(0, 160)}`);
}

export const checks = [
  {id: 'funcs-pane/1 boot: the Forms tab preselects the showcase, both columns render', run: checkBoot},
  {id: 'funcs-pane/2 an edit on the u2 side syncs to the Dart side via the shared call', run: checkSync},
  {id: 'funcs-pane/3 search filters, a click renders both forms', run: checkSearchAndSelect},
  {id: 'funcs-pane/4 Run executes the prepared call with the form values', run: checkRun},
  {id: 'funcs-pane/5 reselection rebinds; Join Tables renders fully in both columns', run: checkRebindAndTableParams},
];
