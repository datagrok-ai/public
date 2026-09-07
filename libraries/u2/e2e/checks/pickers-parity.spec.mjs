/* The pickers-parity wave (plan-pickers.md WO-6): the u2 table field's platform icon rail and
   the column combo's REAL Dart ColumnsGrid dropdown, on the bench's section 4 — permanent
   coverage distilled from the WO-5 live driver. Esc here pins the WO-5 live finding: the Dart
   search box swallows a bubble-phase Esc even when empty (text_input.dart:148-151), so the combo
   owns a capture-phase listener that closes first.

   Console-noise discipline (probed live, wo6/probe-d10.mjs): a VALUE-MOVING single-column pick
   or closed-combo arrow fires the platform's defect-#10 'Cannot fire new event' error once
   (same-frame write into the paired Dart form); same-value confirming picks, Esc dismissals,
   table retargets and the empty-row null pick are clean. The NOISE list must not grow, so every
   pick this file drives is either same-value (the keyboard and pointer pins — value movement is
   already covered call-side by func-form-tables/2's re-defaults) or the null pick; and nothing
   here touches the Forms tab's Records select (the standing lead-pending NOISE call). */
import {ok, openApp, shot} from '../local.mjs';
import {APP} from '../lib.mjs';

const PAGE = '.u2demo-page:has([data-row="batchId"])';
const HOST = '.u2-column-combo-host';

const u2Cell = (name) => `${PAGE} .u2-input-root[data-row="${name}"]`;

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

const settle = (page) => page.evaluate(() => new Promise((resolve) => setTimeout(resolve, 250)));

const summary = (page, name) => page.evaluate((sel) =>
  document.querySelector(`${sel} .u2-columns-summary`)?.textContent ?? null, u2Cell(name));

const popupOpen = (page) => page.evaluate((h) => document.querySelector(h) != null, HOST);

const waitPopupClosed = (page) =>
  page.waitForFunction((h) => document.querySelector(h) == null, HOST, {timeout: 8000});

const openCombo = async (page, name) => {
  await page.evaluate((sel) =>
    document.querySelector(sel)?.scrollIntoView({block: 'center'}), u2Cell(name));
  await page.locator(`${u2Cell(name)} .u2-columns`).click();
  await page.waitForSelector(`${HOST} .d4-column-grid`, {timeout: 8000});
};

export async function fixture(page) {
  await page.evaluate(() => localStorage.setItem('u2demo.tab', 'func-convergence'));
  await openApp(page, APP.package, 'u2DemoApp');
  await page.waitForSelector(`${PAGE} [data-row="metrics"]`, {timeout: 30000});
  // earlier files leave other tables open (or none selected): pin the fixture frame, whose
  // re-default is the stable 'age' — a cross-frame retarget, console-clean
  await waitW3Inputs(page, "df: DataFrame '");
  if (await page.evaluate((sel) =>
    document.querySelector(`${sel} select`)?.value, u2Cell('df')) !== 'fceW3Demog') {
    await page.locator(`${u2Cell('df')} select`).selectOption('fceW3Demog');
    await waitW3Inputs(page, "df: DataFrame 'fceW3Demog'");
  }
  await page.waitForFunction((sel) =>
    document.querySelector(`${sel} .u2-columns-summary`)?.textContent === 'age',
  u2Cell('num'), {timeout: 10000});
  await settle(page);
}

async function checkIconRail(page) {
  await page.evaluate((sel) =>
    document.querySelector(sel)?.scrollIntoView({block: 'center'}), u2Cell('df'));
  const icons = await page.evaluate((sel) =>
    [...document.querySelectorAll(`${sel} .u2-input-options .u2-icon-btn i`)]
      .map((i) => [...i.classList].find((c) => c.startsWith('fa-') && !c.startsWith('fa-l'))), u2Cell('df'));
  await shot(page, 'pickers-parity-1-icons');
  ok('pickers-parity/1a/the-u2-table-field-carries-the-3-platform-icons-in-order',
    icons.length === 3 && icons[0] === 'fa-folder-open' && icons[1] === 'fa-folder-tree' &&
    icons[2] === 'fa-database', JSON.stringify(icons));
}

async function checkQueryDialog(page) {
  const dfBefore = await page.evaluate((sel) =>
    document.querySelector(`${sel} select`)?.value, u2Cell('df'));
  const inputsBefore = await lastStatus(page, 'inputs = ');
  const before = await counters(page);
  await page.locator(`${u2Cell('df')} .u2-icon-btn:has(i.fa-database)`).click();
  await page.waitForFunction(() => DG.Dialog.getOpenDialogs().length > 0, null, {timeout: 8000});
  const titles = await page.evaluate(() => DG.Dialog.getOpenDialogs().map((d) => d.title ?? ''));
  await shot(page, 'pickers-parity-2-query-dialog');
  await page.keyboard.press('Escape');
  await page.waitForFunction(() => DG.Dialog.getOpenDialogs().length === 0, null, {timeout: 8000});
  await settle(page);
  const after = await counters(page);
  ok('pickers-parity/2a/the-database-icon-opens-the-real-dialog-and-esc-leaves-the-call-untouched',
    titles.some((t) => t.includes('Select a database query')) &&
    await page.evaluate((sel) =>
      document.querySelector(`${sel} select`)?.value, u2Cell('df')) === dfBefore &&
    await lastStatus(page, 'inputs = ') === inputsBefore &&
    after.toDart === before.toDart && after.toU2 === before.toU2,
    `titles=${JSON.stringify(titles)} counters ${JSON.stringify(before)}→${JSON.stringify(after)}`);
}

async function checkComboOpenAndEsc(page) {
  const numBefore = await summary(page, 'num');
  const inputsBefore = await lastStatus(page, 'inputs = ');
  const before = await counters(page);
  await openCombo(page, 'num');
  const state = await page.evaluate(({h, cell}) => {
    const host = document.querySelector(h);
    const active = document.activeElement;
    return {
      inOverlay: host?.closest('.u2-overlay') != null,
      grid: host?.querySelector('.d4-column-grid') != null,
      searchFocused: active?.tagName === 'INPUT' && host?.contains(active) === true,
      expanded: document.querySelector(`${cell} .u2-columns`)?.getAttribute('aria-expanded'),
    };
  }, {h: HOST, cell: u2Cell('num')});
  await shot(page, 'pickers-parity-3-grid-open');
  ok('pickers-parity/3a/the-combo-opens-the-real-dart-columnsgrid-in-the-u2-overlay-search-focused',
    state.inOverlay && state.grid && state.searchFocused && state.expanded === 'true',
    JSON.stringify(state));

  // Esc lands in the Dart search box, which swallows it at bubble even when empty — the
  // capture-phase close (the WO-5 live finding's fix) is what this pins
  await page.keyboard.press('Escape');
  await waitPopupClosed(page);
  await settle(page);
  const after = await counters(page);
  const expanded = await page.evaluate((cell) =>
    document.querySelector(`${cell} .u2-columns`)?.getAttribute('aria-expanded'), u2Cell('num'));
  ok('pickers-parity/3b/esc-from-the-platform-search-closes-the-popup-without-a-write',
    await popupOpen(page) === false && expanded === 'false' &&
    await summary(page, 'num') === numBefore &&
    await lastStatus(page, 'inputs = ') === inputsBefore &&
    after.toDart === before.toDart && after.toU2 === before.toU2,
    `num=${numBefore} counters ${JSON.stringify(before)}→${JSON.stringify(after)}`);
}

async function checkKeyboardPick(page) {
  // a CONFIRMING same-value pick: the keyboard path end-to-end with zero call writes (the
  // console-noise discipline above; value movement is func-form-tables/2's re-default pin)
  const current = await summary(page, 'num');
  const before = await counters(page);
  await openCombo(page, 'num');
  await page.keyboard.type(current);
  await settle(page);
  await shot(page, 'pickers-parity-4-search-filtered');
  await page.keyboard.press('Enter');
  await waitPopupClosed(page);
  await settle(page);
  const after = await counters(page);
  ok('pickers-parity/4a/type-to-filter-then-enter-picks-the-unique-match-and-closes',
    await summary(page, 'num') === current &&
    after.toDart === before.toDart && after.toU2 === before.toU2,
    `confirming pick of '${current}': counters ${JSON.stringify(before)}→${JSON.stringify(after)}`);
}

async function checkPointerPick(page) {
  // the one real mouse row-click pin: header 20 + rowHeight 16 (column_grid.dart:384), row 0
  // is the current 'age' — a same-value pick, closing without a write
  const current = await summary(page, 'num');
  await openCombo(page, 'num');
  const canvas = await page.evaluate((h) =>
    [...document.querySelectorAll(`${h} .d4-column-grid canvas`)]
      .map((el) => el.getBoundingClientRect())
      .filter((r) => r.width > 0 && r.height > 0)
      .map((r) => ({left: r.left, top: r.top}))[0] ?? null, HOST);
  await page.mouse.click(canvas.left + 40, canvas.top + 20 + 8);
  await waitPopupClosed(page);
  await settle(page);
  ok('pickers-parity/5a/a-real-mouse-click-on-a-grid-row-picks-and-closes',
    current === 'age' && await summary(page, 'num') === 'age',
    `row 0 clicked at canvas(+40,+28), num=${await summary(page, 'num')}`);
}

async function checkPopupClosesOnRetarget(page) {
  await openCombo(page, 'num');
  await page.locator(`${u2Cell('df')} select`).selectOption('fceW3Alt');
  await waitPopupClosed(page);
  await waitW3Inputs(page, "num: Column 'weight'");
  await settle(page);
  const orphans = await page.evaluate((h) => document.querySelectorAll(h).length, HOST);
  await shot(page, 'pickers-parity-6-retarget');
  ok('pickers-parity/6a/a-table-switch-closes-the-open-popup-and-leaves-no-orphaned-overlay',
    orphans === 0 && await summary(page, 'num') === 'weight',
    `orphan hosts=${orphans} num=${await summary(page, 'num')}`);
  await page.locator(`${u2Cell('df')} select`).selectOption('fceW3Demog');
  await waitW3Inputs(page, "num: Column 'age'");
  await settle(page);
}

async function checkEmptyRowNull(page) {
  // nullable lives on the Forms tab's Age Col (the bench's num/cat are non-nullable by
  // annotation, so their grids carry no empty row); the null pick is console-clean (probed)
  await page.evaluate(() => localStorage.setItem('u2demo.tab', 'funcs'));
  await openApp(page, APP.package, 'u2DemoApp');
  await page.waitForSelector('.u2demo-funcs-ab', {timeout: 30000});
  const ageRoot = () => {
    const form = document.querySelector('.u2demo-funcs-col [data-u2="func-form"]');
    return [...form?.querySelectorAll('.u2-input-root') ?? []]
      .find((r) => r.querySelector('label')?.textContent === 'Age Col') ?? null;
  };
  await page.waitForFunction(`((${ageRoot})()?.querySelector('.u2-columns-summary')?.textContent ?? '') !== ''`,
    null, {timeout: 30000});
  const wasPicked = await page.evaluate(
    `(${ageRoot})()?.querySelector('.u2-columns-summary')?.textContent ?? null`);
  await page.evaluate(`(${ageRoot})()?.querySelector('.u2-columns')?.scrollIntoView({block: 'center'})`);
  await page.evaluate(`(${ageRoot})()?.querySelector('.u2-columns')?.click()`);
  await page.waitForSelector(`${HOST} .d4-column-grid`, {timeout: 8000});
  await page.keyboard.press('ArrowDown'); // the first visible row is the addEmpty row
  await page.keyboard.press('Enter');
  await waitPopupClosed(page);
  await settle(page);
  const now = await page.evaluate(
    `(${ageRoot})()?.querySelector('.u2-columns-summary')?.textContent ?? null`);
  await shot(page, 'pickers-parity-7-empty-row');
  ok('pickers-parity/7a/the-empty-row-picks-null-through-the-real-grid',
    wasPicked !== '' && wasPicked !== null && now === '',
    `ageCol was=${JSON.stringify(wasPicked)} now=${JSON.stringify(now)} (shared call: null landed)`);
}

export const checks = [
  {id: 'pickers-parity/1 the table field carries the platform icon rail', run: checkIconRail},
  {id: 'pickers-parity/2 the database icon opens the real dialog, Esc cancels clean', run: checkQueryDialog},
  {id: 'pickers-parity/3 the combo drops the real ColumnsGrid; Esc closes it without a write', run: checkComboOpenAndEsc},
  {id: 'pickers-parity/4 keyboard: type-to-filter, Enter picks, popup closes', run: checkKeyboardPick},
  {id: 'pickers-parity/5 pointer: a real row click picks', run: checkPointerPick},
  {id: 'pickers-parity/6 a table switch closes the popup — no orphaned overlay', run: checkPopupClosesOnRetarget},
  {id: 'pickers-parity/7 the nullable empty row writes null', run: checkEmptyRowNull},
];
