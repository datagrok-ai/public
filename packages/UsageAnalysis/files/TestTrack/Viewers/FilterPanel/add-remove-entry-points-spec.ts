/* ---
realizes: [filters.cp.add-remove-entry-points]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {cardCount, driveOpenMenuLeaf, expectHeaderCounter, expectHeaderCounterQuiet,
  expectHeaderCounterQuietNow, trueCount} from '../../helpers/filter-panel';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const fullRowCount = 5850;

const SELECT_COLUMNS_DIALOG = '.d4-dialog[name="dialog-Select-columns..."]';
const MULTI_VALUE_DIALOG = '.d4-dialog[name="dialog-Select-column-and-separator-for-multi-value-filter"]';
const FILTER_TO_COLUMN_DIALOG = '.d4-dialog[name="dialog-Filter-to-Column"]';
const MULTI_VALUE_FILTER_TYPE = 'multi-value';
const NO_CATEGORICAL_BALLOON = 'No categorical columns found in the table';

async function orderedCaptions(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name')]
      .map((e) => (e.textContent ?? '').trim()));
}

async function holdTrueCount(page: Page, expected: number, why: string, ms = 3000): Promise<void> {
  const deadline = Date.now() + ms;
  for (;;) {
    expect(await trueCount(page), why).toBe(expected);
    if (Date.now() > deadline) return;
    await page.waitForTimeout(500);
  }
}

async function gridHeaderPoint(page: Page, column: string): Promise<{x: number; y: number}> {
  const p = await page.evaluate((col) => {
    const grid = grok.shell.tv.grid;
    const mainGrid = [...document.querySelectorAll('[name="viewer-Grid"]')]
      .find((g) => !g.closest('.d4-filter'));
    const overlay = mainGrid!.querySelector('[name="overlay"]') as HTMLElement;
    const or = overlay.getBoundingClientRect();
    const gc = grid.columns.byName(col);
    return {x: Math.round(or.left + (gc.left + gc.right) / 2), y: Math.round(or.top + grid.colHeaderHeight / 2)};
  }, column);
  return p;
}

async function dragColumnHeaderToPanel(page: Page, column: string): Promise<void> {
  const src = await gridHeaderPoint(page, column);
  const panel = await page.evaluate(() => {
    const fp = document.querySelector('[name="viewer-Filters"]')!;
    const r = fp.getBoundingClientRect();
    return {cx: Math.round(r.x + r.width / 2), cy: Math.round(r.y + r.height / 2)};
  });
  await page.mouse.move(src.x, src.y);
  await page.mouse.down();
  for (let i = 1; i <= 8; i++) {
    const x = Math.round(src.x + (panel.cx - src.x) * i / 8);
    const y = Math.round(src.y + (panel.cy - src.y) * i / 8);
    await page.mouse.move(x, y, {steps: 3});
    await page.waitForTimeout(30);
  }
  try {
    await page.waitForFunction(() => [...document.querySelectorAll('.d4-drop-zone')]
      .some((z) => z.parentElement === document.body && (z.textContent ?? '').trim() === 'Add filter'),
    null, {timeout: 5000, polling: 100});
  }
  catch {
    await page.mouse.up();
    throw new Error(`dragColumnHeaderToPanel(${column}): the "Add filter" drop zone never appeared within 5s — the panel did not start accepting the drag`);
  }
  const zone = await page.evaluate(() => {
    const dz = [...document.querySelectorAll('.d4-drop-zone')]
      .find((z) => z.parentElement === document.body && (z.textContent ?? '').trim() === 'Add filter')!;
    const r = dz.getBoundingClientRect();
    return {x: Math.round(r.x + r.width / 2), y: Math.round(r.y + r.height / 2)};
  });
  await page.mouse.move(zone.x, zone.y, {steps: 4});
  await page.waitForTimeout(120);
  await page.mouse.up();
  await page.waitForTimeout(700);
}

async function dragCaptionAboveCard(page: Page, fromCaption: string, targetIndex: number): Promise<void> {
  const pts = await page.evaluate(({cap, ti}) => {
    const caps = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name')];
    const src = caps.find((e) => (e.textContent ?? '').trim() === cap) as HTMLElement;
    const cards = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')];
    const target = cards[ti] as HTMLElement;
    if (!src || !target) return null;
    const targetCap = (target.querySelector('.d4-filter-column-name') ?? target) as HTMLElement;
    const sr = src.getBoundingClientRect();
    const tcr = targetCap.getBoundingClientRect();
    return {
      sx: Math.round(sr.x + sr.width / 2), sy: Math.round(sr.y + sr.height / 2),
      tx: Math.round(tcr.x + tcr.width / 2), ty: Math.round(tcr.y + tcr.height / 2),
    };
  }, {cap: fromCaption, ti: targetIndex});
  if (!pts) throw new Error(`dragCaptionAboveCard: caption "${fromCaption}" or card index ${targetIndex} not found`);
  const stripsBefore = await cardStripColours(page);
  expect(stripsBefore.length, `dragCaptionAboveCard(${fromCaption}): the panel exposes no thin card `
    + 'drop strips, so the reorder drag has nothing to wait on and no barrier could be built')
    .toBeGreaterThan(0);
  await page.mouse.move(pts.sx, pts.sy);
  await page.mouse.down();
  for (let i = 1; i <= 8; i++) {
    const x = Math.round(pts.sx + (pts.tx - pts.sx) * i / 8);
    const y = Math.round(pts.sy + (pts.ty - pts.sy) * i / 8);
    await page.mouse.move(x, y, {steps: 3});
    await page.waitForTimeout(30);
  }
  try {
    await page.waitForFunction((before: string[]) => {
      if (!document.body.classList.contains('d4-drag')) return false;
      const now = [...document.querySelectorAll('[name="viewer-Filters"] *')]
        .filter((e) => {
          const r = e.getBoundingClientRect();
          return r.height > 0 && r.height <= 3 && r.width > 8;
        })
        .map((e) => window.getComputedStyle(e).backgroundColor);
      return now.some((c, i) => c !== before[i]);
    }, stripsBefore, {timeout: 5000, polling: 100});
  }
  catch {
    await page.mouse.up();
    throw new Error(`dragCaptionAboveCard(${fromCaption}): the reorder drag never started — within 5s `
      + 'document.body did not carry d4-drag with a card drop strip changing colour. The reorder path '
      + 'is registered with dropIndication:false, so .d4-drop-zone never appears on it and those two '
      + 'are the only indications the product gives.');
  }
  await page.mouse.move(pts.tx, pts.ty, {steps: 4});
  try {
    await page.waitForFunction((before: string[]) => {
      const now = [...document.querySelectorAll('[name="viewer-Filters"] *')]
        .filter((e) => {
          const r = e.getBoundingClientRect();
          return r.height > 0 && r.height <= 3 && r.width > 8;
        })
        .map((e) => window.getComputedStyle(e).backgroundColor);
      return now.some((c, i) => c !== before[i] && /^rgba?\(0, 0, 0/.test(c));
    }, stripsBefore, {timeout: 5000, polling: 100});
  }
  catch {
    await page.mouse.up();
    throw new Error(`dragCaptionAboveCard(${fromCaption}): the pointer reached card index ${targetIndex} `
      + 'but no card drop strip turned black, so the drop point was never armed and a mouseup there '
      + 'would land on nothing.');
  }
  await page.mouse.up();
  await expect.poll(async () => page.evaluate(() => document.body.classList.contains('d4-drag')),
    {timeout: 10_000, intervals: [100, 200, 400],
      message: `dragCaptionAboveCard(${fromCaption}): document.body still carries d4-drag after the `
        + 'release, so the drop never completed'}).toBe(false);
}

async function cardStripColours(page: Page): Promise<string[]> {
  return page.evaluate(() => [...document.querySelectorAll('[name="viewer-Filters"] *')]
    .filter((e) => {
      const r = e.getBoundingClientRect();
      return r.height > 0 && r.height <= 3 && r.width > 8;
    })
    .map((e) => window.getComputedStyle(e).backgroundColor));
}

async function reorderToMatch(page: Page, desiredOrder: string[]): Promise<number> {
  let drags = 0;
  for (let i = 0; i < desiredOrder.length - 1; i++) {
    const current = await orderedCaptions(page);
    if (current[i] === desiredOrder[i]) continue;
    await dragCaptionAboveCard(page, desiredOrder[i], i);
    drags++;
  }
  return drags;
}

async function removeAllViaHamburger(page: Page): Promise<void> {
  await v.drivePanelMenuLeaf(page, 'Filters', null, 'Remove All');
  await expect.poll(async () => cardCount(page),
    {timeout: 20_000, intervals: [300, 600, 1200],
      message: 'the "Remove All" leaf of the panel menu was driven but the panel still carries cards'})
    .toBe(0);
}

async function addViaHeaderCombo(page: Page, column: string): Promise<void> {
  await page.evaluate(() => {
    const header = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header');
    const combo = header?.querySelector('[name="div-column-combobox-"]');
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    const colLabel = combo?.querySelector('.d4-column-selector-column');
    (colLabel || combo)?.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
  });
  await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout: 5000});
  await page.keyboard.press(column[0].toLowerCase());
  await page.waitForFunction(() => !!document.querySelector('input.d4-column-selector-search-input'),
    null, {timeout: 10_000});
  await page.keyboard.press('Control+a');
  await page.keyboard.type(column, {delay: 40});
  expect(await page.evaluate(() =>
    (document.querySelector('input.d4-column-selector-search-input') as HTMLInputElement).value),
  `the header combo's search box does not hold ${column}, so the commit below would add whatever `
    + 'column the fuzzy match happened to leave first').toBe(column);
  await page.keyboard.press('Enter');
  await expect.poll(async () => (await orderedCaptions(page)).includes(column), {
    timeout: 20_000,
    intervals: [500, 1000, 2000],
    message: `no ${column} card came out of the header combo commit`,
  }).toBe(true);
  await expect.poll(async () => page.locator('.d4-column-selector-backdrop').count(), {
    timeout: 10_000,
    message: 'the header combo popup stayed open and would intercept the next gesture',
  }).toBe(0);
}

async function addViaColumnProperties(page: Page, column: string): Promise<void> {
  await page.evaluate((col) => { grok.shell.o = grok.shell.tv.dataFrame.col(col); }, column);
  await page.waitForTimeout(1000);
  await page.evaluate(() => {
    const header = [...document.querySelectorAll('[name="div-section--Filter"]')]
      .find((e) => e.getBoundingClientRect().width > 0) as HTMLElement | undefined;
    const linkVisible = () => [...document.querySelectorAll('label.d4-link-action')]
      .some((e) => (e.textContent ?? '').trim().toLowerCase() === 'add filter' && e.getBoundingClientRect().width > 0);
    if (!linkVisible()) header?.click();
  });
  await page.waitForTimeout(800);
  await page.evaluate(() => {
    const link = [...document.querySelectorAll('label.d4-link-action')]
      .find((e) => (e.textContent ?? '').trim().toLowerCase() === 'add filter' && e.getBoundingClientRect().width > 0) as HTMLElement | undefined;
    link?.click();
  });
  await page.waitForTimeout(900);
}

async function narrowToTopCategory(page: Page, column: string): Promise<number> {
  const c = await page.evaluate(async (col) => {
    const df = grok.shell.tv.dataFrame;
    const fg = grok.shell.tv.getFiltersGroup();
    const counts: Record<string, number> = {};
    const dc = df.col(col);
    for (let i = 0; i < df.rowCount; i++) counts[dc.get(i)] = (counts[dc.get(i)] || 0) + 1;
    const top = Object.entries(counts).sort((a, b) => (b[1] as number) - (a[1] as number))[0][0];
    fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: col, selected: [top]});
    await new Promise((r) => setTimeout(r, 700));
    return df.filter.trueCount;
  }, column);
  expect(typeof c, `${column} narrow: trueCount must be a number`).toBe('number');
  return c;
}

async function openColumnVisibilityDialog(page: Page, column: string): Promise<void> {
  const pt = await gridHeaderPoint(page, column);
  await page.mouse.click(pt.x, pt.y, {button: 'right'});
  await page.waitForTimeout(600);
  await driveOpenMenuLeaf(page, null, 'Order or Hide Columns...');
  await page.locator('.d4-dialog').filter({hasText: 'Order or Hide Columns'}).first()
    .waitFor({timeout: 10_000});
  await page.waitForTimeout(600);
}

const COLUMN_GRID_ROW_PITCH = 28.1;

async function waitForColumnVisibility(page: Page, column: string, expected: boolean): Promise<boolean> {
  const deadline = Date.now() + 8000;
  let visible = await gridColumnVisible(page, column);
  while (visible !== expected && Date.now() < deadline) {
    await page.waitForTimeout(200);
    visible = await gridColumnVisible(page, column);
  }
  return visible;
}

async function setColumnVisibilityInDialog(page: Page, column: string, expected: boolean): Promise<void> {
  const before = await gridColumnVisible(page, column);
  expect(before, `${column}: the dialog check cell is a toggle, so it can only be clicked from the `
    + `opposite state; wanted ${!expected} before the click but the grid reports ${before}`).toBe(!expected);

  const search = page.locator('.d4-dialog input.d4-search-input').first();
  await search.click();
  await search.fill('');
  await page.keyboard.type(column, {delay: 40});
  await page.waitForFunction((c) => (document.querySelector('.d4-dialog input.d4-search-input') as
    HTMLInputElement | null)?.value === c, column, {timeout: 5000, polling: 100});
  await page.evaluate(() => new Promise<void>((r) =>
    requestAnimationFrame(() => requestAnimationFrame(() => r()))));

  const pt = await page.evaluate((pitch) => {
    const dlg = document.querySelector('.d4-dialog');
    const overlay = dlg?.querySelector('[name="viewer-Grid"] [name="overlay"]');
    const allCheck = dlg?.querySelector('input[type="checkbox"]');
    if (!overlay || !allCheck) return null;
    const ov = overlay.getBoundingClientRect();
    const acb = allCheck.getBoundingClientRect();
    return {x: acb.left + acb.width / 2, y: ov.top + pitch / 2};
  }, COLUMN_GRID_ROW_PITCH);
  expect(pt, `${column}: the Order or Hide Columns dialog exposes no ColumnGrid overlay and `
    + 'allCheck pair to derive the check-cell click point from').not.toBeNull();

  await page.mouse.click(pt!.x, pt!.y);

  const after = await waitForColumnVisibility(page, column, expected);
  expect(after, `${column}: the check-cell click at (${pt!.x.toFixed(1)}, ${pt!.y.toFixed(1)}) did `
    + `not land — visible was ${before} before the click and is still ${after} after it`).toBe(expected);
}

async function closeColumnVisibilityDialog(page: Page): Promise<void> {
  await page.evaluate(() => {
    const dlg = [...document.querySelectorAll('.d4-dialog')]
      .find((d) => (d.textContent ?? '').includes('Order or Hide Columns'));
    const close = [...(dlg?.querySelectorAll('button, .ui-btn') ?? [])]
      .find((b) => (b.textContent ?? '').trim().toUpperCase() === 'CLOSE') as HTMLElement | undefined;
    close?.click();
  });
  await expect.poll(async () => page.locator('.d4-dialog').filter({hasText: 'Order or Hide Columns'}).count(),
    {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(0);
}

async function gridColumnVisible(page: Page, column: string): Promise<boolean> {
  return page.evaluate((c: string) => !!grok.shell.tv.grid.columns.byName(c)?.visible, column);
}

async function closePanel(page: Page): Promise<void> {
  await page.locator('[name="viewer-Filters"]').first().hover();
  let clicked = false;
  for (const icon of ['icon-times', 'Close']) {
    try { await v.clickViewerTitlebarIcon(page, 'Filters', icon); clicked = true; break; } catch (_) { /* next */ }
  }
  expect(clicked, 'the Filter Panel exposes no title-bar close control').toBe(true);
  await expect.poll(async () => page.locator('[name="viewer-Filters"]').count(),
    {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(0);
}

async function reopenPanel(page: Page): Promise<void> {
  await page.locator('.d4-ribbon-panel [name="icon-filter"]').first().click();
  await page.locator('[name="viewer-Filters"]').first().waitFor({timeout: 15_000});
  await expect.poll(async () => cardCount(page),
    {timeout: 20_000, intervals: [400, 800, 1500]}).toBeGreaterThan(0);
}

async function columnNames(page: Page): Promise<string[]> {
  return page.evaluate(() => grok.shell.tv.dataFrame.columns.names());
}

async function pickerCheckedCount(page: Page): Promise<number> {
  return page.evaluate((sel) => {
    const dlg = document.querySelector(sel);
    if (!dlg) return -1;
    for (const label of dlg.querySelectorAll('label')) {
      const m = /^(\d+) checked$/.exec((label.textContent ?? '').trim());
      if (m) return Number(m[1]);
    }
    return -1;
  }, SELECT_COLUMNS_DIALOG);
}

async function dialogGone(page: Page, selector: string): Promise<void> {
  await expect.poll(async () => page.locator(selector).count(),
    {timeout: 15_000, intervals: [200, 400, 800],
      message: `the dialog ${selector} is still on screen, so its OK click did not commit`}).toBe(0);
}

async function multiValueFilterColumns(page: Page, filterType: string): Promise<string[]> {
  return page.evaluate((ft) => grok.shell.tv.getFiltersGroup().filters
    .filter((f: any) => f.filterType === ft)
    .map((f: any) => String(f.columnName ?? f.filterColumnName ?? '')), filterType);
}

const MVF_COLUMN = 'MVF_TOKENS';
const MVF_SEPARATOR = ';';
const MVF_TOKEN_A = 'red';
const MVF_TOKEN_B = 'blue';
const MVF_PATTERN: string[][] = [
  ['red', 'green'], ['red'], ['green', 'blue'], ['red', 'green', 'blue'], ['blue'],
  ['red', 'amber'], ['green'], ['red', 'green'], ['amber'], ['blue', 'amber'],
];
const MVF_HEADER_BAND_H = 10;
const MVF_ROW_PITCH = 27;
const MVF_ROW_CENTRE_OFFSET = 13;
const MVF_X_CHECKBOX = 10;

interface TokenDistribution {
  rowCount: number;
  emptyRows: number;
  counts: Record<string, number>;
  tokens: string[];
  pairAnd: number;
  pairOr: number;
}

interface MultiValueCardState { include: string[]; mode: string; trueCount: number; }

async function addMultiValueFixtureColumn(page: Page): Promise<void> {
  await page.evaluate(({name, pattern, sep}) => {
    const df = grok.shell.tv.dataFrame;
    const values: string[] = [];
    for (let i = 0; i < df.rowCount; i++) values.push(pattern[i % pattern.length].join(sep));
    df.columns.add(DG.Column.fromStrings(name, values));
  }, {name: MVF_COLUMN, pattern: MVF_PATTERN, sep: MVF_SEPARATOR});
  await expect.poll(async () => columnNames(page),
    {timeout: 20_000, intervals: [200, 400, 800],
      message: `the "${MVF_COLUMN}" fixture column never appeared on the table, so the dialog below `
        + 'would have nothing multi-value to be pointed at'}).toContain(MVF_COLUMN);
}

async function removeMultiValueFixture(page: Page): Promise<void> {
  await page.evaluate(({name, ft}) => {
    const fg = grok.shell.tv.getFiltersGroup();
    for (const f of [...fg.filters])
      if (f.filterType === ft) fg.remove(f);
    const df = grok.shell.tv.dataFrame;
    if (df.col(name)) df.columns.remove(name);
  }, {name: MVF_COLUMN, ft: MULTI_VALUE_FILTER_TYPE});
}

async function appliedRowFilters(page: Page): Promise<string[]> {
  return page.evaluate(() => [...(grok.shell.tv.dataFrame.rows.filters ?? [])].map((s: any) => String(s)));
}

async function openCleanDemogView(page: Page, path: string): Promise<void> {
  await page.evaluate(async (p) => {
    const previous = grok.shell.v;
    const df = await grok.dapi.files.readCsv(p);
    grok.shell.addTableView(df);
    await new Promise((r) => setTimeout(r, 800));
    previous?.close();
    await new Promise((r) => setTimeout(r, 800));
    grok.shell.tv.getFiltersGroup();
  }, path);
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 30_000});
  await expect.poll(async () => page.evaluate(() =>
    document.querySelectorAll('[name="viewer-Filters"]').length),
  {timeout: 20_000, intervals: [300, 600, 1200],
    message: 'the previous table view must be gone, or every [name="viewer-Filters"] reading below '
      + 'would span two panels at once'}).toBe(1);
  await removeAllViaHamburger(page);
  await expect.poll(async () => cardCount(page),
    {timeout: 20_000, intervals: [300, 600, 1200],
      message: 'the fresh table view\'s panel must be emptied before the Multi Value add, or the '
        + 'card produced below would not be attributable'}).toBe(0);
  const applied = await appliedRowFilters(page);
  expect(applied, `the fresh demog view must start with nothing filtering it; dataFrame.rows.filters `
    + `reads [${applied.join(', ')}], so every row count below would be measuring that filter too`)
    .toEqual([]);
  await expect.poll(async () => trueCount(page),
    {timeout: 20_000, intervals: [300, 600, 1200],
      message: `the fresh demog view must show all ${fullRowCount} rows before any token is ticked`})
    .toBe(fullRowCount);
}

async function tokenDistribution(page: Page): Promise<TokenDistribution> {
  const d = await page.evaluate(({name, sep, a, b}) => {
    const col = grok.shell.tv.dataFrame.col(name);
    if (!col) return null;
    const counts: Record<string, number> = {};
    let emptyRows = 0;
    let pairAnd = 0;
    let pairOr = 0;
    for (let i = 0; i < col.length; i++) {
      const tokens = [...new Set(String(col.get(i) ?? '').split(sep).filter((t: string) => t !== ''))];
      if (tokens.length === 0) emptyRows++;
      for (const t of tokens) counts[t] = (counts[t] ?? 0) + 1;
      const hasA = tokens.indexOf(a) >= 0;
      const hasB = tokens.indexOf(b) >= 0;
      if (hasA && hasB) pairAnd++;
      if (hasA || hasB) pairOr++;
    }
    return {rowCount: col.length, emptyRows, counts, tokens: Object.keys(counts).sort(), pairAnd, pairOr};
  }, {name: MVF_COLUMN, sep: MVF_SEPARATOR, a: MVF_TOKEN_A, b: MVF_TOKEN_B});
  expect(d, `the "${MVF_COLUMN}" column is missing, so no expectation could be derived from its own `
    + 'values and every count below would be comparing against nothing').not.toBeNull();
  return d!;
}

async function multiValueCardState(page: Page): Promise<MultiValueCardState | null> {
  return page.evaluate(({col, ft}) => {
    const f = grok.shell.tv.getFiltersGroup().filters
      .find((x: any) => x.filterType === ft && (x.columnName ?? x.filterColumnName) === col);
    if (!f) return null;
    const st = f.saveState();
    return {
      include: (st.include ?? []).map((s: any) => String(s)),
      mode: String(st.mode ?? ''),
      trueCount: grok.shell.tv.dataFrame.filter.trueCount as number,
    };
  }, {col: MVF_COLUMN, ft: MULTI_VALUE_FILTER_TYPE});
}

async function clickMultiValueRow(page: Page, rowIndex: number): Promise<MultiValueCardState> {
  const res = await page.evaluate(async ({col, ft, ri, cx, band, pitch, centre}) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col);
    const overlay = card?.querySelector('[name="viewer-Grid"] [name="overlay"]') as HTMLElement | null;
    if (!overlay) return null;
    const read = () => {
      const f = grok.shell.tv.getFiltersGroup().filters
        .find((x: any) => x.filterType === ft && (x.columnName ?? x.filterColumnName) === col);
      if (!f) return null;
      const st = f.saveState();
      return {
        include: (st.include ?? []).map((s: any) => String(s)),
        mode: String(st.mode ?? ''),
        trueCount: grok.shell.tv.dataFrame.filter.trueCount as number,
      };
    };
    const before = read();
    if (!before) return null;
    const key = (s: any) => `${s.trueCount}|${s.mode}|${[...s.include].sort().join(',')}`;
    const rect = overlay.getBoundingClientRect();
    const o = {
      bubbles: true, cancelable: true, view: window, button: 0,
      clientX: rect.left + cx, clientY: rect.top + band + pitch * ri + centre,
    };
    overlay.dispatchEvent(new MouseEvent('mousedown', o));
    overlay.dispatchEvent(new MouseEvent('mouseup', o));
    overlay.dispatchEvent(new MouseEvent('click', o));
    let after = before;
    let prev = key(before);
    for (let waited = 0; waited < 10_000; waited += 100) {
      await new Promise((r) => setTimeout(r, 100));
      after = read()!;
      const k = key(after);
      if (k !== key(before) && k === prev) break;
      prev = k;
    }
    return {before, after};
  }, {col: MVF_COLUMN, ft: MULTI_VALUE_FILTER_TYPE, ri: rowIndex, cx: MVF_X_CHECKBOX,
    band: MVF_HEADER_BAND_H, pitch: MVF_ROW_PITCH, centre: MVF_ROW_CENTRE_OFFSET});
  expect(res, `the "${MVF_COLUMN}" card exposes no [name="overlay"] canvas inside its `
    + `[name="viewer-Grid"] body, or no '${MULTI_VALUE_FILTER_TYPE}' filter to read back, so the `
    + `row-${rowIndex} tick could neither be dispatched nor measured`).not.toBeNull();
  return res!.after;
}

async function toggleMultiValueMode(page: Page): Promise<{text: string; state: MultiValueCardState}> {
  const res = await page.evaluate(async ({col, ft}) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col);
    const header = card?.querySelector('.d4-filter-header');
    const el = [...(header?.children ?? [])]
      .find((e) => ['AND', 'OR'].indexOf((e.textContent ?? '').trim()) >= 0) as HTMLElement | undefined;
    if (!el) return null;
    const read = () => {
      const f = grok.shell.tv.getFiltersGroup().filters
        .find((x: any) => x.filterType === ft && (x.columnName ?? x.filterColumnName) === col);
      if (!f) return null;
      const st = f.saveState();
      return {
        include: (st.include ?? []).map((s: any) => String(s)),
        mode: String(st.mode ?? ''),
        trueCount: grok.shell.tv.dataFrame.filter.trueCount as number,
      };
    };
    const before = read();
    if (!before) return null;
    el.click();
    let after = before;
    let prev = `${before.mode}|${before.trueCount}`;
    for (let waited = 0; waited < 10_000; waited += 100) {
      await new Promise((r) => setTimeout(r, 100));
      after = read()!;
      const k = `${after.mode}|${after.trueCount}`;
      if (after.mode !== before.mode && k === prev) break;
      prev = k;
    }
    return {text: (el.textContent ?? '').trim(), before, after};
  }, {col: MVF_COLUMN, ft: MULTI_VALUE_FILTER_TYPE});
  expect(res, `the "${MVF_COLUMN}" card header carries no AND/OR mode control to click, so the mode `
    + 'half of this step could not be driven at all').not.toBeNull();
  return {text: res!.text, state: res!.after};
}

async function boolColumnTrueCount(page: Page, column: string): Promise<number> {
  return page.evaluate((c) => {
    const col = grok.shell.tv.dataFrame.col(c);
    if (!col || col.type !== 'bool') return -1;
    let n = 0;
    for (let i = 0; i < col.length; i++) if (col.get(i) === true) n++;
    return n;
  }, column);
}

async function rebuildPanelFromColumnVisibility(page: Page): Promise<void> {
  await removeAllViaHamburger(page);
  const left = await cardCount(page);
  expect(left, `Remove All left ${left} card(s), so the panel's saved card set is not empty and `
    + 'the reopen below would restore it instead of rebuilding from column visibility').toBe(0);
  await closePanel(page);
  await reopenPanel(page);
}

test('Filter Panel — Add, Reorder, and Remove Entry Points', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  await softStep('Setup — the zero-card baseline: empty the Filter Panel and confirm it carries no cards', async () => {
    await removeAllViaHamburger(page);
    expect(await cardCount(page)).toBe(0);
    expect(await trueCount(page)).toBe(fullRowCount);
    await expectHeaderCounterQuietNow(page,
      'the panel was emptied and nothing filters, so the counter must be hidden or read 0');
  });

  await softStep('Scenario 1 Step 1 — add path (a) header drag: RACE header onto the panel', async () => {
    await dragColumnHeaderToPanel(page, 'RACE');
    const captions = await orderedCaptions(page);
    expect(captions[0], 'RACE card must appear first after the header drag').toBe('RACE');
    expect(await cardCount(page)).toBe(1);
    expect(await trueCount(page)).toBe(fullRowCount);
    await expectHeaderCounterQuietNow(page,
      'a freshly dragged-in RACE card is not filtering yet, so the counter must be hidden or read 0');
  });

  await softStep('Scenario 1 Step 2 — GROK-19516 pinned header drag: pin AGE, drag it onto the panel', async () => {
    const pinBefore = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return {idx: grid.columns.byName('AGE')?.idx, frozen: grid.props.frozenColumns};
    });
    expect(typeof pinBefore.idx, 'the grid exposes no AGE column to pin').toBe('number');
    expect(typeof pinBefore.frozen, 'the grid exposes no frozenColumns count').toBe('number');
    expect(pinBefore.idx < pinBefore.frozen, 'AGE is already pinned, so the drag below would not be a pinned-source drag').toBe(false);
    const agePt = await gridHeaderPoint(page, 'AGE');
    await page.mouse.click(agePt.x, agePt.y, {button: 'right'});
    await page.waitForTimeout(600);
    await driveOpenMenuLeaf(page, 'Pin', 'Pin Column');
    await page.waitForTimeout(700);
    const pinAfter = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return {idx: grid.columns.byName('AGE')?.idx, frozen: grid.props.frozenColumns};
    });
    expect(typeof pinAfter.idx, 'the AGE grid column disappeared while pinning').toBe('number');
    expect(pinAfter.idx < pinAfter.frozen,
      `AGE was not pinned by the context menu (idx ${pinAfter.idx}, frozenColumns ${pinAfter.frozen})`).toBe(true);
    await dragColumnHeaderToPanel(page, 'AGE');
    const captions = await orderedCaptions(page);
    expect(captions[0], 'AGE card must appear first (GROK-19516: pinned drag must produce a card)').toBe('AGE');
    expect(captions[1]).toBe('RACE');
    expect(await cardCount(page)).toBe(2);
    expect(await trueCount(page)).toBe(fullRowCount);
    await expectHeaderCounterQuietNow(page,
      'the pinned AGE drag added a card that is not filtering, so the counter must be hidden or read 0');
  });

  await softStep('Scenario 1 Step 3 — add path (b) panel header picker adds SEX first', async () => {
    await addViaHeaderCombo(page, 'SEX');
    const captions = await orderedCaptions(page);
    expect(captions[0]).toBe('SEX');
    expect(await cardCount(page)).toBe(3);
    expect(await trueCount(page)).toBe(fullRowCount);
    await expectHeaderCounterQuietNow(page,
      'the header picker added SEX without a criterion, so the counter must be hidden or read 0');
  });

  await softStep('Scenario 1 Step 4 — add path (c) column properties Add filter link adds DIS_POP first', async () => {
    await addViaColumnProperties(page, 'DIS_POP');
    const captions = await orderedCaptions(page);
    expect(captions[0], 'DIS_POP card must appear first via the column-properties Add filter link').toBe('DIS_POP');
    expect(await cardCount(page)).toBe(4);
    expect(await trueCount(page)).toBe(fullRowCount);
    await expectHeaderCounterQuietNow(page,
      'the column-properties link added DIS_POP without a criterion, so the counter must be hidden or read 0');
  });

  await softStep('Scenario 1 Step 5 — add path (d) context menu Expression adds first', async () => {
    await v.drivePanelMenuLeaf(page, 'Filters', 'Add Filter', 'Expression');
    await expect.poll(async () => orderedCaptions(page),
      {timeout: 20_000, intervals: [300, 600, 1200],
        message: 'the Add Filter > Expression leaf was driven but no Expression card appeared'})
      .toContain('Expression');
    const captions = await orderedCaptions(page);
    expect(captions[0]).toBe('Expression');
    expect(await cardCount(page)).toBe(5);
    const hasExprCard = await page.evaluate(() =>
      !!document.querySelector('[name="viewer-Filters"] .d4-filter .d4-expression-filter'));
    expect(hasExprCard).toBe(true);
    expect(await trueCount(page)).toBe(fullRowCount);
    await expectHeaderCounterQuietNow(page,
      'the Expression card carries no criterion yet, so the counter must be hidden or read 0');
  });

  let trueCountRace = 0;

  await softStep('Scenario 1 Step 6 — RACE card starts filtering', async () => {
    trueCountRace = await narrowToTopCategory(page, 'RACE');
    expect(trueCountRace).toBeGreaterThan(0);
    expect(trueCountRace).toBeLessThan(fullRowCount);
    await expectHeaderCounter(page, '1',
      'the RACE card is the only one filtering, so the counter must read 1');
    expect(await cardCount(page)).toBe(5);
  });

  await softStep('Scenario 1 Step 7 — GROK-18765 disable of a column-properties card (DIS_POP)', async () => {
    const trueCountRaceDispop = await narrowToTopCategory(page, 'DIS_POP');
    expect(trueCountRaceDispop).toBeLessThan(trueCountRace);
    await expectHeaderCounter(page, '2',
      'RACE and DIS_POP are both filtering, so the counter must read 2');

    const disabled = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const cards = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')] as HTMLElement[];
      const card = cards.find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === 'DIS_POP')!;
      const cb = (card.querySelector('input[type="checkbox"].ui-input-editor')
        ?? card.querySelector('input[type="checkbox"]')) as HTMLInputElement;
      cb.click();
      await new Promise((r) => setTimeout(r, 700));
      const f = fg.filters.find((x: any) => (x.filterColumnName ?? x.columnName) === 'DIS_POP');
      return {
        trueCount: grok.shell.tv.dataFrame.filter.trueCount,
        isFiltering: f?.isFiltering,
        hasDisabledClass: card.classList.contains('d4-filter-disabled'),
        checked: cb.checked,
      };
    });
    expect(disabled.trueCount, 'disabling DIS_POP must return trueCount to the RACE-only value').toBe(trueCountRace);
    expect(disabled.isFiltering).toBe(false);
    expect(disabled.hasDisabledClass).toBe(true);
    expect(disabled.checked).toBe(false);
    await expectHeaderCounter(page, '1',
      'the DIS_POP card was disabled, so only RACE is left filtering and the counter must read 1');
    expect(await cardCount(page)).toBe(5);

    const rechecked = await page.evaluate(async () => {
      const cards = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')] as HTMLElement[];
      const card = cards.find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === 'DIS_POP')!;
      const cb = (card.querySelector('input[type="checkbox"].ui-input-editor')
        ?? card.querySelector('input[type="checkbox"]')) as HTMLInputElement;
      if (!cb.checked) cb.click();
      await new Promise((r) => setTimeout(r, 700));
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(rechecked).toBe(trueCountRaceDispop);
    await expectHeaderCounter(page, '2',
      're-checking DIS_POP puts two cards back into filtering, so the counter must read 2');
  });

  const consoleErrors: string[] = [];
  page.on('console', (msg) => { if (msg.type() === 'error') consoleErrors.push(msg.text()); });
  const AMBIENT = 'Permissions policy violation: compute-pressure';
  const errorSet = () => new Set(consoleErrors.filter((t) => !t.includes(AMBIENT)));

  let orderBeforeReorder: string[] = [];
  let trueCountBeforeReorder = 0;
  let errorsBeforeReorder = new Set<string>();

  await softStep('Scenario 2 Step 1 — GROK-12955 reorder drag: bottom card to the top', async () => {
    orderBeforeReorder = await orderedCaptions(page);
    trueCountBeforeReorder = await trueCount(page);
    errorsBeforeReorder = errorSet();
    const bottom = orderBeforeReorder[orderBeforeReorder.length - 1];

    await dragCaptionAboveCard(page, bottom, 0);

    expect(await page.evaluate(() => 1 + 1)).toBe(2);

    const orderAfter = await orderedCaptions(page);
    const expectedOrder = [bottom, ...orderBeforeReorder.filter((c) => c !== bottom)];
    expect(orderAfter, 'reorder must move the bottom caption to position 0, others keep relative order').toEqual(expectedOrder);
    expect(await cardCount(page)).toBe(orderBeforeReorder.length);
    expect(await trueCount(page), 'reorder must not change filtering').toBe(trueCountBeforeReorder);

    const newErrors = [...errorSet()].filter((t) => !errorsBeforeReorder.has(t));
    expect(newErrors, `no new console error text after reorder; got: ${newErrors.join(' | ')}`).toEqual([]);
  });

  await softStep('Scenario 2 Step 2 — reorder revert', async () => {
    const drags = await reorderToMatch(page, orderBeforeReorder);
    expect(drags, 'the revert completed without dragging anything, so the order check below would '
      + 'pass on a no-op rather than on a restored order').toBeGreaterThan(0);

    const orderAfter = await orderedCaptions(page);
    expect(orderAfter, 'revert must restore the exact pre-Step-1 caption order').toEqual(orderBeforeReorder);
    expect(await cardCount(page)).toBe(orderBeforeReorder.length);
    expect(await trueCount(page), 'revert must not change filtering').toBe(trueCountBeforeReorder);
    const newErrors = [...errorSet()].filter((t) => !errorsBeforeReorder.has(t));
    expect(newErrors, `no new console error text after revert; got: ${newErrors.join(' | ')}`).toEqual([]);
  });

  await softStep('Scenario 3 Step 1 — remove via card X icon (Expression card)', async () => {
    const before = await cardCount(page);
    const beforeCount = await trueCount(page);
    await page.evaluate(async () => {
      const cards = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')] as HTMLElement[];
      const exprCard = cards.find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === 'Expression')!;
      (exprCard.querySelector('[name="icon-times"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 600));
    });
    const captions = await orderedCaptions(page);
    expect(captions).not.toContain('Expression');
    expect(await cardCount(page)).toBe(before - 1);
    expect(await trueCount(page), 'the Expression card was not filtering, so removing it cannot '
      + 'move the row count').toBe(beforeCount);
  });

  await softStep('Scenario 3 Step 2 — remove others keeps the actively-filtering cards', async () => {
    const before = await trueCount(page);
    const captionsBefore = await orderedCaptions(page);
    expect(captionsBefore).toContain('RACE');
    expect(captionsBefore).toContain('DIS_POP');

    const counterOpened = await page.evaluate(() => {
      const counter = document.querySelector(
        '[name="viewer-Filters"] .d4-filter-group-header .d4-filter-indicator') as HTMLElement | null;
      if (!counter) return false;
      counter.click();
      return true;
    });
    expect(counterOpened, 'the panel header carries no active-filter counter to open the '
      + '"Remove others" menu from').toBe(true);
    const cardsBefore = await cardCount(page);
    await driveOpenMenuLeaf(page, null, 'Remove others');
    await expect.poll(async () => cardCount(page),
      {timeout: 20_000, intervals: [300, 600, 1200],
        message: `"Remove others" was driven but the panel still carries all ${cardsBefore} cards`})
      .toBeLessThan(cardsBefore);

    const captionsAfter = await orderedCaptions(page);
    expect(captionsAfter, 'RACE stays after Remove others (it is filtering)').toContain('RACE');
    expect(captionsAfter, 'DIS_POP stays after Remove others (it is filtering)').toContain('DIS_POP');
    expect(captionsAfter).not.toContain('AGE');
    expect(captionsAfter).not.toContain('SEX');
    await holdTrueCount(page, before, 'removing non-filtering cards cannot move the filtered row '
      + 'count, and it must not move over the seconds after the removal either');
  });

  await softStep('Scenario 3 Step 3 — Remove All empties the panel with no confirmation dialog', async () => {
    await removeAllViaHamburger(page);
    expect(await page.evaluate(() =>
      document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name').length)).toBe(0);
    expect(await cardCount(page), 'Remove All left a card behind that carries no caption, so the '
      + 'caption count above cannot see it').toBe(0);
    expect(await trueCount(page)).toBe(fullRowCount);
    await expectHeaderCounterQuiet(page,
      'Remove All took every card away, so nothing filters and the counter must be hidden or read 0');
    expect(await page.evaluate(() => !!document.querySelector('.d4-dialog [name="button-OK"]'))).toBe(false);
  });

  const HIDDEN = ['RACE', 'AGE'];
  let persistedTrueCount = 0;

  await softStep('Scenario 4 Step 1 — the cards exist before the columns are hidden: the panel auto-creates a card per eligible column', async () => {
    await rebuildPanelFromColumnVisibility(page);
    const captions = await orderedCaptions(page);
    for (const c of HIDDEN)
      expect(captions, `${c} must have a card before it is hidden`).toContain(c);
    expect(await trueCount(page)).toBe(fullRowCount);
  });

  await softStep('Scenario 4 Step 2 — the hide gesture really takes: hide the two columns through the Order or Hide Columns dialog', async () => {
    for (const column of HIDDEN) {
      await openColumnVisibilityDialog(page, column);
      await setColumnVisibilityInDialog(page, column, false);
      await closeColumnVisibilityDialog(page);
      expect(await gridColumnVisible(page, column),
        `${column} is still visible in the grid — the dialog gesture did not take`).toBe(false);
    }
  });

  await softStep('Scenario 4 Step 3 — hidden columns have no cards, the others are untouched', async () => {
    await rebuildPanelFromColumnVisibility(page);
    const captions = await orderedCaptions(page);
    for (const c of HIDDEN)
      expect(captions, `${c} is hidden in the grid, so it must have no filter card`).not.toContain(c);
    expect(captions, 'hiding two columns emptied the whole panel').not.toEqual([]);
    expect(captions, 'SEX was not hidden and must still have a card').toContain('SEX');
    expect(await trueCount(page), 'hiding a column must not filter any row').toBe(fullRowCount);
  });

  await softStep('Scenario 4 Step 4 — showing the columns again restores their cards', async () => {
    for (const column of HIDDEN) {
      await openColumnVisibilityDialog(page, 'SEX');
      await setColumnVisibilityInDialog(page, column, true);
      await closeColumnVisibilityDialog(page);
      expect(await gridColumnVisible(page, column),
        `${column} did not become visible again`).toBe(true);
    }
    await rebuildPanelFromColumnVisibility(page);
    const captions = await orderedCaptions(page);
    for (const c of HIDDEN)
      expect(captions, `${c} is visible again, so the rebuilt panel must carry its card; got `
        + `[${captions.join(', ')}]`).toContain(c);
    expect(await trueCount(page)).toBe(fullRowCount);
  });

  await softStep('Scenario 5 Step 1 — three distinguishable card states: narrow RACE, disable AGE, remove SEX', async () => {
    const filtered = await narrowToTopCategory(page, 'RACE');
    expect(filtered).toBeGreaterThan(0);
    expect(filtered).toBeLessThan(fullRowCount);
    persistedTrueCount = filtered;

    const state = await page.evaluate(async () => {
      const cards = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')] as HTMLElement[];
      const byName = (n: string) => cards.find((c) =>
        (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === n);
      const age = byName('AGE')!;
      const cb = (age.querySelector('input[type="checkbox"].ui-input-editor')
        ?? age.querySelector('input[type="checkbox"]')) as HTMLInputElement;
      cb.click();
      await new Promise((r) => setTimeout(r, 700));
      (byName('SEX')!.querySelector('[name="icon-times"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 700));
      return {ageChecked: cb.checked, count: grok.shell.tv.dataFrame.filter.trueCount};
    });
    expect(state.ageChecked, 'the AGE card did not become unchecked').toBe(false);
    expect(state.count, 'disabling AGE and removing SEX must not move the RACE-only count')
      .toBe(persistedTrueCount);
    expect(await orderedCaptions(page)).not.toContain('SEX');
    await expectHeaderCounter(page, '1',
      'only the RACE card is filtering after the AGE disable and the SEX removal, so the counter must read 1');
  });

  await softStep('Scenario 5 Step 2 — closing the panel releases its filtering', async () => {
    await closePanel(page);
    await expect.poll(async () => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount),
      {timeout: 10_000, intervals: [300, 600, 1200]}).toBe(fullRowCount);
  });

  await softStep('Scenario 5 Step 3 — reopening restores per-card state: criterion, disabled card and removal', async () => {
    await reopenPanel(page);
    const captions = await orderedCaptions(page);
    expect(captions, 'the RACE card must come back').toContain('RACE');
    expect(captions, 'the AGE card was disabled, not removed, so it must come back').toContain('AGE');
    expect(captions, 'the SEX card was removed and must NOT be recreated on reopen').not.toContain('SEX');

    const restored = await page.evaluate(() => {
      const cards = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')] as HTMLElement[];
      const age = cards.find((c) =>
        (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === 'AGE');
      const cb = (age?.querySelector('input[type="checkbox"].ui-input-editor')
        ?? age?.querySelector('input[type="checkbox"]')) as HTMLInputElement | undefined;
      return {ageCheckboxFound: !!cb, ageChecked: cb?.checked, count: grok.shell.tv.dataFrame.filter.trueCount};
    });
    expect(restored.ageCheckboxFound, 'the reopened AGE card carries no enable/disable checkbox').toBe(true);
    expect(restored.ageChecked, 'the AGE card came back enabled — its disabled state was not kept').toBe(false);
    expect(restored.count, 'the RACE criterion did not come back with its exact row count')
      .toBe(persistedTrueCount);
    await expectHeaderCounter(page, '1',
      'the reopened panel restored the RACE criterion alone, so the counter must read 1');
  });

  let baselineTrueCount = -1;

  await softStep('Scenario 6 Step 1 — Select Columns... drives the card set both ways', async () => {
    await removeAllViaHamburger(page);
    expect(await cardCount(page), 'the picker is measured from an empty panel, so Remove All must '
      + 'leave zero cards or the "All" result below would not be attributable').toBe(0);
    // Remove All drops the cards first and RELEASES their filters a moment later — reading
    // the baseline immediately captures the still-filtered count, which then never matches
    // the (correctly unfiltered) count the picker leaves behind
    // Remove All empties the PANEL but does not release a criterion the previous scenario
    // restored — the dataframe filter survives with zero cards present, so the baseline
    // must be taken from an explicitly cleared filter or the picker's no-op reads as a change
    await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.filter.setAll(true));
    baselineTrueCount = await v.pollStable(() => trueCount(page), (a, b) => a === b, 5000, 200);
    expect(baselineTrueCount, 'the row count at the start of this scenario must be a positive '
      + `number to compare the picker's effect against; got ${baselineTrueCount}`).toBeGreaterThan(0);

    const pickedColumns = await columnNames(page);
    expect(pickedColumns.length, `the table exposes no columns for the picker to offer; got `
      + `[${pickedColumns.join(', ')}]`).toBeGreaterThan(0);

    await v.drivePanelMenuLeaf(page, 'Filters', null, 'Select Columns...');
    await page.locator(SELECT_COLUMNS_DIALOG).first().waitFor({timeout: 15_000});
    expect(await pickerCheckedCount(page), 'the picker must open showing the panel it was opened '
      + 'from — 0 checked for an empty panel; -1 means no "N checked" label was found at all').toBe(0);

    await page.locator(`${SELECT_COLUMNS_DIALOG} [name="label-All"]`).click();
    await expect.poll(async () => pickerCheckedCount(page),
      {timeout: 10_000, intervals: [200, 400, 800],
        message: `the "All" link must check every one of the ${pickedColumns.length} columns `
          + `[${pickedColumns.join(', ')}]`}).toBe(pickedColumns.length);

    await page.locator(`${SELECT_COLUMNS_DIALOG} [name="button-OK"]`).click();
    await dialogGone(page, SELECT_COLUMNS_DIALOG);
    await expect.poll(async () => cardCount(page),
      {timeout: 20_000, intervals: [300, 600, 1200],
        message: `the picker committed ${pickedColumns.length} checked columns, so the panel must `
          + 'carry that many cards'}).toBe(pickedColumns.length);

    const captions = await orderedCaptions(page);
    expect([...captions].sort(), 'every picked column must get a card and nothing else may appear; '
      + `got [${captions.join(', ')}] for columns [${pickedColumns.join(', ')}]`)
      .toEqual([...pickedColumns].sort());
    await expect.poll(async () => trueCount(page),
      {timeout: 20_000, intervals: [300, 600, 1200],
        message: `the ${pickedColumns.length} cards the picker committed are all unconfigured, so `
          + `choosing which columns get cards must leave the ${baselineTrueCount} shown rows alone`})
      .toBe(baselineTrueCount);

    await v.drivePanelMenuLeaf(page, 'Filters', null, 'Select Columns...');
    await page.locator(SELECT_COLUMNS_DIALOG).first().waitFor({timeout: 15_000});
    expect(await pickerCheckedCount(page), 'reopened, the picker must read back the card set it '
      + `just produced (${pickedColumns.length} checked)`).toBe(pickedColumns.length);

    await page.locator(`${SELECT_COLUMNS_DIALOG} [name="label-None"]`).click();
    await expect.poll(async () => pickerCheckedCount(page),
      {timeout: 10_000, intervals: [200, 400, 800],
        message: 'the "None" link must clear every check'}).toBe(0);

    await page.locator(`${SELECT_COLUMNS_DIALOG} [name="button-OK"]`).click();
    await dialogGone(page, SELECT_COLUMNS_DIALOG);
    await expect.poll(async () => cardCount(page),
      {timeout: 20_000, intervals: [300, 600, 1200],
        message: 'committing an empty check set must take every card away again'}).toBe(0);
    expect(await orderedCaptions(page), 'the panel must be empty after the "None" commit').toEqual([]);
    await expect.poll(async () => trueCount(page),
      {timeout: 20_000, intervals: [300, 600, 1200],
        message: `taking every card away through the picker must leave the ${baselineTrueCount} `
          + 'shown rows alone — none of those cards was filtering'}).toBe(baselineTrueCount);
  });

  await softStep('Scenario 6 Step 2 — Add Filter > Multi Value... really filters a multi-value column', async () => {
    expect(await cardCount(page), 'the Multi Value add is measured from the empty panel Step 1 left').toBe(0);
    const before = await multiValueFilterColumns(page, MULTI_VALUE_FILTER_TYPE);
    expect(before, `the panel already carries multi-value filters [${before.join(', ')}], so a new `
      + 'one would not be attributable to this gesture').toEqual([]);

    await openCleanDemogView(page, datasetPath);

    const columnsBefore = await columnNames(page);
    expect(columnsBefore, `the table already carries a "${MVF_COLUMN}" column, so the fixture below `
      + `would not be this step's own; columns: [${columnsBefore.join(', ')}]`).not.toContain(MVF_COLUMN);

    try {
      await addMultiValueFixtureColumn(page);
      const dist = await tokenDistribution(page);
      expect(dist.rowCount, `the fixture must cover every one of demog's ${fullRowCount} rows; it `
        + `covers ${dist.rowCount}`).toBe(fullRowCount);
      expect(dist.emptyRows, `${dist.emptyRows} fixture rows carry no token at all, so missing-value `
        + 'semantics would muddy every row count below').toBe(0);
      const countA = dist.counts[MVF_TOKEN_A] ?? 0;
      const countB = dist.counts[MVF_TOKEN_B] ?? 0;
      expect(countA, `token "${MVF_TOKEN_A}" must select a strict subset to be worth asserting; it `
        + `selects ${countA} of ${fullRowCount} rows (distribution ${JSON.stringify(dist.counts)})`)
        .toBeGreaterThan(0);
      expect(countA, `token "${MVF_TOKEN_A}" selects all ${countA} rows, so a filter that did `
        + 'nothing would satisfy the single-token check').toBeLessThan(fullRowCount);
      expect(countB, `token "${MVF_TOKEN_B}" must select a strict subset; it selects ${countB} rows`)
        .toBeGreaterThan(0);
      expect(countB, `token "${MVF_TOKEN_B}" selects every row`).toBeLessThan(fullRowCount);
      expect(dist.pairAnd, `"${MVF_TOKEN_A}" and "${MVF_TOKEN_B}" never co-occur, so an AND that `
        + 'returned nothing would be indistinguishable from a broken filter').toBeGreaterThan(0);
      expect(dist.pairAnd, `the AND set (${dist.pairAnd}) must be strictly smaller than both single `
        + `sets (${MVF_TOKEN_A} ${countA}, ${MVF_TOKEN_B} ${countB}) or the mode check is vacuous`)
        .toBeLessThan(Math.min(countA, countB));
      expect(dist.pairOr, `the OR set (${dist.pairOr}) must be strictly larger than both single sets `
        + `(${MVF_TOKEN_A} ${countA}, ${MVF_TOKEN_B} ${countB})`).toBeGreaterThan(Math.max(countA, countB));
      expect(dist.pairOr, `the OR set (${dist.pairOr}) must still exclude some row of ${fullRowCount}`)
        .toBeLessThan(fullRowCount);

      await v.drivePanelMenuLeaf(page, 'Filters', 'Add Filter', 'Multi Value...');
      await page.locator(MULTI_VALUE_DIALOG).first().waitFor({timeout: 15_000});

      const dialog = await page.evaluate((sel) => {
        const d = document.querySelector(sel)!;
        return {
          column: (d.querySelector('.d4-column-selector-column')?.textContent ?? '').trim(),
          okDisabled: !!d.querySelector('[name="button-OK"]')?.classList.contains('disabled'),
          separatorInvalid: !!d.querySelector('[name="input-Separator"]')?.classList.contains('d4-invalid'),
          balloons: [...document.querySelectorAll('[class*="balloon"]')].map((b) => (b.textContent ?? '').trim()),
        };
      }, MULTI_VALUE_DIALOG);
      expect(dialog.column, 'demog carries no MultiValueSeparator column, so the add must fall back to '
        + 'a categorical one — an empty selector means it resolved no column at all').not.toBe('');
      expect(dialog.balloons.filter((b) => b.includes(NO_CATEGORICAL_BALLOON)),
        `the fallback found "${dialog.column}", so the "${NO_CATEGORICAL_BALLOON}" short-circuit must `
        + `not have fired; balloons on screen: [${dialog.balloons.join(' | ')}]`).toEqual([]);
      expect(dialog.column, `the dialog already opened on "${MVF_COLUMN}", so the column pick below `
        + 'would change nothing and prove nothing').not.toBe(MVF_COLUMN);
      expect(dialog.okDisabled, 'the dialog must refuse to commit before a separator is entered').toBe(true);
      expect(dialog.separatorInvalid, 'the empty separator must be marked invalid').toBe(true);

      await page.evaluate((sel) => {
        const d = document.querySelector(sel)!;
        document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
        (d.querySelector('.d4-column-selector-column') as HTMLElement | null)
          ?.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
      }, MULTI_VALUE_DIALOG);
      await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
        null, {timeout: 10_000, polling: 100});
      await page.keyboard.type(MVF_COLUMN.toLowerCase(), {delay: 40});
      await page.keyboard.press('ArrowDown');
      await page.keyboard.press('Enter');
      await expect.poll(async () => page.evaluate((sel) =>
        (document.querySelector(`${sel} .d4-column-selector-column`)?.textContent ?? '').trim(), MULTI_VALUE_DIALOG),
      {timeout: 15_000, intervals: [200, 400, 800],
        message: `the dialog must now be pointed at "${MVF_COLUMN}"; it opened on "${dialog.column}", and `
          + 'a dialog that silently kept that column would build a card over the wrong values'})
        .toBe(MVF_COLUMN);

      await page.locator(`${MULTI_VALUE_DIALOG} [name="input-Separator"]`).click();
      await page.keyboard.type(MVF_SEPARATOR);
      await expect.poll(async () => page.evaluate((sel) =>
        !document.querySelector(`${sel} [name="button-OK"]`)!.classList.contains('disabled'), MULTI_VALUE_DIALOG),
      {timeout: 10_000, intervals: [200, 400, 800],
        message: `entering the "${MVF_SEPARATOR}" separator must enable OK; it stayed disabled`}).toBe(true);

      await page.locator(`${MULTI_VALUE_DIALOG} [name="button-OK"]`).click();
      await dialogGone(page, MULTI_VALUE_DIALOG);
      await expect.poll(async () => cardCount(page),
        {timeout: 20_000, intervals: [300, 600, 1200],
          message: 'the Multi Value commit must add exactly one card to the empty panel'}).toBe(1);

      const after = await multiValueFilterColumns(page, MULTI_VALUE_FILTER_TYPE);
      expect(after, `the new card must really be a "${MULTI_VALUE_FILTER_TYPE}" filter on `
        + `"${MVF_COLUMN}", not a card that merely took that column; got [${after.join(', ')}]`)
        .toEqual([MVF_COLUMN]);
      expect(await orderedCaptions(page), 'the only caption must be the column the dialog was pointed at')
        .toEqual([MVF_COLUMN]);
      await expect.poll(async () => trueCount(page),
        {timeout: 20_000, intervals: [300, 600, 1200],
          message: 'a multi-value filter with nothing selected must exclude no further row, so the '
            + `count must stay at ${fullRowCount}`}).toBe(fullRowCount);

      const opened = await multiValueCardState(page);
      expect(opened, `the committed "${MVF_COLUMN}" multi-value filter cannot be read back, so no `
        + 'tick below could be measured').not.toBeNull();
      expect(opened!.include, `a freshly committed card must include nothing; it includes `
        + `[${opened!.include.join(', ')}]`).toEqual([]);
      expect(opened!.mode, `the card must open in AND mode for the mode switch below to be a change; `
        + `it opened in "${opened!.mode}"`).toBe('AND');

      const rowTokens: string[] = [];
      for (let row = 0; row < dist.tokens.length; row++) {
        const ticked = await clickMultiValueRow(page, row);
        expect(ticked.include.length, `ticking row ${row} of the "${MVF_COLUMN}" card must select `
          + `exactly one token; it selected [${ticked.include.join(', ')}] (mapped so far `
          + `[${rowTokens.join(', ')}])`).toBe(1);
        rowTokens.push(ticked.include[0]);
        const unticked = await clickMultiValueRow(page, row);
        expect(unticked.include, `unticking row ${row} ("${rowTokens[row]}") must leave the card `
          + `empty again; it left [${unticked.include.join(', ')}]`).toEqual([]);
      }
      expect([...rowTokens].sort(), `the card must render exactly one row per distinct token of `
        + `"${MVF_COLUMN}"; rows 0..${dist.tokens.length - 1} map to [${rowTokens.join(', ')}] while the `
        + `column's own values split on "${MVF_SEPARATOR}" into [${dist.tokens.join(', ')}]`)
        .toEqual(dist.tokens);

      const rowA = rowTokens.indexOf(MVF_TOKEN_A);
      const rowB = rowTokens.indexOf(MVF_TOKEN_B);
      expect(rowA, `the card renders no row for "${MVF_TOKEN_A}"; rows are [${rowTokens.join(', ')}]`)
        .toBeGreaterThanOrEqual(0);
      expect(rowB, `the card renders no row for "${MVF_TOKEN_B}"; rows are [${rowTokens.join(', ')}]`)
        .toBeGreaterThanOrEqual(0);

      const single = await clickMultiValueRow(page, rowA);
      expect(single.include, `ticking row ${rowA} must select exactly ["${MVF_TOKEN_A}"]; it selected `
        + `[${single.include.join(', ')}]`).toEqual([MVF_TOKEN_A]);
      expect(single.trueCount, `with only "${MVF_TOKEN_A}" ticked the grid must show exactly the `
        + `${countA} rows whose own "${MVF_COLUMN}" value contains that token; it shows `
        + `${single.trueCount}`).toBe(countA);

      const pairAnd = await clickMultiValueRow(page, rowB);
      expect([...pairAnd.include].sort(), `ticking row ${rowB} must ADD "${MVF_TOKEN_B}" rather than `
        + `replace "${MVF_TOKEN_A}"; the card now includes [${pairAnd.include.join(', ')}]`)
        .toEqual([MVF_TOKEN_A, MVF_TOKEN_B].sort());
      expect(pairAnd.mode, 'the second tick must not have changed the mode').toBe('AND');
      expect(pairAnd.trueCount, `in AND mode both tokens must be required, so the count must be the `
        + `${dist.pairAnd} rows carrying BOTH "${MVF_TOKEN_A}" and "${MVF_TOKEN_B}"; it is `
        + `${pairAnd.trueCount} (single-token counts ${countA} and ${countB})`).toBe(dist.pairAnd);
      expect(pairAnd.trueCount, `the AND count ${pairAnd.trueCount} must be strictly below the `
        + `single-token count ${countA}, or adding the second token did nothing`).toBeLessThan(single.trueCount);

      const or = await toggleMultiValueMode(page);
      expect(or.text, `clicking the card's mode control must flip its label to OR; it reads `
        + `"${or.text}"`).toBe('OR');
      expect(or.state.mode, `the filter must report OR mode after the click; it reports `
        + `"${or.state.mode}"`).toBe('OR');
      expect([...or.state.include].sort(), 'switching the mode must not change which tokens are ticked')
        .toEqual([MVF_TOKEN_A, MVF_TOKEN_B].sort());
      expect(or.state.trueCount, `in OR mode either token must qualify, so the count must be the `
        + `${dist.pairOr} rows carrying "${MVF_TOKEN_A}" OR "${MVF_TOKEN_B}"; it is `
        + `${or.state.trueCount}`).toBe(dist.pairOr);
      expect(or.state.trueCount, `the OR count ${or.state.trueCount} must be strictly above the AND `
        + `count ${pairAnd.trueCount}, or the mode switch did nothing`).toBeGreaterThan(pairAnd.trueCount);

      const back = await toggleMultiValueMode(page);
      expect(back.text, `clicking the mode control again must flip the label back to AND; it reads `
        + `"${back.text}"`).toBe('AND');
      expect(back.state.trueCount, `switching back to AND must restore the ${dist.pairAnd}-row `
        + `intersection; it shows ${back.state.trueCount}`).toBe(dist.pairAnd);

      const dropB = await clickMultiValueRow(page, rowB);
      expect(dropB.include, `unticking row ${rowB} must leave only ["${MVF_TOKEN_A}"]; it left `
        + `[${dropB.include.join(', ')}]`).toEqual([MVF_TOKEN_A]);
      expect(dropB.trueCount, `back to one token the count must be "${MVF_TOKEN_A}"'s own ${countA} `
        + `rows; it is ${dropB.trueCount}`).toBe(countA);

      const cleared = await clickMultiValueRow(page, rowA);
      expect(cleared.include, `unticking row ${rowA} must leave the card including nothing; it left `
        + `[${cleared.include.join(', ')}]`).toEqual([]);
      expect(cleared.trueCount, `a multi-value card with nothing ticked must exclude no row, so all `
        + `${fullRowCount} demog rows must be back; ${cleared.trueCount} are shown`).toBe(fullRowCount);
    }
    finally {
      await removeMultiValueFixture(page);
    }
    expect(await columnNames(page), `the "${MVF_COLUMN}" fixture this step added must be gone again; `
      + `columns before were [${columnsBefore.join(', ')}]`).toEqual(columnsBefore);
  });

  await softStep('Scenario 6 Step 3 — Filter to Column... writes the filter into a boolean column', async () => {
    const narrowed = await narrowToTopCategory(page, 'SEX');
    expect(narrowed, 'the SEX narrowing must leave some rows, or the column below could be all-false '
      + 'and still match').toBeGreaterThan(0);
    expect(narrowed, 'the SEX narrowing must exclude some rows, or an all-true column would match')
      .toBeLessThan(fullRowCount);

    const columnsBefore = await columnNames(page);
    let created = '';
    try {
      await v.drivePanelMenuLeaf(page, 'Filters', null, 'Filter to Column...');
      await page.locator(FILTER_TO_COLUMN_DIALOG).first().waitFor({timeout: 15_000});
      const proposed = await page.evaluate((sel) =>
        (document.querySelector(`${sel} [name="input-Name"]`) as HTMLInputElement | null)?.value ?? '',
      FILTER_TO_COLUMN_DIALOG);
      expect(proposed, 'the dialog names the new column after the active filter; an empty proposal '
        + 'means it read no filter to save').not.toBe('');
      expect(columnsBefore, `the table already carries a "${proposed}" column, so its presence after `
        + `the commit would prove nothing; columns: [${columnsBefore.join(', ')}]`).not.toContain(proposed);

      await page.locator(`${FILTER_TO_COLUMN_DIALOG} [name="button-OK"]`).click();
      created = proposed;
      await dialogGone(page, FILTER_TO_COLUMN_DIALOG);
      await expect.poll(async () => columnNames(page),
        {timeout: 20_000, intervals: [300, 600, 1200],
          message: `no "${proposed}" column appeared on the table after the commit`}).toContain(proposed);

      const type = await page.evaluate((c) => grok.shell.tv.dataFrame.col(c)?.type ?? '', proposed);
      expect(type, `"${proposed}" must be a boolean column; got type "${type}"`).toBe('bool');

      const trues = await boolColumnTrueCount(page, proposed);
      expect(trues, `"${proposed}" could not be counted (-1 means it is missing or not boolean), so `
        + 'the comparison below would be meaningless').toBeGreaterThanOrEqual(0);
      expect(trues, `"${proposed}" must be true on exactly the ${narrowed} rows the panel was showing `
        + `when the command ran; the column itself carries ${trues} true values`).toBe(narrowed);
      await holdTrueCount(page, narrowed, 'saving the filter into a column must not change what the '
        + 'panel filters, and it must not drift over the seconds after the commit either');
    }
    finally {
      if (created !== '')
        await page.evaluate((c) => {
          const df = grok.shell.tv.dataFrame;
          if (df.col(c)) df.columns.remove(c);
        }, created);
    }
    expect(await columnNames(page), 'the column this step added must be gone again; columns before '
      + `were [${columnsBefore.join(', ')}]`).toEqual(columnsBefore);
  });

  v.finishSpec();
});
