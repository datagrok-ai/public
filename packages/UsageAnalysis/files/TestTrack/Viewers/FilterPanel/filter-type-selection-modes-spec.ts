/* ---
realizes: [filters.cp.filter-type-and-selection-modes]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {addCardViaColumnSelector, cardCount, expectHeaderCounterNow, headerCounterTarget,
  trueCount} from '../../helpers/filter-panel';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const fullRowCount = 5850;

async function orderedCaptions(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name')]
      .map((e) => (e.textContent ?? '').trim()));
}

async function cardFilterType(page: Page, column: string): Promise<string | null> {
  return page.evaluate((col) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => ((c.querySelector('.d4-filter-column-name') as HTMLElement | null)?.textContent ?? '').trim() === col);
    if (!card) return null;
    if (card.querySelector('[name="viewer-Grid"]')) return 'categorical';
    if (card.querySelector('[name="viewer-Histogram"]')) return 'histogram';
    return null;
  }, column);
}

async function cardBody(page: Page, column: string): Promise<{grid: boolean; histo: boolean; menuIcon: boolean}> {
  return page.evaluate((col) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col);
    return {
      grid: !!card?.querySelector('[name="viewer-Grid"]'),
      histo: !!card?.querySelector('[name="viewer-Histogram"]'),
      menuIcon: !!card?.querySelector('[name="icon-font-icon-menu"]'),
    };
  }, column);
}

async function cardCanvasBox(page: Page, column: string): Promise<{w: number; h: number} | null> {
  return page.evaluate((col) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col);
    const canvas = card?.querySelector('canvas') as HTMLCanvasElement | null;
    if (!canvas) return null;
    const r = canvas.getBoundingClientRect();
    return {w: r.width, h: r.height};
  }, column);
}

async function clickCardTypeSwitch(page: Page, column: string, icon: 'icon-list' | 'icon-signal'): Promise<void> {
  await page.evaluate(({col, ic}) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col);
    if (!card) throw new Error(`the panel carries no ${col} card to switch the filter type on`);
    const target = card.querySelector(`[name="${ic}"]`) as HTMLElement | null;
    if (!target) {
      const offered = [...card.querySelectorAll('[name^="icon-"]')].map((e) => e.getAttribute('name')).join(' | ');
      throw new Error(`the ${col} card header carries no [name="${ic}"] type-switch icon; icons present: ${offered}`);
    }
    target.click();
  }, {col: column, ic: icon});
  await page.waitForTimeout(1200);
}

async function cardIndicatorMenu(page: Page, column: string,
  opts: {hover?: string[]; path?: string[]} = {}): Promise<string[]> {
  return page.evaluate(async ({col, hover, path}) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col);
    if (!card) throw new Error(`the panel carries no ${col} card to open an indicator menu on`);
    const indicator = card.querySelector('.d4-filter-indicator') as HTMLElement | null;
    if (!indicator) throw new Error(`the ${col} card carries no .d4-filter-indicator to open its menu with`);
    indicator.click();

    const norm = (s: string | null | undefined) => (s ?? '').trim().toLowerCase();
    const items = () => {
      const popup = [...document.querySelectorAll('.d4-menu-popup')].pop();
      return [...(popup?.querySelectorAll('.d4-menu-item-label') ?? [])];
    };
    const texts = () => items().map((i) => (i.textContent ?? '').trim());
    const find = (t: string) =>
      items().find((i) => norm(i.textContent) === norm(t))?.closest('.d4-menu-item') ?? null;
    const waitForItem = async (t: string) => {
      const deadline = Date.now() + 5000;
      let it = find(t);
      while (!it && Date.now() < deadline) {
        await new Promise((r) => setTimeout(r, 100));
        it = find(t);
      }
      return it;
    };
    const openDeadline = Date.now() + 5000;
    while (items().length === 0 && Date.now() < openDeadline)
      await new Promise((r) => setTimeout(r, 100));
    if (items().length === 0)
      throw new Error(`the ${col} card indicator click opened no menu popup`);

    const seen = texts();
    const merge = () => { for (const t of texts()) if (!seen.includes(t)) seen.push(t); };
    const hoverItem = async (item: Element) => {
      const b = item.getBoundingClientRect();
      for (const type of ['mouseover', 'mousemove'])
        item.dispatchEvent(new MouseEvent(type, {bubbles: true, clientX: b.x + 5, clientY: b.y + 5}));
      await new Promise((r) => setTimeout(r, 700));
    };

    for (const g of hover ?? []) {
      const item = await waitForItem(g);
      if (!item)
        throw new Error(`${col} indicator menu: group "${g}" not found; visible: ${texts().join(' | ')}`);
      await hoverItem(item);
      merge();
    }

    for (let i = 0; i < (path ?? []).length; i++) {
      const label = path![i];
      const item = await waitForItem(label);
      if (!item)
        throw new Error(`${col} indicator menu: "${label}" not found; visible: ${texts().join(' | ')}`);
      if (i < path!.length - 1) {
        await hoverItem(item);
        merge();
      }
      else {
        (item as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 800));
      }
    }

    document.body.click();
    await new Promise((r) => setTimeout(r, 300));
    return seen;
  }, {col: column, hover: opts.hover ?? [], path: opts.path ?? []});
}

async function setCardMode(page: Page, column: string, mode: 'Radio' | 'Multi-Select'): Promise<void> {
  await cardIndicatorMenu(page, column, {path: ['Mode', mode]});
}

async function runBatchOp(page: Page, column: string, label: string): Promise<number> {
  await cardIndicatorMenu(page, column, {path: [label]});
  return trueCount(page);
}

const HEADER_BAND_H = 10;
const ROW_PITCH = 27;
const ROW_CENTRE_OFFSET = 13;
const X_CHECKBOX = 10;
const X_NAME = 60;

const rowCentreY = (rowIndex: number): number => HEADER_BAND_H + ROW_PITCH * rowIndex + ROW_CENTRE_OFFSET;

type CardState = {selected: string[] | null; count: number};
type RowClick = {rowIndex: number; y: number; before: CardState; after: CardState};

const fmt = (s: CardState): string => `selected=[${(s.selected ?? []).join(',')}] trueCount=${s.count}`;

async function clickCategoryRow(page: Page, column: string, rowIndex: number,
  target: 'checkbox' | 'name'): Promise<RowClick> {
  const y = rowCentreY(rowIndex);
  const res = await page.evaluate(async ({col, cx, cy}) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col);
    const overlay = card?.querySelector('[name="viewer-Grid"] [name="overlay"]') as HTMLElement | null;
    if (!overlay) return null;
    const read = () => {
      const st = grok.shell.tv.getFiltersGroup().getStates(col, 'categorical');
      const sel = st?.[0]?.selected;
      return {
        selected: (Array.isArray(sel) ? sel.map((s: any) => String(s)) : null) as string[] | null,
        count: grok.shell.tv.dataFrame.filter.trueCount as number,
      };
    };
    const key = (s: {selected: string[] | null; count: number}) =>
      `${s.count}|${s.selected === null ? '<null>' : s.selected.join(',')}`;
    const before = read();
    const rect = overlay.getBoundingClientRect();
    const o = {
      bubbles: true, cancelable: true, view: window, button: 0,
      clientX: rect.left + cx,
      clientY: rect.top + cy,
    };
    overlay.dispatchEvent(new MouseEvent('mousedown', o));
    overlay.dispatchEvent(new MouseEvent('mouseup', o));
    overlay.dispatchEvent(new MouseEvent('click', o));
    let after = before;
    let prevKey = key(before);
    for (let waited = 0; waited < 4000; waited += 100) {
      await new Promise((r) => setTimeout(r, 100));
      after = read();
      const k = key(after);
      if (k !== key(before) && k === prevKey) break;
      prevKey = k;
    }
    return {before, after};
  }, {col: column, cx: target === 'checkbox' ? X_CHECKBOX : X_NAME, cy: y});
  expect(res, `the ${column} card must expose a categorical [name="overlay"] body to click`).not.toBeNull();
  return {rowIndex, y, before: res!.before, after: res!.after};
}

async function cardRenderedCategories(page: Page, column: string, maxRows: number): Promise<string[]> {
  const rendered: string[] = [];
  for (let row = 0; row < maxRows; row++) {
    const click = await clickCategoryRow(page, column, row, 'name');
    const after = click.after.selected;
    const before = click.before.selected;
    const landed = after !== null && after.length === 1 &&
      (before === null || before.length !== 1 || before[0] !== after[0]);
    if (landed && !rendered.includes(after![0]))
      rendered.push(after![0]);
  }
  return rendered;
}

type RowPick = {row: number; y: number; label: string; count: number; click: RowClick};

async function pickCategoryRow(page: Page, column: string, rows: number[], counts: Record<string, number>,
  accept: (label: string, count: number) => boolean, probes: string[]): Promise<RowPick | null> {
  for (const row of rows) {
    const click = await clickCategoryRow(page, column, row, 'name');
    const label = click.after.selected?.length === 1 ? click.after.selected[0] : '';
    const count = counts[label] ?? 0;
    probes.push(`row ${row} (y=${click.y}) -> ${label === '' ? fmt(click.after) : `'${label}' (${count} rows)`}`);
    if (label !== '' && count > 0 && accept(label, count))
      return {row, y: click.y, label, count, click};
  }
  return null;
}

async function cardSelectedCategories(page: Page, column: string): Promise<string[] | null> {
  return page.evaluate((col) => {
    const st = grok.shell.tv.getFiltersGroup().getStates(col, 'categorical');
    const sel = st?.[0]?.selected;
    return Array.isArray(sel) ? sel.map((s: any) => String(s)) : null;
  }, column);
}

async function heldTrueCount(page: Page, expected: number, why: string, holdMs = 2500): Promise<void> {
  const deadline = Date.now() + holdMs;
  for (;;) {
    expect(await trueCount(page), why).toBe(expected);
    if (Date.now() >= deadline) return;
    await page.waitForTimeout(250);
  }
}

async function inCardSearch(page: Page, column: string, fragment: string): Promise<number> {
  return page.evaluate(async ({col, frag}) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col) as HTMLElement;
    (card.querySelector('[name="icon-search"]') as HTMLElement | null)?.click();
    await new Promise((r) => setTimeout(r, 500));
    const input = card.querySelector('input.d4-search-input[placeholder="Search"]') as HTMLInputElement;
    input.focus();
    const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    setter.call(input, frag);
    input.dispatchEvent(new Event('input', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 1200));
    return grok.shell.tv.dataFrame.filter.trueCount;
  }, {col: column, frag: fragment});
}

async function clearInCardSearch(page: Page, column: string): Promise<void> {
  await page.evaluate(async (col) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col) as HTMLElement;
    const input = card.querySelector('input.d4-search-input[placeholder="Search"]') as HTMLInputElement;
    if (input) {
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
      setter.call(input, '');
      input.dispatchEvent(new Event('input', {bubbles: true}));
    }
    await new Promise((r) => setTimeout(r, 800));
  }, column);
}

async function pasteInCardSearch(page: Page, column: string, values: string[], trailingSep: boolean): Promise<number> {
  return page.evaluate(async ({col, vals, trail}) => {
    const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col) as HTMLElement;
    (card.querySelector('[name="icon-search"]') as HTMLElement | null)?.click();
    await new Promise((r) => setTimeout(r, 500));
    const input = card.querySelector('input.d4-search-input[placeholder="Search"]') as HTMLInputElement;
    const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    setter.call(input, '');
    input.dispatchEvent(new Event('input', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 600));
    input.focus();
    const text = vals.join('\n') + (trail ? '\n' : '');
    const dt = new DataTransfer();
    dt.setData('text', text);
    input.dispatchEvent(new ClipboardEvent('paste', {bubbles: true, cancelable: true, clipboardData: dt}));
    await new Promise((r) => setTimeout(r, 1500));
    return grok.shell.tv.dataFrame.filter.trueCount;
  }, {col: column, vals: values, trail: trailingSep});
}

test('Filter Panel — Filter type switching and category selection modes', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  const consoleErrors: string[] = [];
  page.on('console', (msg) => { if (msg.type() === 'error') consoleErrors.push(msg.text()); });
  const AMBIENT = 'Permissions policy violation: compute-pressure';
  const errorSet = () => new Set(consoleErrors.filter((t) => !t.includes(AMBIENT)));

  await softStep('Setup: empty the panel and add an AGE histogram card', async () => {
    await v.drivePanelMenuLeaf(page, 'Filters', null, 'Remove All');
    await expect.poll(async () => cardCount(page), {timeout: 20_000, intervals: [400, 800, 1500]}).toBe(0);
    expect(await trueCount(page), 'removing every card leaves the table unfiltered').toBe(fullRowCount);
    await addCardViaColumnSelector(page, 'AGE');
    expect(await orderedCaptions(page)).toEqual(['AGE']);
    expect(await cardFilterType(page, 'AGE'), 'AGE is numerical, so its default card is a histogram')
      .toBe('histogram');
    expect(await trueCount(page)).toBe(fullRowCount);
    await expectHeaderCounterNow(page, null, 'an untouched card filters nothing, so the counter is off screen at Setup');
  });

  await softStep('Step 1 — switch AGE to categorical; menu icon present, count unchanged', async () => {
    await clickCardTypeSwitch(page, 'AGE', 'icon-list');
    const body = await cardBody(page, 'AGE');
    expect(await cardFilterType(page, 'AGE'), 'card is categorical after the switch').toBe('categorical');
    expect(body.grid, 'categorical card body is a viewer-Grid, not a histogram').toBe(true);
    expect(body.histo).toBe(false);
    expect(body.menuIcon, 'GROK-19915: filter menu icon present on rebuilt categorical header').toBe(true);
    expect(await trueCount(page), 'switch alone does not filter').toBe(fullRowCount);
    await expectHeaderCounterNow(page, null, 'no category is checked yet, so the counter stays off screen');
  });

  let catA = '';
  let catB = '';
  let countA = 0;
  let countB = 0;
  let rowA = -1;
  let rowB = -1;
  let trueCountTwoCats = 0;

  await softStep('Step 2 — check two categories by clicking their rows; count drops, counter reads 1', async () => {
    const ageCounts: Record<string, number> = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('AGE');
      const n: Record<string, number> = {};
      for (let i = 0; i < df.rowCount; i++) {
        const label = String(col.get(i));
        n[label] = (n[label] ?? 0) + 1;
      }
      return n;
    });

    const probesA: string[] = [];
    const pickA = await pickCategoryRow(page, 'AGE', [0, 1, 2, 3, 4, 5], ageCounts, () => true, probesA);
    expect(pickA, `no row in the scanned range yielded a category with rows — a blank category or a click that never reached the canvas leaves every count comparison unfalsifiable. Probes: ${probesA.join('; ')}`)
      .not.toBeNull();
    rowA = pickA!.row;
    const first = pickA!.click;
    expect(first.after.selected,
      `the card must expose its checked-category list after the click; pre=${fmt(first.before)} post=${fmt(first.after)}`)
      .not.toBeNull();
    expect(first.after.selected!.length,
      `clicking a category name checks exactly that one category; pre=${fmt(first.before)} post=${fmt(first.after)}`)
      .toBe(1);
    catA = pickA!.label;
    countA = pickA!.count;
    expect(countA, `the clicked category '${catA}' must carry rows, or the count comparisons prove nothing`)
      .toBeGreaterThan(0);
    expect(first.after.count,
      `the exclusive select on '${catA}' leaves its own count ${countA}; pre=${fmt(first.before)} post=${fmt(first.after)}`)
      .toBe(countA);
    expect(first.after.count,
      `the exclusive select must narrow the table below ${fullRowCount}; post=${fmt(first.after)}`)
      .toBeLessThan(fullRowCount);

    const probesB: string[] = [];
    const rowsB = [1, 2, 3, 4, 5].map((offset) => rowA + offset);
    const pickB = await pickCategoryRow(page, 'AGE', rowsB, ageCounts,
      (label, count) => label !== catA && count !== countA, probesB);
    expect(pickB, `a click on a row below row ${rowA} must reach a second category whose own count differs from '${catA}' (${countA}) — a selection that never moves means the gesture did not land on a category row. Probes: ${probesB.join('; ')}`)
      .not.toBeNull();
    catB = pickB!.label;
    countB = pickB!.count;
    rowB = pickB!.row;

    const both = await clickCategoryRow(page, 'AGE', rowA, 'checkbox');
    expect(both.after.selected,
      `the card must expose its checked-category list; pre=${fmt(both.before)} post=${fmt(both.after)}`)
      .not.toBeNull();
    expect(both.after.count,
      `the checkbox click must not leave the card unfiltered at ${fullRowCount} — that is the invert() flip, not an added category; pre=${fmt(both.before)} post=${fmt(both.after)}`)
      .not.toBe(fullRowCount);
    expect(both.after.selected!.length,
      `the checkbox click ADDS in Multi-Select — '${catA}' and '${catB}' both checked, not one; pre=${fmt(both.before)} post=${fmt(both.after)}`)
      .toBe(2);
    expect([...both.after.selected!].sort(),
      `the two checked categories are exactly the two that were clicked ('${catA}', '${catB}'); post=${fmt(both.after)}`)
      .toEqual([catA, catB].sort());
    trueCountTwoCats = both.after.count;
    expect(trueCountTwoCats,
      `two checked categories give the sum ${countA} + ${countB} = ${countA + countB}; post=${fmt(both.after)}`)
      .toBe(countA + countB);
    expect(trueCountTwoCats, `two checked categories must stay below the full ${fullRowCount}`).toBeLessThan(fullRowCount);
    await expectHeaderCounterNow(page, '1', 'one card is filtering, so the active-filter counter is on screen reading 1');
  });

  await softStep('Step 3 — Radio mode collapses to one category; batch items absent', async () => {
    await setCardMode(page, 'AGE', 'Radio');
    const labels = await cardIndicatorMenu(page, 'AGE', {hover: ['Mode']});
    expect(labels, 'Select all absent in Radio mode').not.toContain('Select all');
    expect(labels, 'Deselect all absent in Radio mode').not.toContain('Deselect all');
    expect(labels, 'Invert all absent in Radio mode').not.toContain('Invert all');
    expect(labels).toContain('Mode');
    expect(labels).toContain('Radio');
    const collapsed = await trueCount(page);
    expect(collapsed, 'Radio mode keeps only the first checked category (its own count)').toBe(countA);
    expect(collapsed).toBeLessThan(trueCountTwoCats);
    const selected = await cardSelectedCategories(page, 'AGE');
    expect(selected, 'the card must expose its checked-category list').not.toBeNull();
    expect(selected!.length, 'Radio mode leaves exactly one category checked').toBe(1);
    expect(selected![0], 'the surviving category is the first of the two checked in Step 2').toBe(catA);
  });

  await softStep('Step 4 — clicking a second category in Radio replaces, not unions', async () => {
    expect(rowB, `Step 2 must have located the second category row before the Radio click can address it; rowA=${rowA} rowB=${rowB}`)
      .toBeGreaterThan(0);
    const replaced = await clickCategoryRow(page, 'AGE', rowB, 'checkbox');

    expect(replaced.before.selected, `Radio mode must enter this step with a checked category; pre=${fmt(replaced.before)}`)
      .not.toBeNull();
    expect(replaced.before.selected!, `the Radio pre-state holds '${catA}' alone; pre=${fmt(replaced.before)}`).toEqual([catA]);
    expect(replaced.before.count,
      `the Radio pre-state must already be narrowing to '${catA}' (${countA} rows), not sitting at the full ${fullRowCount}; pre=${fmt(replaced.before)}`)
      .toBe(countA);

    expect(replaced.after.selected, `the card must expose its checked-category list; post=${fmt(replaced.after)}`).not.toBeNull();
    expect(replaced.after.selected!.join(','),
      `the Radio click on row ${rowB} (y=${replaced.y}) must move the checked category off '${catA}'; pre=${fmt(replaced.before)} post=${fmt(replaced.after)}`)
      .not.toBe(replaced.before.selected!.join(','));
    expect(replaced.after.count,
      `the Radio click must move the row count off the pre-click ${replaced.before.count}; pre=${fmt(replaced.before)} post=${fmt(replaced.after)}`)
      .not.toBe(replaced.before.count);

    expect(replaced.after.selected!.length,
      `'${catB}' replaced '${catA}' — one category checked, not two; post=${fmt(replaced.after)}`).toBe(1);
    expect(replaced.after.selected![0],
      `the checked category is the clicked one, '${catB}'; post=${fmt(replaced.after)}`).toBe(catB);
    expect(replaced.after.count,
      `Radio replacement leaves '${catB}' own count ${countB}; a count of ${fullRowCount} here with the selection moved is the unexplained Radio anomaly, not a pass; pre=${fmt(replaced.before)} post=${fmt(replaced.after)}`)
      .toBe(countB);
    expect(replaced.after.count,
      `Radio replaces rather than unions — not ${countA} + ${countB} = ${countA + countB}; post=${fmt(replaced.after)}`)
      .not.toBe(countA + countB);
    expect(replaced.after.count,
      `the card must still be filtering, not flipped back to the full ${fullRowCount}; post=${fmt(replaced.after)}`)
      .not.toBe(fullRowCount);
    await expectHeaderCounterNow(page, '1', 'the card is still filtering, so the counter stays on screen reading 1');
  });

  await softStep('Step 5 — back to Multi-Select; count unchanged, batch items return', async () => {
    await setCardMode(page, 'AGE', 'Multi-Select');
    await heldTrueCount(page, countB,
      `the mode switch alone must not move the row set off '${catB}' (${countB} rows), and must not move it later either`);
    await expectHeaderCounterNow(page, '1', 'the card still filters after the mode switch');
    const labels = await cardIndicatorMenu(page, 'AGE', {hover: ['Mode']});
    expect(labels, 'Select all returns in Multi-Select mode').toContain('Select all');
    expect(labels, 'Deselect all returns in Multi-Select mode').toContain('Deselect all');
    expect(labels, 'Invert all returns in Multi-Select mode').toContain('Invert all');
    expect(await orderedCaptions(page)).toContain('AGE');
  });

  await softStep('Step 6 — in-card search narrows the main dataframe', async () => {
    const preSearch = await runBatchOp(page, 'AGE', 'Select all');
    expect(preSearch).toBe(fullRowCount);
    const afterSearch = await inCardSearch(page, 'AGE', '5');
    expect(afterSearch, 'in-card search moves the main dataframe').toBeLessThan(preSearch);
    expect(afterSearch).toBeGreaterThan(0);

    const matchContent = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('AGE');
      const matching = new Set<string>();
      const passing = new Set<string>();
      let passingWithout = 0; let excludedWith = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const label = String(col.get(i));
        const has = label.includes('5');
        if (has) matching.add(label);
        if (df.filter.get(i)) {
          passing.add(label);
          if (!has) passingWithout++;
        }
        else if (has) excludedWith++;
      }
      return {
        matched: matching.size,
        passing: passing.size,
        total: col.categories.length,
        passingWithout,
        excludedWith,
      };
    });
    expect(matchContent.matched, 'the fragment must match at least one AGE category').toBeGreaterThan(0);
    expect(matchContent.matched, 'the fragment must not match every AGE category').toBeLessThan(matchContent.total);
    expect(matchContent.passing, 'the matched-category count in the card equals the categories carrying the fragment')
      .toBe(matchContent.matched);
    expect(matchContent.passingWithout, 'a row whose AGE label lacks the fragment passed the search filter').toBe(0);
    expect(matchContent.excludedWith, 'a row whose AGE label carries the fragment was excluded').toBe(0);
  });

  const pasteVals = ['21', '46', '63'];

  await softStep('Step 7 — multi-line paste filters to the union of pasted values', async () => {
    await clearInCardSearch(page, 'AGE');
    await runBatchOp(page, 'AGE', 'Select all');
    expect(await trueCount(page), 'pre-paste state is the full row set').toBe(fullRowCount);

    const counts = await page.evaluate((vals) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('AGE');
      const n: Record<string, number> = {};
      for (const v of vals) n[v] = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const label = String(col.get(i));
        if (n[label] !== undefined) n[label]++;
      }
      return n;
    }, pasteVals);
    for (const v of pasteVals)
      expect(counts[v], `pasted value ${v} must be a category with rows, or the union proves nothing`)
        .toBeGreaterThan(0);
    const pasteUnion = pasteVals.reduce((sum, v) => sum + counts[v], 0);
    const largestSingle = Math.max(...pasteVals.map((v) => counts[v]));
    expect(pasteUnion, 'the union must be a real subset of the table, or the check is unfalsifiable')
      .toBeLessThan(fullRowCount);

    const withTrailing = await pasteInCardSearch(page, 'AGE', pasteVals, true);
    expect(withTrailing, 'paste filters to the union of the three pasted values').toBe(pasteUnion);
    expect(withTrailing, 'union is not the full row count').not.toBe(fullRowCount);
    expect(withTrailing, 'union is strictly greater than the largest single pasted value').toBeGreaterThan(largestSingle);
    expect(withTrailing, 'union is strictly less than the full row count').toBeLessThan(fullRowCount);
    await expectHeaderCounterNow(page, '1', 'the card filters after the paste, so the counter reads 1');

    const withoutTrailing = await pasteInCardSearch(page, 'AGE', pasteVals, false);
    expect(withoutTrailing, 'trailing separator changes nothing — same union count').toBe(withTrailing);
    expect(withoutTrailing).not.toBe(0);
    expect(withoutTrailing).not.toBe(fullRowCount);
  });

  await softStep('Step 8 — batch operations produce their documented row sets', async () => {
    await clearInCardSearch(page, 'AGE');
    const selectAll = await runBatchOp(page, 'AGE', 'Select all');
    const counterSelectAll = await headerCounterTarget(page);
    const deselectAll = await runBatchOp(page, 'AGE', 'Deselect all');
    const counterDeselectAll = await headerCounterTarget(page);
    const invert = await runBatchOp(page, 'AGE', 'Invert all');
    const counterInvert = await headerCounterTarget(page);
    expect(selectAll, 'Select all keeps every row').toBe(fullRowCount);
    expect(deselectAll, 'Deselect all excludes every row').toBe(0);
    expect(invert, 'Invert all is the complement of Deselect all').toBe(fullRowCount - deselectAll);
    expect(invert, 'inverting the empty selection restores the full row set').toBe(selectAll);
    expect(counterSelectAll.present,
      'the header counter element left the DOM after Select all, so its "off screen" reading below ' +
      'would be satisfied by a counter that stopped rendering entirely').toBe(true);
    expect(counterDeselectAll.present,
      'the header counter element left the DOM after Deselect all, so it reads neither a value nor ' +
      '"hidden"').toBe(true);
    expect(counterInvert.present,
      'the header counter element left the DOM after Invert all, so its "off screen" reading below ' +
      'would be satisfied by a counter that stopped rendering entirely').toBe(true);
    expect(counterSelectAll.visible,
      `no category is excluded after Select all, so the counter is off screen; it read ${JSON.stringify(counterSelectAll.text)}`)
      .toBe(false);
    expect(counterDeselectAll.visible,
      'the card filters after Deselect all, so the counter is on screen').toBe(true);
    expect(counterDeselectAll.text, 'the counter reads 1 after Deselect all').toBe('1');
    expect(counterInvert.visible,
      `Invert all restores the full set, so the counter is off screen again; it read ${JSON.stringify(counterInvert.text)}`)
      .toBe(false);
  });

  await softStep('Step 9 — switch back to histogram reverts the categorical selection', async () => {
    await runBatchOp(page, 'AGE', 'Select all');
    await clickCardTypeSwitch(page, 'AGE', 'icon-signal');
    expect(await cardFilterType(page, 'AGE'), 'card is a histogram again after the revert').toBe('histogram');
    const body = await cardBody(page, 'AGE');
    expect(body.histo, 'histogram card body present after revert').toBe(true);
    expect(body.grid).toBe(false);
    expect(await trueCount(page), 'round-trip revert returns to the full row count').toBe(fullRowCount);
    await expectHeaderCounterNow(page, null, 'nothing filters after the revert, so the counter is off screen');
  });

  let raceCatCountBefore = 0;
  let raceProbeRows = 0;

  await softStep('Scenario 2 setup — add RACE, Select all, record category count', async () => {
    const captionsBefore = await orderedCaptions(page);
    await addCardViaColumnSelector(page, 'RACE');
    const captionsAfter = await orderedCaptions(page);
    expect(captionsAfter, 'the RACE card joined the panel').toContain('RACE');
    expect(captionsAfter.length, 'adding RACE grows the card count by exactly one — nothing else was rebuilt')
      .toBe(captionsBefore.length + 1);
    expect(await cardFilterType(page, 'RACE'), 'RACE is a string column, so its default card is categorical')
      .toBe('categorical');
    await runBatchOp(page, 'RACE', 'Select all');
    raceCatCountBefore = await page.evaluate(() => grok.shell.tv.dataFrame.col('RACE').categories.length);
    expect(raceCatCountBefore, 'demog RACE has four categories').toBe(4);
    raceProbeRows = raceCatCountBefore + 1;
    expect(await trueCount(page)).toBe(fullRowCount);

    const renderedBefore = await cardRenderedCategories(page, 'RACE', raceProbeRows);
    expect(renderedBefore.sort(),
      `the RACE card's own rendered rows, read by clicking each row and reading back the category the card reports checked, must be exactly the four RACE categories before anything is deleted — otherwise the same reading after the undo would prove nothing; got [${renderedBefore.join(',')}]`)
      .toEqual(await page.evaluate(() => [...grok.shell.tv.dataFrame.col('RACE').categories].map(String).sort()));
    await runBatchOp(page, 'RACE', 'Select all');
    expect(await trueCount(page), 'the row-probe leaves the table unfiltered again').toBe(fullRowCount);
  });

  await softStep('Step 10 — delete Asian rows, undo, RACE card re-lists Asian', async () => {
    const del = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.init((i: number) => df.col('RACE').get(i) === 'Asian');
      await new Promise((r) => setTimeout(r, 300));
      const selected = df.selection.trueCount;
      const rowsBefore = df.rowCount;
      const fns = DG.Func.find({name: 'CmdRemoveSelectedRows'});
      await fns[0].prepare({}).call();
      await new Promise((r) => setTimeout(r, 1200));
      const dfa = grok.shell.tv.dataFrame;
      return {selected, rowsBefore, rowsAfter: dfa.rowCount, catsAfter: dfa.col('RACE').categories};
    });
    expect(del.selected, 'Asian rows selected before delete').toBeGreaterThan(0);
    expect(del.rowsAfter, 'deleted rows are removed from the table').toBe(del.rowsBefore - del.selected);
    expect(del.catsAfter, 'the RACE column itself no longer carries the deleted Asian category').not.toContain('Asian');

    const renderedAfterDelete = await cardRenderedCategories(page, 'RACE', raceProbeRows);
    expect(renderedAfterDelete,
      `the RACE card's own rendered rows must drop Asian once its rows are gone; card rows: [${renderedAfterDelete.join(',')}]`)
      .not.toContain('Asian');
    await runBatchOp(page, 'RACE', 'Select all');

    const gridPoint = await page.evaluate(() => {
      const mainGrid = [...document.querySelectorAll('[name="viewer-Grid"]')].find((g) => !g.closest('.d4-filter'));
      const overlay = mainGrid!.querySelector('[name="overlay"]') as HTMLElement;
      const r = overlay.getBoundingClientRect();
      return {x: Math.round(r.left + r.width / 2), y: Math.round(r.top + 40)};
    });
    await page.mouse.click(gridPoint.x, gridPoint.y);
    await page.keyboard.press('Control+z');
    await page.waitForTimeout(1500);

    const restored = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowCount: df.rowCount, cats: [...df.col('RACE').categories].map(String)};
    });
    expect(restored.rowCount, 'Ctrl+Z restores the deleted rows').toBe(fullRowCount);
    expect(restored.cats, 'the RACE column carries Asian again after the undo').toContain('Asian');

    const renderedAfterUndo = await cardRenderedCategories(page, 'RACE', raceProbeRows);
    expect(renderedAfterUndo,
      `GROK-19537: the RACE card's OWN rendered rows must list Asian again after Ctrl+Z. This reads the card, not the column: each card row is clicked and the category the card reports checked is read back, so a card that kept its stale three-row list — the exact shape of the defect — cannot reach Asian and fails here. Card rows: [${renderedAfterUndo.join(',')}]`)
      .toContain('Asian');
    expect(renderedAfterUndo.length,
      `the card's rendered row count returns to the pre-deletion category count; card rows: [${renderedAfterUndo.join(',')}]`)
      .toBe(raceCatCountBefore);

    await runBatchOp(page, 'RACE', 'Select all');
    expect(await trueCount(page), 'the restored table is unfiltered again').toBe(fullRowCount);
  });

  let countBeforeProbeAdded = 0;

  await softStep('Step 11 — zero-variance column filter renders and clicks without console errors', async () => {
    const errorsBefore = errorSet();
    await v.drivePanelMenuLeaf(page, 'Filters', null, 'Remove All');
    await expect.poll(async () => cardCount(page), {timeout: 20_000, intervals: [400, 800, 1500]}).toBe(0);
    countBeforeProbeAdded = await trueCount(page);
    expect(countBeforeProbeAdded, 'the table is unfiltered before the probe column is added').toBe(fullRowCount);

    await page.evaluate(async () => {
      await grok.shell.tv.dataFrame.columns.addNewCalculated('probe_constant', '1');
      await new Promise((r) => setTimeout(r, 800));
    });
    await addCardViaColumnSelector(page, 'probe_constant');
    expect(await orderedCaptions(page)).toContain('probe_constant');

    const canvas = await cardCanvasBox(page, 'probe_constant');
    expect(canvas,
      'GROK-19897: the zero-variance card must paint a canvas body — a card whose body never rendered cannot be clicked, and a step that skipped its own click would report an empty console delta for free')
      .not.toBeNull();
    expect(canvas!.w, 'the zero-variance card canvas must have a non-zero width to click').toBeGreaterThan(0);
    expect(canvas!.h, 'the zero-variance card canvas must have a non-zero height to click').toBeGreaterThan(0);

    const clicked = await page.evaluate(async () => {
      const card = [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
        .find((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === 'probe_constant');
      const canvas = card?.querySelector('canvas') as HTMLCanvasElement | null;
      if (!canvas) throw new Error('the probe_constant card lost its canvas between the presence check and the click');
      const r = canvas.getBoundingClientRect();
      const o = {bubbles: true, cancelable: true, clientX: r.left + r.width / 2, clientY: r.top + r.height / 2};
      for (const type of ['mousedown', 'mouseup', 'click'])
        canvas.dispatchEvent(new MouseEvent(type, o));
      await new Promise((r2) => setTimeout(r2, 1500));
      return true;
    });
    expect(clicked, 'the zero-variance card body was actually clicked').toBe(true);

    const newErrors = [...errorSet()].filter((t) => !errorsBefore.has(t));
    expect(newErrors, `GROK-19897: no new console error on the zero-variance column; got: ${newErrors.join(' | ')}`).toEqual([]);
    const tc = await trueCount(page);
    expect([0, countBeforeProbeAdded],
      `every probe_constant row carries the same value, so a click on its card either keeps all ${countBeforeProbeAdded} rows or excludes them all — ${tc} is neither, which means the card responded with a criterion its own data cannot produce`)
      .toContain(tc);
  });

  await softStep('Step 12 — remove the probe column; its card leaves the panel', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.dataFrame.columns.remove('probe_constant');
      await new Promise((r) => setTimeout(r, 1000));
    });
    expect(await page.evaluate(() => !grok.shell.tv.dataFrame.col('probe_constant')),
      'probe_constant column removed from the table').toBe(true);
    expect(await orderedCaptions(page), 'probe_constant card absent from the panel').not.toContain('probe_constant');
    await heldTrueCount(page, countBeforeProbeAdded,
      `removing the probe column returns filtering to the ${countBeforeProbeAdded} rows that were passing before the probe column was ever added, and holds there`);
  });

  await softStep('Teardown', async () => {
    await page.evaluate(async () => {
      try { grok.shell.tv?.dataFrame?.columns?.remove('probe_constant'); } catch (_) {}
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
    });
  });

  v.finishSpec();
});
