import {expect, Page} from '@playwright/test';
import * as v from '../../helpers/viewers';
import {cardCount, trueCount} from '../../helpers/filter-panel';

declare const grok: any;

export const FULL = 5850;
export const RACE_CATEGORY = 'Black';
export const CAT_ROW_TOP = 10;
export const CAT_ROW_PITCH = 27;
export const CAT_ROW_CENTRE = 13;

export async function cardCaptions(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
      .map((c) => c.textContent?.trim() ?? ''));
}

export async function visibleCardCaptions(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .filter((card) => (card as HTMLElement).offsetParent !== null)
      .map((card) => card.querySelector('.d4-filter-column-name')?.textContent?.trim() ?? ''));
}

export async function checkboxCensus(page: Page): Promise<{cards: number; boxes: number; checked: number}> {
  return page.evaluate(() => {
    const cards = Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'));
    const boxes = cards
      .map((c) => c.querySelector('input[type="checkbox"].ui-input-editor') as HTMLInputElement | null)
      .filter((cb) => cb !== null) as HTMLInputElement[];
    return {cards: cards.length, boxes: boxes.length, checked: boxes.filter((cb) => cb.checked).length};
  });
}

export async function filterState(page: Page, column: string, type: string):
    Promise<{selected: string[] | null; min: number | null; max: number | null} | null> {
  return page.evaluate(({column, type}) => {
    const states = grok.shell.tv.getFiltersGroup().getStates(column, type) as any[];
    if (!states || states.length === 0) return null;
    const s = states[0];
    return {
      selected: Array.isArray(s.selected) ? s.selected.map((c: any) => String(c)) : null,
      min: typeof s.min === 'number' ? s.min : null,
      max: typeof s.max === 'number' ? s.max : null,
    };
  }, {column, type});
}

// Waits for the row filter to leave `changedFrom` and settle, with the panel painted; a layout
// apply or a project open replaces the panel, so readiness is part of the wait.
export async function waitForPanelSettled(page: Page,
  opts: {changedFrom?: number; timeoutMs?: number} = {}): Promise<number> {
  const timeoutMs = opts.timeoutMs ?? 30_000;
  const ready = await page.evaluate((cap) => (window as any).__poll(() =>
    !!document.querySelector('[name="viewer-Filters"]') &&
      document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length > 0 &&
      (grok.shell.tv?.dataFrame?.rowCount ?? 0) > 0, (ok: boolean) => ok, cap, 50), timeoutMs);
  if (!ready) throw new Error(`waitForPanelSettled: the panel never painted within ${timeoutMs}ms`);
  const count = () => page.evaluate(() => grok.shell.tv?.dataFrame?.filter?.trueCount ?? -1);
  const from = opts.changedFrom;
  if (from !== undefined) {
    const moved = await v.pollValue(count, (c) => c !== from, timeoutMs, 50);
    if (moved === from) {
      throw new Error(`waitForPanelSettled: the panel never settled within ${timeoutMs}ms ` +
        `(last row count ${moved}, waiting for it to leave ${from})`);
    }
  }
  return page.evaluate(() => (window as any).__settledFor(
    () => grok.shell.tv.dataFrame.filter.trueCount, 250, 3000, 25));
}

export async function driveHeaderSearch(page: Page, text: string): Promise<{
  visibleBefore: string[]; visibleAfter: string[]; countBefore: number; countAfter: number; typed: string;
}> {
  const visibleBefore = await visibleCardCaptions(page);
  const countBefore = await trueCount(page);
  const input = page.locator('[name="viewer-Filters"] input.d4-search-input[placeholder="Search filters"]');
  if (!(await input.isVisible())) {
    const opened = await page.evaluate(() => {
      const icon = document.querySelector(
        '[name="viewer-Filters"] .d4-filter-group-header [name="icon-search"]') as HTMLElement | null;
      if (!icon) return false;
      icon.click();
      return true;
    });
    if (!opened)
      throw new Error('driveHeaderSearch: the panel header carries no search icon — the search cannot be opened');
  }
  await input.waitFor({state: 'visible', timeout: 15_000});
  await input.click();
  await page.keyboard.press('Control+a');
  await page.keyboard.type(text, {delay: 15});
  const typed = await v.pollValue(() => input.inputValue(), (val) => val === text, 800, 50);
  const wasVisible = visibleBefore.join(' ');
  await v.pollValue(async () => (await visibleCardCaptions(page)).join(' '),
    (now) => now !== wasVisible, 800, 50);
  return {
    visibleBefore,
    visibleAfter: await visibleCardCaptions(page),
    countBefore,
    countAfter: await trueCount(page),
    typed,
  };
}

export async function collapseHeaderSearch(page: Page): Promise<void> {
  const input = page.locator('[name="viewer-Filters"] input.d4-search-input[placeholder="Search filters"]');
  if (!(await input.isVisible()))
    throw new Error('collapseHeaderSearch: the header search box is not open, so it cannot be closed');
  const wasVisible = (await visibleCardCaptions(page)).join(' ');
  await page.evaluate(() => {
    const icon = document.querySelector(
      '[name="viewer-Filters"] .d4-filter-group-header [name="icon-search"]') as HTMLElement | null;
    if (!icon) throw new Error('collapseHeaderSearch: the panel header carries no search icon');
    icon.click();
  });
  await input.waitFor({state: 'hidden', timeout: 10_000});
  await v.pollValue(async () => (await visibleCardCaptions(page)).join(' '),
    (now) => now !== wasVisible, 800, 50);
}

export async function headerSearchState(page: Page): Promise<{value: string; visible: boolean} | null> {
  return page.evaluate(() => {
    const el = document.querySelector(
      '[name="viewer-Filters"] input.d4-search-input[placeholder="Search filters"]') as HTMLInputElement | null;
    return el === null ? null : {value: el.value, visible: el.offsetParent !== null};
  });
}

export async function removeAllViaPanelMenu(page: Page): Promise<void> {
  await v.drivePanelMenuLeaf(page, 'Filters', null, 'Remove All');
  await expect.poll(async () => cardCount(page),
    {timeout: 20_000, intervals: [30, 60, 120, 250, 500, 1000],
      message: 'the panel menu\'s "Remove All" leaf was driven but the panel still carries cards'})
    .toBe(0);
}

// The two-filter state of Step 10 / Step 13, set through the filter group's own API.
export async function establishTwoFilterState(page: Page): Promise<number> {
  return page.evaluate(async (category: string) => {
    const w = window as any;
    const DG = w.DG;
    const fg = grok.shell.tv.getFiltersGroup();
    const full = grok.shell.tv.dataFrame.filter.trueCount;
    const afterRace = await w.__filtered(() =>
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: [category]}), 600, full);
    return w.__filtered(() =>
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 18, max: 60}), 700, afterRace);
  }, RACE_CATEGORY);
}

// The row is addressed by the card's own row pitch, and only once the card body is tall enough
// to paint it — a card that just landed from a drag has no rows yet.
export async function categoryRowPoint(page: Page, column: string, category: string) {
  return page.evaluate(async ({column, category, top, pitch, centre}) => {
    const cats: string[] = grok.shell.tv.dataFrame.col(column).categories;
    const idx = cats.indexOf(category);
    if (idx < 0) return null;
    const overlayOf = () => Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === column)
      ?.querySelector('[name="viewer-Grid"] [name="overlay"]') as HTMLElement | null | undefined;
    const needed = top + pitch * (idx + 1);
    const canvasOf = () => Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .find((c) => c.querySelector('.d4-filter-column-name')?.textContent?.trim() === column)
      ?.querySelector('[name="viewer-Grid"] canvas[name="canvas"]') as HTMLCanvasElement | null | undefined;
    // the row exists once the card body is tall enough AND its band carries paint: a click that
    // reaches the overlay before the grid has drawn the row finds nothing under it
    const rowPainted = () => {
      const overlay = overlayOf();
      const cv = canvasOf();
      if (!overlay || !cv || overlay.getBoundingClientRect().height < needed || cv.width === 0) return false;
      const ctx = cv.getContext('2d');
      if (!ctx) return false;
      const r = cv.getBoundingClientRect();
      const sy = cv.height / r.height;
      const y0 = Math.round((top + pitch * idx) * sy);
      const h = Math.max(1, Math.round(pitch * sy));
      if (y0 + h > cv.height) return false;
      const d = ctx.getImageData(0, y0, cv.width, h).data;
      for (let i = 0; i < d.length; i += 4)
        if (d[i] < 200 || d[i + 1] < 200 || d[i + 2] < 200) return true;
      return false;
    };
    if (!(await (window as any).__poll(rowPainted, (ok: boolean) => ok, 3000, 50))) return null;
    const rect = overlayOf()!.getBoundingClientRect();
    return {x: rect.left + 12, y: rect.top + top + pitch * idx + centre};
  }, {column, category, top: CAT_ROW_TOP, pitch: CAT_ROW_PITCH, centre: CAT_ROW_CENTRE});
}
