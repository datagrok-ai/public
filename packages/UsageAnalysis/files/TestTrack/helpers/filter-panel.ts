import {Page, expect} from '@playwright/test';

export const RESET_CRITERIA_ICON =
  '[name="viewer-Filters"] .d4-filter-group-header [name="icon-arrow-rotate-left"]';

export const HEADER_COUNTER =
  '[name="viewer-Filters"] .d4-filter-group-header .d4-filter-indicator';

export async function trueCount(page: Page): Promise<number> {
  const c = await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.filter.trueCount);
  expect(typeof c, 'trueCount must be a number, not null/undefined').toBe('number');
  return c;
}

export async function cardCount(page: Page): Promise<number> {
  return page.evaluate(() =>
    document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length);
}

// The panel header's reset icon clears the cards' criteria and KEEPS the cards. It is not
// viewers.resetFilters, which wipes the filter through the API and REMOVES every card; call
// sites assert that distinction. `via: 'dom'` dispatches the synthetic element click some
// callers rely on inside their own page.evaluate; the default is a real Playwright click.
export async function clickResetCriteriaIcon(
  page: Page,
  opts: {via?: 'locator' | 'dom'} = {},
): Promise<void> {
  if (opts.via === 'dom') {
    await page.evaluate(() => {
      const header = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header')!;
      (header.querySelector('[name="icon-arrow-rotate-left"]') as HTMLElement).click();
    });
    return;
  }
  await page.locator(RESET_CRITERIA_ICON).first().click();
}

export interface CounterTarget {
  present: boolean; visible: boolean; display: string; visibility: string;
  text: string; x: number; y: number; w: number; h: number;
}

// `present` is reported separately from `visible` on purpose: the panel hides the counter block
// while nothing filters, and a reader that conflates "absent from the DOM" with "hidden" passes
// on a counter that stopped rendering altogether.
export async function headerCounterTarget(page: Page): Promise<CounterTarget> {
  return page.evaluate((sel) => {
    const el = document.querySelector(sel) as HTMLElement | null;
    if (!el)
      return {present: false, visible: false, display: '', visibility: '', text: '', x: -1, y: -1, w: 0, h: 0};
    const cs = window.getComputedStyle(el);
    const r = el.getBoundingClientRect();
    return {
      present: true,
      visible: cs.display !== 'none' && cs.visibility !== 'hidden' && r.width > 0 && r.height > 0,
      display: cs.display,
      visibility: cs.visibility,
      text: el.textContent?.trim() ?? '',
      x: r.left + r.width / 2,
      y: r.top + r.height / 2,
      w: r.width,
      h: r.height,
    };
  }, HEADER_COUNTER);
}

// Readiness barrier for "the counter BECOMES x". Never use it for "the counter must not move":
// a settle poll on a no-change claim waits out exactly the defect it is meant to catch — use
// expectHeaderCounterNow / expectHeaderCounterQuietNow there.
export async function settledCounter(page: Page,
  ok: (t: CounterTarget) => boolean): Promise<CounterTarget> {
  let t = await headerCounterTarget(page);
  const deadline = Date.now() + 15_000;
  while (!ok(t) && Date.now() < deadline) {
    await page.waitForTimeout(250);
    t = await headerCounterTarget(page);
  }
  return t;
}

function assertCounter(t: CounterTarget, expected: string, why: string): void {
  expect(t.present, `${why} — the header active-filter counter is not in the DOM at all, so it ` +
    'reads neither a value nor "hidden"').toBe(true);
  expect(t.visible, `${why} — the counter block is present but off screen (display=${t.display}, ` +
    `visibility=${t.visibility}, box=${t.w}x${t.h}, text=${JSON.stringify(t.text)})`).toBe(true);
  expect(t.text, why).toBe(expected);
}

function assertCounterQuiet(t: CounterTarget, why: string): void {
  expect(t.present, `${why} — the counter element is gone from the DOM, which is neither "hidden ` +
    'while nothing filters" nor "reads 0"').toBe(true);
  expect(!t.visible || t.text === '0', `${why} — the counter is on screen reading ` +
    `${JSON.stringify(t.text)} (display=${t.display}, visibility=${t.visibility}, ` +
    `box=${t.w}x${t.h})`).toBe(true);
}

export async function expectHeaderCounter(page: Page, expected: string, why: string): Promise<void> {
  assertCounter(await settledCounter(page, (c) => c.present && c.visible && c.text === expected), expected, why);
}

// Single immediate read, for claims of the form "the counter must NOT move" / "must still be
// off screen". `expected === null` means the counter must be present and off screen.
export async function expectHeaderCounterNow(page: Page, expected: string | null, why: string): Promise<void> {
  const t = await headerCounterTarget(page);
  if (expected === null) {
    expect(t.present, `${why} — the header active-filter counter is not in the DOM at all, so a ` +
      'counter that stopped rendering entirely would satisfy every "hidden" check').toBe(true);
    expect(t.visible, `${why}; the counter block reads ${JSON.stringify(t.text)} (display=${t.display}, ` +
      `visibility=${t.visibility}, box=${t.w}x${t.h}) and must be off screen`).toBe(false);
    return;
  }
  assertCounter(t, expected, why);
}

export async function expectHeaderCounterQuiet(page: Page, why: string): Promise<void> {
  assertCounterQuiet(await settledCounter(page, (c) => c.present && (!c.visible || c.text === '0')), why);
}

export async function expectHeaderCounterQuietNow(page: Page, why: string): Promise<void> {
  assertCounterQuiet(await headerCounterTarget(page), why);
}

export async function cardCaptions(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name')]
      .map((e) => (e.textContent ?? '').trim()));
}

// The panel header's column selector: opens the picker, types the column into its search box,
// commits with Enter. The backdrop-closed poll at the end is load-bearing — a picker left open
// swallows the next gesture the caller makes.
export async function addCardViaColumnSelector(page: Page, column: string): Promise<void> {
  expect(await cardCaptions(page),
    `the panel must not already carry a ${column} card before the column selector adds it`)
    .not.toContain(column);
  const combo = await page.evaluate(() => {
    const el = document.querySelector(
      '[name="viewer-Filters"] .d4-filter-group-header [name="div-column-combobox-"]');
    if (!el) return null;
    const r = el.getBoundingClientRect();
    if (r.width === 0 || r.height === 0) return null;
    return {x: r.left + r.width / 2, y: r.top + r.height / 2};
  });
  expect(combo, 'the panel header exposes no column selector to add a filter with').not.toBeNull();
  await page.mouse.click(combo!.x, combo!.y);
  await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout: 15_000});
  await page.keyboard.press('a');
  await page.waitForFunction(() => !!document.querySelector('input.d4-column-selector-search-input'),
    null, {timeout: 10_000});
  await page.keyboard.press('Control+a');
  await page.keyboard.type(column, {delay: 40});
  expect(await page.evaluate(() =>
    (document.querySelector('input.d4-column-selector-search-input') as HTMLInputElement).value),
  `the column picker search box does not hold ${column}, so Enter would commit some other pick`)
    .toBe(column);
  await page.keyboard.press('Enter');
  await expect.poll(async () => (await cardCaptions(page)).includes(column), {
    timeout: 20_000,
    intervals: [500, 1000, 2000],
    message: `no ${column} card came out of the column-picker commit`,
  }).toBe(true);
  await expect.poll(async () => page.locator('.d4-column-selector-backdrop').count(), {
    timeout: 10_000,
    message: 'the column picker popup stayed open and would intercept the next gesture',
  }).toBe(0);
}

// Drives one leaf of the popup menu that is ALREADY open (the last one on screen), optionally
// hovering a group first. The caller owns opening the menu.
export async function driveOpenMenuLeaf(page: Page, group: string | null, leaf: string): Promise<void> {
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 15_000});
  await page.evaluate(async ({g, l}: {g: string | null; l: string}) => {
    const norm = (s: string | null | undefined) => (s ?? '').trim().toLowerCase();
    const labels = () => {
      const popup = [...document.querySelectorAll('.d4-menu-popup')].pop();
      return [...(popup?.querySelectorAll('.d4-menu-item-label') ?? [])];
    };
    const visible = () => labels().map((i) => (i.textContent ?? '').trim()).join(' | ');
    const find = (text: string) =>
      labels().find((i) => norm(i.textContent) === norm(text))?.closest('.d4-menu-item') ?? null;
    const waitFor = async (text: string) => {
      const deadline = Date.now() + 5000;
      let item = find(text);
      while (!item && Date.now() < deadline) {
        await new Promise((r) => setTimeout(r, 100));
        item = find(text);
      }
      return item;
    };
    if (g !== null) {
      const groupItem = await waitFor(g);
      if (!groupItem)
        throw new Error(`open menu: group "${g}" not found in the last popup; visible: ${visible()}`);
      const b = groupItem.getBoundingClientRect();
      for (const type of ['mouseover', 'mousemove'])
        groupItem.dispatchEvent(new MouseEvent(type, {bubbles: true, clientX: b.x + 5, clientY: b.y + 5}));
    }
    const leafItem = await waitFor(l);
    if (!leafItem)
      throw new Error(`open menu: leaf "${l}" not found in the last popup; visible: ${visible()}`);
    (leafItem as HTMLElement).click();
  }, {g: group, l: leaf});
}
