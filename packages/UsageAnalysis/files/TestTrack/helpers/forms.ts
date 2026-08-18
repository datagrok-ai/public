import {Page, expect} from '@playwright/test';

export const HOST = '[name="viewer-Forms"]';

export const ORDINARY = `${HOST} #vlist .d4-multi-form-form:not(.temp)`;

export const CURRENT = `${ORDINARY}:has(.d4-multi-form-form-indicator-current-row)`;

export const PINNED_PANE = `${HOST} .d4-multi-form-pinned-forms`;

export const PINNED = `${PINNED_PANE} .d4-multi-form-form`;

export async function drawnLabelNames(page: Page): Promise<string[]> {
  return page.evaluate((host) => Array.from(document.querySelectorAll(
    `${host} .d4-multi-form-header .d4-multi-form-column-name div[name^="div-"]`))
    .map((e) => (e.textContent ?? '').trim()).filter((n) => n.length > 0), HOST);
}

export async function balloonCount(page: Page): Promise<number> {
  return page.locator('.d4-balloon').count();
}

export async function cardFieldValue(
  page: Page, cardIndex: number, column: string, cardSel: string = ORDINARY,
): Promise<string | null> {
  return page.evaluate(({sel, i, col}) => {
    const card = document.querySelectorAll(sel)[i] as HTMLElement | undefined;
    const el = card?.querySelector(`[column="${col}"]`) as HTMLInputElement | null;
    return el ? el.value : null;
  }, {sel: cardSel, i: cardIndex, col: column});
}

export async function fieldValuesByPosition(
  page: Page, column: string, cardSel: string = ORDINARY,
): Promise<(string | null)[]> {
  return page.evaluate(({sel, col}) => Array.from(document.querySelectorAll(sel))
    .map((c) => ((c as HTMLElement).querySelector(`[column="${col}"]`) as HTMLInputElement)?.value ?? null),
  {sel: cardSel, col: column});
}

export async function cardIndexByValue(
  page: Page, column: string, value: string, cardSel: string = ORDINARY,
): Promise<number> {
  return page.evaluate(({sel, col, val}) => {
    const cards = Array.from(document.querySelectorAll(sel));
    for (let i = 0; i < cards.length; i++) {
      const el = (cards[i] as HTMLElement).querySelector(`[column="${col}"]`) as HTMLInputElement | null;
      if (el && el.value === val) return i;
    }
    return -1;
  }, {sel: cardSel, col: column, val: value});
}

export async function withConsoleErrorCount(
  page: Page, fn: () => Promise<void>, settleMs = 400, texts?: string[],
): Promise<number> {
  let count = 0;
  const handler = (msg: {type(): string; text(): string}) => {
    if (msg.type() === 'error') { count++; if (texts) texts.push(msg.text()); }
  };
  page.on('console', handler);
  try {
    await fn();
    await page.waitForTimeout(settleMs);
  } finally {
    page.off('console', handler);
  }
  return count;
}

export async function waitForOrderStable(
  page: Page, opts: {column?: string; cardSel?: string; timeoutMs?: number} = {},
): Promise<void> {
  const column = opts.column ?? 'USUBJID';
  const cardSel = opts.cardSel ?? ORDINARY;
  let prev: string | null = null;
  await expect.poll(async () => {
    const cur = await page.evaluate(({sel, col}) => JSON.stringify(
      Array.from(document.querySelectorAll(sel))
        .map((c) => ((c as HTMLElement).querySelector(`[column="${col}"]`) as HTMLInputElement)?.value ?? null)
        .filter((x) => x !== null)), {sel: cardSel, col: column});
    const stable = prev !== null && cur === prev;
    prev = cur;
    return stable;
  }, {timeout: opts.timeoutMs ?? 20_000, intervals: [250, 250, 250, 300, 500]}).toBe(true);
}
