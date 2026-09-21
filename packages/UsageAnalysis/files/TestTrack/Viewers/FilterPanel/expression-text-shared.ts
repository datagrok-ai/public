import {expect, Page} from '@playwright/test';
import * as v from '../../helpers/viewers';
import {cardCount} from '../../helpers/filter-panel';

declare const grok: any;

export async function filterSummaryText(page: Page): Promise<string> {
  return page.evaluate(() => {
    const el = grok.shell.tv.getFiltersGroup().getFilterSummary();
    return (el?.textContent ?? '').trim();
  });
}

export async function removeAllCards(page: Page): Promise<void> {
  await v.drivePanelMenuLeaf(page, 'Filters', null, 'Remove All');
  await expect.poll(async () => cardCount(page), {timeout: 20_000, intervals: [400, 800, 1500]}).toBe(0);
}

export async function andOrControl(page: Page, column: string | null, toggle: boolean): Promise<string> {
  return page.evaluate(async ({col, doToggle}) => {
    const name = col === null ? 'expression' : col;
    const card = col === null
      ? (document.querySelector('.d4-expression-filter') as HTMLElement | null)?.closest('.d4-filter')
      : [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
        .find((c) => ((c.querySelector('.d4-filter-column-name'))?.textContent ?? '').trim() === col);
    if (!card) throw new Error(`the panel carries no ${name} card to read its AND/OR control on`);
    const header = card.querySelector('.d4-filter-header');
    if (!header) throw new Error(`the ${name} card carries no .d4-filter-header`);
    const read = () => [...header.querySelectorAll('div')]
      .find((d) => { const t = (d.textContent ?? '').trim(); return t === 'OR' || t === 'AND'; }) as HTMLElement | undefined;
    const before = read();
    if (!before) throw new Error(`the ${name} card header carries no AND/OR control to read or click`);
    const was = (before.textContent ?? '').trim();
    if (!doToggle) return was;
    before.click();
    return (window as any).__poll(() => (read()?.textContent ?? '').trim(),
      (t: string) => t !== '' && t !== was, 1500, 50);
  }, {col: column, doToggle: toggle});
}
