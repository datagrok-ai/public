import {expect, Page} from '@playwright/test';

declare const grok: any;

export const FULL = 5850;

export async function installViewResolver(page: Page): Promise<void> {
  await page.evaluate(() => {
    (window as any).__tv = (vn: string) => {
      const tvs: any[] = [];
      for (const x of grok.shell.tableViews) tvs.push(x);
      const hits = tvs.filter((t: any) => t.name === vn);
      if (hits.length !== 1)
        throw new Error(`TableView "${vn}" matched ${hits.length} views (open: ${tvs.map((t: any) => t.name).join(' | ')})`);
      return hits[0];
    };
    (window as any).__cards = (vn: string): HTMLElement[] => {
      const root = (window as any).__tv(vn).root as HTMLElement;
      return Array.from(root.querySelectorAll('[name="viewer-Filters"] .d4-filter')) as HTMLElement[];
    };
    (window as any).__card = (vn: string, caption: string): HTMLElement => {
      const hits = (window as any).__cards(vn).filter((x: HTMLElement) =>
        x.querySelector('.d4-filter-column-name')?.textContent?.trim() === caption);
      if (hits.length !== 1)
        throw new Error(`the Filter Panel of view "${vn}" paints ${hits.length} cards captioned "${caption}", expected exactly one`);
      return hits[0];
    };
  });
}

export async function activateView(page: Page, viewName: string): Promise<void> {
  await page.evaluate((vn: string) => { grok.shell.v = (window as any).__tv(vn); }, viewName);
  await expect.poll(async () => page.evaluate((vn: string) => {
    const el = ((window as any).__tv(vn).root as HTMLElement)
      .querySelector('[name="viewer-Filters"]') as HTMLElement | null;
    return el != null && el.offsetParent !== null;
  }, viewName), {
    message: `the Filter Panel of view "${viewName}" never became the laid-out one after switching to it — `
      + 'every :visible panel read below would have addressed the other view',
    timeout: 15_000, intervals: [30, 60, 120, 250, 500, 1000],
  }).toBe(true);
}

export async function clickCardCheckboxIn(page: Page, viewName: string, column: string): Promise<void> {
  await activateView(page, viewName);
  await page.evaluate(({vn, col}) => {
    const box = (window as any).__card(vn, col)
      .querySelector('input[type="checkbox"].ui-input-editor') as HTMLElement | null;
    if (!box)
      throw new Error(`the "${col}" card in view "${vn}" carries no enable/disable checkbox — the toggle gesture never happened`);
    box.click();
  }, {vn: viewName, col: column});
}

export async function trueCountOf(page: Page, viewName: string): Promise<number> {
  return page.evaluate((vn: string) => (window as any).__tv(vn).dataFrame.filter.trueCount, viewName);
}
