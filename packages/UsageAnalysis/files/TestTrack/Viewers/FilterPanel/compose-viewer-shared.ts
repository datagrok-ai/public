import {Page} from '@playwright/test';
import * as v from '../../helpers/viewers';

declare const grok: any;

export const FULL = 5850;
export const PANEL_CATEGORY = 'Asian';

export async function raceSelectedCategories(page: Page): Promise<string[] | null> {
  return page.evaluate(() => {
    const states = grok.shell.tv.getFiltersGroup().getStates('RACE', 'categorical') as any[];
    if (!states || states.length === 0) return null;
    const s = states[0];
    return Array.isArray(s.selected) ? s.selected.map((c: any) => String(c)) : null;
  });
}

export async function seedPanelCriterion(page: Page, category: string): Promise<number> {
  await page.evaluate(() => { grok.shell.tv.dataFrame.selection.setAll(false); });
  await v.resetFilters(page, {clearScatterFilter: true});
  const {filteredCount} = await v.applyCategoricalFilter(page, 'RACE', [category]);
  return filteredCount;
}

export async function addViewer(page: Page, type: string): Promise<void> {
  await page.evaluate(async (t: string) => {
    const w = window as any;
    const viewer = grok.shell.tv.addViewer(t);
    await w.__poll(() => Array.from(grok.shell.tv.viewers).includes(viewer),
      (there: boolean) => there, 1500, 25);
    await w.__eventFired(`viewer:${t}.onViewerRendered`, 1500).catch(() => {});
  }, type);
}

export async function viewerCanvasRect(page: Page, type: string):
    Promise<{x: number; y: number; w: number; h: number} | null> {
  return page.evaluate((t: string) => {
    const vw = grok.shell.tv.viewers.find((x: any) => x.type === t);
    if (!vw) return null;
    const cv = vw.root.querySelector('canvas[name="canvas"]') || vw.root.querySelector('canvas');
    if (!cv) return null;
    const r = cv.getBoundingClientRect();
    return {x: r.x, y: r.y, w: r.width, h: r.height};
  }, type);
}

export async function zoomScatterPlot(page: Page, rect: {x: number; y: number; w: number; h: number}): Promise<void> {
  const x1 = rect.x + rect.w * 0.30, y1 = rect.y + rect.h * 0.30;
  const x2 = rect.x + rect.w * 0.70, y2 = rect.y + rect.h * 0.70;
  await page.mouse.move(x1, y1);
  await page.mouse.down();
  await page.mouse.move((x1 + x2) / 2, (y1 + y2) / 2, {steps: 2});
  await page.mouse.move(x2, y2, {steps: 3});
  await page.mouse.up();
}
