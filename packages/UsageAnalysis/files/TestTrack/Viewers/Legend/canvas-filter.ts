import {Page} from '@playwright/test';
import * as v from '../../helpers/viewers';

export interface CanvasFilterResult {
  before: number;
  totalFiltered: number;
  survivors: number;
  clickedAt: {x: number; y: number};
}

/**
 * Puts the row filter back to every row: the Filter Panel entries, the scatter plot's own filter
 * expression and the click-to-filter state a bar, pie or trellis viewer holds. A viewer drops its
 * click filter when its onClick leaves Filter, which is the only API way to clear it.
 */
export async function clearClickFilters(page: Page): Promise<void> {
  await v.installEventWaits(page);
  await page.evaluate(async () => {
    const w = window as any;
    const tv = w.grok.shell.tv;
    const df = tv.dataFrame;
    const fg = tv.getFiltersGroup();
    for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
    for (const x of tv.viewers) {
      if (x.type === 'Scatter plot') try { x.props.filter = ''; } catch (_) {}
      if (['Bar chart', 'Pie chart', 'Trellis plot'].includes(x.type) && x.props.onClick === 'Filter')
        x.props.onClick = x.type === 'Trellis plot' ? 'None' : 'Select';
    }
    df.filter.setAll(true);
    await w.__poll(() => df.filter.trueCount, (n: number) => n === df.rowCount, 1500, 25);
  });
}

/**
 * Clears every filter, arms click-to-filter on the viewer and clicks ONE painted bar or pie
 * slice: the saturated pixel (with saturated neighbours, so not an anti-aliased edge) nearest
 * the canvas centre that the canvas itself receives at that point. Grid lines, labels and the
 * background are grey or white, so saturation alone separates a bar or slice from the rest.
 * Resolves once the row filter has moved off its pre-click count; throws if it never does.
 */
export async function clickCanvasFilter(
  page: Page, opts: {viewerType: 'Bar chart' | 'Pie chart'; column: string},
): Promise<CanvasFilterResult> {
  await clearClickFilters(page);
  const target = await page.evaluate(async (vt) => {
    const w = window as any;
    const tv = w.grok.shell.tv;
    const df = tv.dataFrame;
    const viewer = tv.viewers.find((x: any) => x.type === vt);
    viewer.props.onClick = 'Filter';
    await w.__quiet(`viewer:${vt}.onViewerRendered`, 200, 1500);

    const cv = (Array.from(viewer.root.querySelectorAll('canvas')) as HTMLCanvasElement[])
      .sort((a, b) => b.width * b.height - a.width * a.height)[0];
    const r = cv.getBoundingClientRect();
    const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
    const saturated = (px: number, py: number) => {
      if (px < 0 || py < 0 || px >= cv.width || py >= cv.height) return false;
      const i = (py * cv.width + px) * 4;
      return img[i + 3] > 200 && Math.max(img[i], img[i + 1], img[i + 2]) - Math.min(img[i], img[i + 1], img[i + 2]) > 40;
    };
    const painted = (px: number, py: number) => saturated(px, py) &&
      saturated(px - 2, py) && saturated(px + 2, py) && saturated(px, py - 2) && saturated(px, py + 2);
    const sx = cv.width / r.width;
    const sy = cv.height / r.height;
    const cx = r.width / 2;
    const cy = r.height / 2;
    const maxR = Math.hypot(cx, cy);
    for (let rad = 0; rad <= maxR; rad += 3) {
      const n = rad === 0 ? 1 : 24;
      for (let k = 0; k < n; k++) {
        const a = (2 * Math.PI * k) / n;
        const x = cx + rad * Math.cos(a);
        const y = cy + rad * Math.sin(a);
        if (x < 3 || y < 3 || x > r.width - 3 || y > r.height - 3) continue;
        if (!painted(Math.round(x * sx), Math.round(y * sy))) continue;
        const clientX = r.left + x;
        const clientY = r.top + y;
        if (document.elementFromPoint(clientX, clientY) !== cv) continue;
        return {x: clientX, y: clientY, before: df.filter.trueCount};
      }
    }
    return null;
  }, opts.viewerType);
  if (!target) throw new Error(`${opts.viewerType}: no painted bar or slice found on the canvas`);

  await page.mouse.click(target.x, target.y);

  const after = await page.evaluate(async ({col, before}) => {
    const w = window as any;
    const df = w.grok.shell.tv.dataFrame;
    const totalFiltered: number = await w.__moved(() => df.filter.trueCount, before, 3000);
    const c = df.col(col);
    const counts: Record<string, number> = {};
    for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) counts[c.get(i)] = (counts[c.get(i)] ?? 0) + 1;
    return {totalFiltered, survivors: Object.keys(counts).length};
  }, {col: opts.column, before: target.before});
  if (after.totalFiltered === target.before)
    throw new Error(`${opts.viewerType}: click at (${Math.round(target.x)}, ${Math.round(target.y)}) did not change the row filter (${target.before} rows)`);
  return {before: target.before, ...after, clickedAt: {x: target.x, y: target.y}};
}
