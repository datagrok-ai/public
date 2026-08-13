/* ---
sub_features_covered: [legend.corner.pie, legend.corner.bar-chart, legend.trellis.default,
  legend.mini-icon.z-index, legend.mini-icon.small, legend.tooltip.hover-highlight]
--- */
// Automatic legend placement per viewer type, plus the mini icon and the tooltip legend.
// Pins: the pie chart's auto legend takes a free corner (marking its bounding box instead
// of the circle used to cover the very spot its offsets aim for); the bar chart's never
// lands on the bars (the candidate rect must be built in viewer coordinates, not chart-box
// ones); a freshly added trellis shows its legend without any property change (its
// requestLegendFrame must schedule a real render); the mini icon stays clickable above the
// line chart's full-width split strip; a too-small viewer shows the 20px icon; and
// hovering a tooltip legend item highlights the rows.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

async function legendState(page: any, viewerType: string) {
  return await page.evaluate((vt: string) => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === vt);
    if (!x) return {noViewer: true};
    const root = x.root.querySelector('[name="legend"]') as HTMLElement;
    if (!root) return {noLegend: true};
    const vr = x.root.getBoundingClientRect();
    const r = root.getBoundingClientRect();
    return {mode: root.getAttribute('data-legend-mode'), slot: root.getAttribute('data-legend-slot'),
      items: Number(root.getAttribute('data-legend-items')),
      box: {x: Math.round(r.x - vr.x), y: Math.round(r.y - vr.y),
        w: Math.round(r.width), h: Math.round(r.height)},
      viewer: {w: Math.round(vr.width), h: Math.round(vr.height)}};
  }, viewerType);
}

test('Pie corner, trellis default, bar chart placement (SPGI)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Pie chart auto placement takes a free corner beside the pie', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Pie chart', {categoryColumnName: 'Stereo Category'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Pie chart', 900, 600);
    await page.waitForTimeout(1500);
    const s = await legendState(page, 'Pie chart');
    expect(s.mode, 'pie legend stuck docked').toBe('corner');
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((q: any) => q.type === 'Pie chart')?.close();
    });
    await page.waitForTimeout(500);
  });

  await softStep('A freshly added trellis has its legend without any property change', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Trellis plot');
    });
    await page.waitForTimeout(4000);
    await v.resizeViewer(page, 'Trellis plot', 1000, 700);
    await page.waitForTimeout(2000);
    const s = await legendState(page, 'Trellis plot');
    expect(s.noLegend, 'trellis legend absent on default add').toBeFalsy();
    expect(s.mode).toBe('docked');
    expect(s.items).toBeGreaterThan(0);
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((q: any) => q.type === 'Trellis plot')?.close();
    });
    await page.waitForTimeout(500);
  });

  await softStep('Bar chart auto placement avoids the full-width top bar', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Bar chart', {
        splitColumnName: 'Stereo Category', stackColumnName: 'Primary Series Name'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Bar chart', 950, 700);
    await page.waitForTimeout(1500);
    const s = await legendState(page, 'Bar chart');
    expect(s.mode).toBe('corner');
    // with this data the top bar spans the full width, so the right-top corner is occupied;
    // the bottom bars are short, and the ladder must land there instead
    expect(s.slot, 'corner over the full-width top bar').not.toBe('rightTop');
  });

  expect(errors, `page errors: ${errors.join(' | ')}`).toHaveLength(0);
  v.finishSpec('Auto placement failures');
});

test('Mini icon clickability and sizes, tooltip hover highlight', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await page.evaluate(() => {
    const w = (window as any).grok.shell.windows;
    w.showContextPanel = false; w.showProperties = false; w.showHelp = false;
    const tv = (window as any).grok.shell.tv;
    tv.addViewer('Line chart', {xColumnName: 'Competition assay Date',
      yColumnNames: ['CAST Idea ID'],
      splitColumnNames: ['Primary Series Name', 'Scaffold Names', 'Series', 'Stereo Category']});
  });
  await page.waitForTimeout(2000);
  await v.resizeViewer(page, 'Line chart', 1100, 500);
  await page.waitForTimeout(800);
  await v.setViewerProps(page, 'Line chart', [{set: {legendPosition: 'RightTop'}, wait: 1500}]);
  await v.waitForLegendIdle(page, 'Line chart');

  await softStep('Collapsed at normal size: the 30px icon is clickable over the split strip', async () => {
    await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      legend.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      (document.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement).click();
    });
    await page.waitForTimeout(1200);
    const icon = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Line chart');
      const el = document.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      const vr = x.root.getBoundingClientRect();
      // hover the chart so the split-selector strip (full-width, z-index 1) is visible
      x.root.querySelector('.d4-layout-center')?.dispatchEvent(
        new MouseEvent('mouseenter', {bubbles: true}));
      const hit = document.elementFromPoint(r.x + r.width / 2, r.y + r.height / 2);
      return {w: Math.round(r.width), h: Math.round(r.height),
        display: el.style.display, insideViewer: r.x >= vr.x && r.right <= vr.right + 1,
        hitIsIcon: el === hit || el.contains(hit)};
    });
    expect(icon.display).toBe('block');
    expect(icon.w, 'user-collapsed icon must be the 30px one').toBe(30);
    expect(icon.hitIsIcon, 'mini icon buried under the split strip').toBe(true);
  });

  await softStep('A viewer too small for the legend shows the 20px icon', async () => {
    // the dock row's height is pinned by the grid above it, so shrink the width instead —
    // one dimension under minChartBox is enough for the small-icon state
    await v.resizeViewer(page, 'Line chart', 240, 400);
    await page.waitForTimeout(1500);
    const icon = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Line chart');
      const vr = x.root.getBoundingClientRect();
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      const el = document.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      return {mode: legend.getAttribute('data-legend-mode'), display: el.style.display,
        w: Math.round(r.width), viewer: {w: Math.round(vr.width), h: Math.round(vr.height)}};
    });
    expect(icon.mode).toBe('miniIcon');
    expect(icon.display).toBe('block');
    expect(icon.w, 'small-viewer icon must shrink to 20px').toBe(20);
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((q: any) => q.type === 'Line chart')?.close();
    });
    await page.waitForTimeout(500);
  });

  await softStep('Hovering a tooltip legend item highlights its rows', async () => {
    await page.evaluate(() => {
      for (const s of Array.from(document.querySelectorAll('style')))
        if (s.textContent?.includes('.d4-tooltip'))
          s.textContent = s.textContent.replace(/\.d4-tooltip\s*{[^}]*}/g, '');
      (window as any).grok.shell.tv.addViewer('Scatter plot', {
        colorColumnName: 'Primary Series Name'});
    });
    await page.waitForTimeout(2500);
    await v.resizeViewer(page, 'Scatter plot', 700, 500);
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'RightTop'}, wait: 1200}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      legend.dispatchEvent(new MouseEvent('mouseenter', {bubbles: false}));
      (document.querySelector('[name="icon-hide-corner-legend"]') as HTMLElement).click();
      const w = window as any;
      w.__hoverGroupEvents = 0;
      w.grok.shell.tv.dataFrame.onMouseOverRowGroupChanged.subscribe(() => w.__hoverGroupEvents++);
    });
    await page.waitForTimeout(1500);
    const mini = await page.evaluate(() => {
      const el = document.querySelector('[name="mini-legend-icon"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    await page.mouse.move(mini.x - 40, mini.y + 40);
    await page.mouse.move(mini.x, mini.y, {steps: 5});
    await page.waitForTimeout(1500);
    const item = await page.evaluate(() => {
      const legend = document.querySelector('[name="legend"]') as HTMLElement;
      if (!legend.closest('.d4-tooltip')) return null;
      const el = legend.querySelector('[name="legend-item"]') as HTMLElement;
      const r = el.getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    expect(item, 'tooltip legend did not appear').not.toBeNull();
    await page.mouse.move(item!.x, item!.y, {steps: 4});
    await page.waitForTimeout(800);
    const events = await page.evaluate(() => (window as any).__hoverGroupEvents);
    expect(events, 'hovering a tooltip legend item set no mouse-over group').toBeGreaterThan(0);
    await page.mouse.move(200, 900);
    await page.waitForTimeout(600);
  });

  expect(errors, `page errors: ${errors.join(' | ')}`).toHaveLength(0);
  v.finishSpec('Mini icon / tooltip failures');
});
