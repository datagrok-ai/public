import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

interface TooltipSnapshot { display: string; text: string; cols: string[] }

async function readTooltipState(page: import('@playwright/test').Page): Promise<TooltipSnapshot> {
  return page.evaluate(() => {
    const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
    return {
      display: tt ? getComputedStyle(tt).display : 'missing',
      text: tt?.textContent ?? '',
      cols: tt ? Array.from(tt.querySelectorAll('.d4-row-tooltip-table tr td:first-child')).map((td) => td.textContent ?? '') : [],
    };
  });
}

async function captureViewerTooltipBySweep(
  page: import('@playwright/test').Page,
  viewerSelector: string,
  canvasIdx: 'data' | 'overlay',
  sweep: {dx: number; dy: number}[],
): Promise<TooltipSnapshot> {
  const box = await page.evaluate(({sel, idx}) => {
    const v = document.querySelector(sel)!;
    const canvases = v.querySelectorAll('canvas');
    const c = (idx === 'data' ? canvases[0] : canvases[canvases.length - 1]) as HTMLCanvasElement;
    const r = c.getBoundingClientRect();
    return {left: r.left, top: r.top, width: r.width, height: r.height};
  }, {sel: viewerSelector, idx: canvasIdx});
  for (const {dx, dy} of sweep) {
    await page.evaluate(() => { try { (window as any).ui?.tooltip?.hide?.(); } catch (_) {} });
    await page.mouse.move(1, 1);
    await page.waitForTimeout(150);
    await page.mouse.move(box.left + dx - 25, box.top + dy - 25);
    await page.waitForTimeout(120);
    await page.mouse.move(box.left + dx, box.top + dy);
    await page.waitForTimeout(1100);
    const snap = await readTooltipState(page);
    if (snap.display === 'block' && snap.cols.length > 0)
      return snap;
  }
  return readTooltipState(page);
}

async function captureScatterPlotTooltip(page: import('@playwright/test').Page): Promise<TooltipSnapshot> {
  const sweep = [
    {dx: 287, dy: 222},   // approx center; scatter plots are dense at center
    {dx: 250, dy: 200},
    {dx: 320, dy: 250},
    {dx: 200, dy: 180},
    {dx: 350, dy: 280},
  ];
  return captureViewerTooltipBySweep(page, '[name="viewer-Scatter-plot"]', 'data', sweep);
}

async function captureBoxPlotTooltip(page: import('@playwright/test').Page): Promise<TooltipSnapshot> {
  // Sweep the first category column where markers cluster
  const sweep: {dx: number; dy: number}[] = [];
  for (let dx = 90; dx <= 180; dx += 12)
    for (let dy = 80; dy <= 320; dy += 18)
      sweep.push({dx, dy});
  return captureViewerTooltipBySweep(page, '[name="viewer-Box-plot"]', 'overlay', sweep);
}

async function captureGridTooltip(page: import('@playwright/test').Page): Promise<TooltipSnapshot> {
  const sweep = [
    {dx: 200, dy: 200},
    {dx: 300, dy: 250},
    {dx: 100, dy: 150},
  ];
  return captureViewerTooltipBySweep(page, '[name="viewer-Grid"]', 'data', sweep);
}

test('Viewers: uniform default tooltip', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (grok as any).shell.settings.showFiltersIconsConstantly = true;
    (grok as any).shell.windows.simpleMode = true;
    (grok as any).shell.closeAll();
    const df = await (grok as any).dapi.files.readCsv('System:DemoFiles/demog.csv');
    (grok as any).shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Open a scatter plot and a box plot', async () => {
    await page.evaluate(() => {
      (document.querySelector('[name="icon-scatter-plot"]') as HTMLElement).click();
    });
    await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 10000});
    await page.evaluate(() => {
      (document.querySelector('[name="icon-box-plot"]') as HTMLElement).click();
    });
    await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
    await page.waitForTimeout(2000);
    const types = await page.evaluate(() => (grok as any).shell.tv.viewers.map((v: any) => v.type));
    expect(types).toEqual(expect.arrayContaining(['Grid', 'Scatter plot', 'Box plot']));
  });

  await softStep('Hover viewers (incl. main grid): tooltip is same column set in same order', async () => {
    const spTip = await captureScatterPlotTooltip(page);
    expect(spTip.display).toBe('block');
    expect(spTip.cols.length).toBeGreaterThan(0);

    const bpTip = await captureBoxPlotTooltip(page);
    expect(bpTip.display).toBe('block');
    expect(bpTip.cols.length).toBeGreaterThan(0);

    const gridTip = await captureGridTooltip(page);
    // Scenario expects the grid to also produce a tooltip with the same columns
    expect(gridTip.display).toBe('block');
    expect(gridTip.cols.length).toBeGreaterThan(0);

    // Same set of columns
    const spSet = [...spTip.cols].sort();
    const bpSet = [...bpTip.cols].sort();
    const gridSet = [...gridTip.cols].sort();
    expect(bpSet).toEqual(spSet);
    expect(gridSet).toEqual(spSet);

    // Same order
    expect(bpTip.cols).toEqual(spTip.cols);
    expect(gridTip.cols).toEqual(spTip.cols);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((s) => `  - ${s.step}: ${s.error}`).join('\n'));
});
