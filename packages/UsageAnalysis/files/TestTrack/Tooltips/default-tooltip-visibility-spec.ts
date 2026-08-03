import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function openTooltipSubmenu(page: import('@playwright/test').Page, viewerSelector: string) {
  await page.evaluate(async (sel) => {
    const viewer = document.querySelector(sel)!;
    const canvas = viewer.querySelector('canvas')! as HTMLElement;
    const r = canvas.getBoundingClientRect();
    const cx = r.left + r.width / 2;
    const cy = r.top + r.height / 2;
    for (const type of ['mousedown', 'mouseup', 'contextmenu'])
      canvas.dispatchEvent(new MouseEvent(type, {bubbles: true, cancelable: true, button: 2, buttons: 2, clientX: cx, clientY: cy}));
    await new Promise((r) => setTimeout(r, 350));
    const tooltipGroup = Array.from(document.querySelectorAll('[name="div-Tooltip"]'))
      .find((el) => (el as HTMLElement).offsetParent !== null) as HTMLElement;
    const trect = tooltipGroup.getBoundingClientRect();
    const opts = (x: number, y: number) => ({bubbles: true, cancelable: true, clientX: x, clientY: y});
    tooltipGroup.dispatchEvent(new MouseEvent('mouseover', opts(trect.left + 1, trect.top + trect.height / 2)));
    tooltipGroup.dispatchEvent(new MouseEvent('mouseenter', opts(trect.left + 1, trect.top + trect.height / 2)));
    for (let i = 1; i <= 8; i++) {
      const x = trect.left + (trect.width * i) / 8;
      tooltipGroup.dispatchEvent(new MouseEvent('mousemove', opts(x, trect.top + trect.height / 2)));
      await new Promise((r) => setTimeout(r, 30));
    }
    await new Promise((r) => setTimeout(r, 400));
  }, viewerSelector);
}

test('Viewers: default tooltip visibility', async ({page}) => {
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
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await softStep('Enable Show Visible Columns In Tooltip on grid', async () => {
    await page.evaluate(() => {
      const grid = document.querySelector('[name="viewer-Grid"]')!;
      (grid.querySelector('[name="icon-font-icon-settings"]') as HTMLElement).click();
    });
    await page.waitForTimeout(300);
    await page.evaluate(() => {
      const row = document.querySelector('[name="prop-show-visible-columns-in-tooltip"]')!;
      const checkbox = row.querySelector('input[type="checkbox"]') as HTMLInputElement;
      if (!checkbox.checked) checkbox.click();
    });
    await page.waitForTimeout(200);
    const val = await page.evaluate(() => (grok as any).shell.tv.grid.props.showVisibleColumnsInTooltip);
    expect(val).toBe(true);
  });

  await softStep('Open Scatter plot, Box plot, Histogram, Line Chart, Bar chart, Trellis plot', async () => {
    const icons = [
      'icon-scatter-plot', 'icon-box-plot', 'icon-histogram',
      'icon-line-chart', 'icon-bar-chart', 'icon-trellis-plot',
    ];
    for (const name of icons) {
      await page.evaluate((n) => (document.querySelector(`[name="${n}"]`) as HTMLElement).click(), name);
      await page.waitForTimeout(200);
    }
    await page.waitForTimeout(1500);
    const types = await page.evaluate(() =>
      Array.from((grok as any).shell.tv.viewers).map((v: any) => v.type));
    expect(types).toContain('Scatter plot');
    expect(types).toContain('Box plot');
    expect(types).toContain('Histogram');
    expect(types).toContain('Line chart');
    expect(types).toContain('Bar chart');
    expect(types).toContain('Trellis plot');
  });

  await softStep('Right-click viewer → Tooltip > Hide; dataframe tooltip-visibility becomes false', async () => {
    await openTooltipSubmenu(page, '[name="viewer-Scatter-plot"]');
    const hideVisible = await page.evaluate(() =>
      !!Array.from(document.querySelectorAll('[name="div-Tooltip---Hide"]'))
        .find((el) => (el as HTMLElement).offsetParent !== null));
    expect(hideVisible).toBe(true);
    await page.evaluate(() => {
      const hide = Array.from(document.querySelectorAll('[name="div-Tooltip---Hide"]'))
        .find((el) => (el as HTMLElement).offsetParent !== null) as HTMLElement;
      hide.click();
    });
    await page.waitForTimeout(400);
    const tooltipVisibility = await page.evaluate(() =>
      (grok as any).shell.t.tags.get('.tooltip-visibility'));
    expect(tooltipVisibility).toBe('false');
  });

  await softStep('Right-click viewer → Tooltip > Show Custom (replaces Hide); tooltip-visibility becomes true', async () => {
    await openTooltipSubmenu(page, '[name="viewer-Scatter-plot"]');
    const showCustomVisible = await page.evaluate(() =>
      !!Array.from(document.querySelectorAll('[name="div-Tooltip---Show-Custom"]'))
        .find((el) => (el as HTMLElement).offsetParent !== null));
    const hideVisible = await page.evaluate(() =>
      !!Array.from(document.querySelectorAll('[name="div-Tooltip---Hide"]'))
        .find((el) => (el as HTMLElement).offsetParent !== null));
    expect(showCustomVisible).toBe(true);
    expect(hideVisible).toBe(false);
    await page.evaluate(() => {
      const sc = Array.from(document.querySelectorAll('[name="div-Tooltip---Show-Custom"]'))
        .find((el) => (el as HTMLElement).offsetParent !== null) as HTMLElement;
      sc.click();
    });
    await page.waitForTimeout(400);
    const tooltipVisibility = await page.evaluate(() =>
      (grok as any).shell.t.tags.get('.tooltip-visibility'));
    expect(tooltipVisibility).toBe('true');
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((s) => `  - ${s.step}: ${s.error}`).join('\n'));
});
