import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const EXPECTED_TOOLTIP_ITEMS = ['Hide', 'Edit...', 'Use as Group Tooltip', 'Remove Group Tooltip'];

async function readTooltipMenu(page: import('@playwright/test').Page, viewerSelector: string) {
  // Close any previously-open menu first
  await page.keyboard.press('Escape').catch(() => {});
  await page.mouse.click(2, 2).catch(() => {});
  await page.waitForTimeout(300);

  // Use a real right-click via page.mouse so Dart's full menu pipeline runs
  const v = page.locator(viewerSelector);
  await v.waitFor({timeout: 10_000});
  const canvasOrViewer = page.locator(`${viewerSelector} canvas`).first();
  const target = (await canvasOrViewer.count()) > 0 ? canvasOrViewer : v;
  const box = await target.boundingBox();
  if (!box) return {err: 'no bounding box'};
  const cx = box.x + box.width / 2;
  const cy = box.y + box.height / 2;
  await page.mouse.move(cx, cy);
  await page.mouse.click(cx, cy, {button: 'right'});
  await page.locator('.d4-menu-popup').waitFor({timeout: 5_000}).catch(() => {});
  await page.waitForTimeout(400);

  // Hover the Tooltip group so Dart mounts the submenu DOM, then read labels
  return page.evaluate(async () => {
    const popups = document.querySelectorAll('.d4-menu-popup');
    const popup = popups[popups.length - 1] as HTMLElement | undefined;
    if (!popup) return {err: 'no popup'};

    const tooltipGroup = Array.from(popup.querySelectorAll('.d4-menu-item.d4-menu-group'))
      .find((it) => (it.querySelector(':scope > .d4-menu-item-label')?.textContent ?? '').trim() === 'Tooltip') as HTMLElement | undefined;
    if (!tooltipGroup) {
      const topGroups = Array.from(popup.children).filter((el) => el.classList.contains('d4-menu-item'))
        .map((it) => (it.querySelector(':scope > .d4-menu-item-label')?.textContent ?? '').trim());
      return {err: 'no Tooltip section', topGroups};
    }

    // Dispatch hover events so the submenu is mounted/expanded
    tooltipGroup.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
    tooltipGroup.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    tooltipGroup.dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
    await new Promise((r) => setTimeout(r, 400));

    // Match the MCP semantics: any descendant submenu container under the Tooltip group.
    const containers = tooltipGroup.querySelectorAll('.d4-menu-item-container');
    if (containers.length === 0) return {err: 'no submenu container', html: tooltipGroup.outerHTML.slice(0, 500)};
    // The first container is the direct submenu of "Tooltip"; subsequent ones belong to nested groups.
    const labels = Array.from(containers[0].querySelectorAll(':scope > .d4-menu-item > .d4-menu-item-label'))
      .map((l) => (l.textContent ?? '').trim());
    return {labels};
  });
}

test('Viewers: actions in the context menu', async ({page}) => {
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

  await softStep('Add Histogram, Line Chart, Bar Chart, Trellis plot', async () => {
    await page.evaluate(async () => {
      const icons = ['icon-histogram', 'icon-line-chart', 'icon-bar-chart', 'icon-trellis-plot'];
      for (const name of icons) {
        const el = document.querySelector(`[name="${name}"]`) as HTMLElement | null;
        if (el) el.click();
        await new Promise((r) => setTimeout(r, 400));
      }
      await new Promise((r) => setTimeout(r, 800));
    });
    await page.locator('[name="viewer-Histogram"]').waitFor({timeout: 10_000});
    await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 10_000});
    await page.locator('[name="viewer-Bar-chart"]').waitFor({timeout: 10_000});
    await page.locator('[name="viewer-Trellis-plot"]').waitFor({timeout: 10_000});
    const count = await page.evaluate(() => (grok as any).shell.tv.viewers.length);
    expect(count).toBe(5);
  });

  const viewers: Array<{name: string; selector: string}> = [
    {name: 'Grid', selector: '[name="viewer-Grid"]'},
    {name: 'Histogram', selector: '[name="viewer-Histogram"]'},
    {name: 'Line chart', selector: '[name="viewer-Line-chart"]'},
    {name: 'Bar chart', selector: '[name="viewer-Bar-chart"]'},
    {name: 'Trellis plot', selector: '[name="viewer-Trellis-plot"]'},
  ];

  for (const v of viewers) {
    await softStep(`${v.name}: Tooltip section contains the four expected items`, async () => {
      const res = await readTooltipMenu(page, v.selector) as {err?: string; labels?: string[]};
      expect(res.err, `${v.name}: ${res.err ?? ''}`).toBeUndefined();
      for (const expected of EXPECTED_TOOLTIP_ITEMS)
        expect(res.labels, `${v.name}: missing "${expected}"`).toContain(expected);
    });
  }

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((s) => `  - ${s.step}: ${s.error}`).join('\n'));
});
