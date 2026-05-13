import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Helper: hover off the grid then onto the grid canvas, wait past the 200ms
// tooltip debounce, and return the tooltip text content (or '' when hidden).
async function readTooltipAfterHover(page: import('@playwright/test').Page): Promise<{display: string; text: string}> {
  const overlay = page.locator('[name="viewer-Grid"] canvas').last();
  const box = await overlay.boundingBox();
  if (!box) throw new Error('Grid canvas not visible');
  await page.evaluate(() => (window as any).ui?.tooltip?.hide?.());
  await page.mouse.move(0, 0);
  await page.waitForTimeout(200);
  // Move off the grid first to force mouseLeave, then back onto the center to fire mouseEnter
  await page.mouse.move(box.x - 50, box.y + box.height / 2);
  await page.waitForTimeout(150);
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.waitForTimeout(1200);
  return page.evaluate(() => {
    const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
    if (!tt) return {display: 'missing', text: ''};
    return {display: getComputedStyle(tt).display, text: tt.textContent ?? ''};
  });
}

test('Grid: include visible columns in tooltip', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (grok as any).shell.settings.showFiltersIconsConstantly = true;
    (grok as any).shell.windows.simpleMode = true;
    (grok as any).shell.closeAll();
    const df = await (grok as any).dapi.files.readCsv('System:DemoFiles/energy_uk.csv');
    (grok as any).shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // The scenario's expected tooltip behavior matches `showTooltip = 'inherit from table'`.
  // The grid's actual default is 'show custom tooltip', which short-circuits the tooltip
  // pipeline when `rowTooltip` is empty (regardless of `showVisibleColumnsInTooltip`) —
  // tooltip.dart#getRowTooltip lines 525-537 returns null. Switching mode lets the rest
  // of the scenario observe the documented `showVisibleColumnsInTooltip` filter behavior.
  await page.evaluate(() => {
    const tv: any = (grok as any).shell.tv;
    tv.grid.props.showTooltip = 'inherit from table';
  });

  await softStep('Open grid properties and find Show Visible Columns In Tooltip', async () => {
    await page.evaluate(() => {
      const grid = document.querySelector('[name="viewer-Grid"]')!;
      (grid.querySelector('[name="icon-font-icon-settings"]') as HTMLElement).click();
    });
    await page.locator('[name="prop-category-tooltip"]').waitFor({timeout: 5000});
    await page.evaluate(() => {
      const cat = document.querySelector('[name="prop-category-tooltip"]') as HTMLElement | null;
      const row = document.querySelector('[name="prop-show-visible-columns-in-tooltip"]') as HTMLElement | null;
      if (row && getComputedStyle(row).display === 'none' && cat) cat.click();
    });
    await page.locator('[name="prop-show-visible-columns-in-tooltip"]').waitFor({state: 'attached', timeout: 5000});
    const defaultVal = await page.evaluate(() => (grok as any).shell.tv.grid.props.showVisibleColumnsInTooltip);
    expect(defaultVal).toBe(false);
  });

  await softStep('Unchecked + all columns visible: no tooltip on cell hover', async () => {
    await page.evaluate(() => {
      const tv: any = (grok as any).shell.tv;
      tv.grid.props.showVisibleColumnsInTooltip = false;
      tv.grid.columns.byName('value').width = 162;
      tv.grid.invalidate?.();
    });
    await page.waitForTimeout(400);
    const tt = await readTooltipAfterHover(page);
    expect(tt.display).toBe('none');
    expect(tt.text).toBe('');
  });

  await softStep('Unchecked + last column clipped: tooltip shows the clipped column', async () => {
    await page.evaluate(() => {
      const tv: any = (grok as any).shell.tv;
      tv.grid.columns.byName('value').width = 1100;
      tv.grid.invalidate?.();
    });
    await page.waitForTimeout(400);
    const tt = await readTooltipAfterHover(page);
    expect(tt.display).toBe('block');
    // Only the clipped column ('target') should appear when showVisibleColumnsInTooltip is off.
    expect(tt.text).toContain('target');
    expect(tt.text).not.toContain('source');
  });

  await softStep('Enable Show Visible Columns In Tooltip', async () => {
    await page.evaluate(() => {
      if (!document.querySelector('[name="prop-show-visible-columns-in-tooltip"]')) {
        const grid = document.querySelector('[name="viewer-Grid"]')!;
        (grid.querySelector('[name="icon-font-icon-settings"]') as HTMLElement).click();
      }
    });
    await page.locator('[name="prop-category-tooltip"]').waitFor({timeout: 5000});
    await page.evaluate(() => {
      const cat = document.querySelector('[name="prop-category-tooltip"]') as HTMLElement | null;
      const row = document.querySelector('[name="prop-show-visible-columns-in-tooltip"]') as HTMLElement | null;
      if (row && getComputedStyle(row).display === 'none' && cat) cat.click();
    });
    await page.waitForTimeout(200);
    await page.evaluate(() => {
      const row = document.querySelector('[name="prop-show-visible-columns-in-tooltip"]') as HTMLElement | null;
      const checkbox = row?.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
      if (checkbox) checkbox.click();
      else (grok as any).shell.tv.grid.props.showVisibleColumnsInTooltip = true;
    });
    await page.waitForTimeout(300);
    const propValue = await page.evaluate(() => (grok as any).shell.tv.grid.props.showVisibleColumnsInTooltip);
    expect(propValue).toBe(true);
  });

  await softStep('Checked: tooltip shows the same visible columns in both states', async () => {
    // All columns visible
    await page.evaluate(() => {
      const tv: any = (grok as any).shell.tv;
      tv.grid.columns.byName('value').width = 162;
      tv.grid.invalidate?.();
    });
    await page.waitForTimeout(400);
    const ttAll = await readTooltipAfterHover(page);
    expect(ttAll.display).toBe('block');
    expect(ttAll.text).toContain('value');
    expect(ttAll.text).toContain('source');
    expect(ttAll.text).toContain('target');

    // Last column clipped
    await page.evaluate(() => {
      const tv: any = (grok as any).shell.tv;
      tv.grid.columns.byName('value').width = 1100;
      tv.grid.invalidate?.();
    });
    await page.waitForTimeout(400);
    const ttClipped = await readTooltipAfterHover(page);
    expect(ttClipped.display).toBe('block');
    expect(ttClipped.text).toContain('value');
    expect(ttClipped.text).toContain('source');
    expect(ttClipped.text).toContain('target');
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((s) => `  - ${s.step}: ${s.error}`).join('\n'));
});
