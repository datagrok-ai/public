import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('DiffStudio Cyclic Models: PK-PD clickers and tooltips', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
    await page.waitForFunction(() => {
      try { return typeof grok !== 'undefined' && typeof grok.shell.closeAll === 'function'; }
      catch { return false; }
    }, {timeout: 45000});
  }

  // Setup: launch DiffStudio and load PK-PD from Library
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;
    await grok.functions.call('DiffStudio:runDiffStudio');
    await new Promise(r => setTimeout(r, 3000));

    // Click open model combo
    const combo = document.querySelector('.diff-studio-ribbon-widget') as HTMLElement;
    if (combo) combo.click();
    await new Promise(r => setTimeout(r, 800));

    // Click Library
    const items = document.querySelectorAll('[role="menuitem"]');
    const lib = Array.from(items).find(i => i.textContent?.trim().startsWith('Library'));
    if (lib) (lib as HTMLElement).click();
    await new Promise(r => setTimeout(r, 800));

    // Select PK-PD
    const items2 = document.querySelectorAll('[role="menuitem"]');
    const pkpd = Array.from(items2).find(i => {
      const t = i.textContent?.trim();
      return t === 'PK-PD' || t === 'pk-pd';
    });
    if (pkpd) (pkpd as HTMLElement).click();
    await new Promise(r => setTimeout(r, 5000));
  });

  // Step 1: Verify PK-PD model loaded
  await softStep('Open DiffStudio and load PK-PD from Library', async () => {
    const info = await page!.evaluate(() => {
      const tv = grok.shell.tv;
      const df = tv?.dataFrame;
      return {name: tv?.name, rows: df?.rowCount, cols: df?.columns?.length};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.cols).toBeGreaterThanOrEqual(5);
  });

  // Step 2: Check Multiaxis and Facet tabs
  await softStep('Check Multiaxis and Facet plots are updated', async () => {
    const result = await page!.evaluate(() => {
      const allEls = Array.from(document.querySelectorAll('*'));
      const hasMultiaxis = allEls.some(el => el.textContent?.trim() === 'Multiaxis' && el.children.length === 0);
      const hasFacet = allEls.some(el => el.textContent?.trim() === 'Facet' && el.children.length === 0);
      const canvasCount = document.querySelectorAll('canvas').length;
      const lineCharts = document.querySelectorAll('[name^="viewer-Line-chart"]').length;
      return {hasMultiaxis, hasFacet, canvasCount, lineCharts};
    });
    expect(result.hasMultiaxis).toBe(true);
    expect(result.hasFacet).toBe(true);
    expect(result.canvasCount).toBeGreaterThanOrEqual(2);
    expect(result.lineCharts).toBeGreaterThanOrEqual(4);
  });

  // Step 3: Modify Count using clickers
  await softStep('Modify Count input using clickers', async () => {
    const result = await page!.evaluate(async () => {
      const initialRows = grok.shell.tv?.dataFrame?.rowCount;

      // Find count input and click plus clicker
      const countInput = Array.from(document.querySelectorAll('input')).find(i => {
        const parent = i.closest('.ui-input-root');
        const label = parent?.querySelector('.ui-label');
        return label?.textContent?.trim().toLowerCase() === 'count';
      }) as HTMLInputElement;
      const parent = countInput?.closest('.ui-input-root');
      const plusIcon = parent?.querySelector('.ui-input-plus, [name="icon-plus"]') as HTMLElement;
      if (plusIcon) plusIcon.click();
      await new Promise(r => setTimeout(r, 3000));

      const afterPlusRows = grok.shell.tv?.dataFrame?.rowCount;
      const afterPlusValue = countInput?.value;

      // Click minus to go back
      const minusIcon = parent?.querySelector('.ui-input-minus, [name="icon-minus"]') as HTMLElement;
      if (minusIcon) minusIcon.click();
      await new Promise(r => setTimeout(r, 3000));

      const afterMinusRows = grok.shell.tv?.dataFrame?.rowCount;
      const afterMinusValue = countInput?.value;

      return {initialRows, afterPlusRows, afterPlusValue, afterMinusRows, afterMinusValue};
    });
    // Plus: 10→11 should increase rows (1210 → 1331)
    expect(result.afterPlusRows).toBeGreaterThan(result.initialRows!);
    expect(result.afterPlusValue).toBe('11');
    // Minus: 11→10 should restore rows
    expect(result.afterMinusRows).toBe(result.initialRows);
    expect(result.afterMinusValue).toBe('10');
  });

  // Step 4: Check tooltips on Begin, End, Step inputs
  await softStep('Check tooltips on Begin, End, Step inputs', async () => {
    const tooltips = await page!.evaluate(async () => {
      const results: {name: string; tooltip: string}[] = [];
      const hosts = ['input-host-begin', 'input-host-end', 'input-host-step'];
      for (const hostName of hosts) {
        const host = document.querySelector(`[name="${hostName}"]`);
        const input = host?.querySelector('input');
        if (input) {
          input.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
          input.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
          await new Promise(r => setTimeout(r, 1500));
          const tooltip = document.querySelector('.d4-tooltip');
          results.push({
            name: hostName,
            tooltip: tooltip?.textContent?.trim() || ''
          });
          // Move away to dismiss tooltip
          document.body.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
          await new Promise(r => setTimeout(r, 300));
        }
      }
      return results;
    });

    for (const t of tooltips) {
      expect(t.tooltip.length, `Tooltip for ${t.name} should not be empty`).toBeGreaterThan(0);
    }
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
