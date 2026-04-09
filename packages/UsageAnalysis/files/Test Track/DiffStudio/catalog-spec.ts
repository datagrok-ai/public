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

test('DiffStudio Catalog: Library, Model Hub, Run, Modify', async () => {
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

  // Setup
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;
  });

  // Step 1: Open DiffStudio and load PK-PD from Library
  await softStep('Open DiffStudio and load PK-PD from Library', async () => {
    await page!.evaluate(async () => {
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

    const info = await page!.evaluate(() => {
      const tv = grok.shell.tv;
      const df = tv?.dataFrame;
      return {name: tv?.name, rows: df?.rowCount, cols: df?.columns?.length};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.cols).toBeGreaterThanOrEqual(5);
  });

  // Step 2: Save to Model Hub
  await softStep('Save to Model Hub', async () => {
    const result = await page!.evaluate(async () => {
      const icon = document.querySelector('.diff-studio-ribbon-save-to-model-catalog-icon') as HTMLElement;
      if (icon) icon.click();
      await new Promise(r => setTimeout(r, 3000));
      return {iconFound: !!icon};
    });
    expect(result.iconFound).toBe(true);
  });

  // Step 3: Access Model Hub
  await softStep('Access Model Hub', async () => {
    const viewName = await page!.evaluate(async () => {
      const v = await grok.functions.call('Compute2:modelCatalog');
      grok.shell.v = v;
      await new Promise(r => setTimeout(r, 3000));
      return grok.shell.v?.name;
    });
    expect(viewName).toBe('Model Hub');
  });

  // Step 4: Refresh catalog
  await softStep('Refresh Model Hub catalog', async () => {
    const hasPKPD = await page!.evaluate(async () => {
      const syncIcons = document.querySelectorAll('i.fa-sync');
      for (const el of syncIcons) {
        const rect = el.getBoundingClientRect();
        if (rect.top > 60 && rect.top < 120 && rect.width > 0) {
          (el as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 3000));
      const labels = document.querySelectorAll('.d4-link-label');
      return Array.from(labels).some(l => l.textContent?.trim() === 'PK-PD');
    });
    expect(hasPKPD).toBe(true);
  });

  // Step 5: Run PK-PD from catalog
  await softStep('Run PK-PD model from Model Hub', async () => {
    await page!.evaluate(async () => {
      // Select PK-PD card
      const labels = document.querySelectorAll('span.d4-link-label');
      for (const el of labels) {
        if (el.textContent?.trim() === 'PK-PD') {
          (el as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 1000));

      // Run via JS API (double-click does not open model)
      const obj = grok.shell.o;
      if (obj && typeof obj.apply === 'function')
        await obj.apply();
    });

    // Wait for the new view to appear and switch to it
    await page!.waitForTimeout(5000);
    const info = await page!.evaluate(() => {
      const views = Array.from(grok.shell.views);
      const tableViews = views.filter(v => v.type === 'TableView');
      const newest = tableViews[tableViews.length - 1];
      if (newest) grok.shell.v = newest;
      const tv = grok.shell.tv;
      const df = tv?.dataFrame;
      return {rows: df?.rowCount, cols: df?.columns?.length, name: tv?.name, viewCount: views.length};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.cols).toBeGreaterThanOrEqual(5);
  });

  // Step 6: Modify inputs and verify results
  await softStep('Modify inputs and verify results update', async () => {
    // Ensure the PK-PD TableView is active and wait for render
    await page!.evaluate(async () => {
      const views = Array.from(grok.shell.views);
      const tableViews = views.filter(v => v.type === 'TableView');
      const newest = tableViews[tableViews.length - 1];
      if (newest) {
        grok.shell.v = newest;
        // Click the view tab handle to bring it to front
        const handle = document.querySelector(`[name="view-handle: ${newest.name}"]`) as HTMLElement;
        if (handle) handle.click();
      }
      await new Promise(r => setTimeout(r, 2000));
    });

    // Find and modify the dose input (uses name="input-dose" attribute)
    const result = await page!.evaluate(async () => {
      const doseInput = document.querySelector('input[name="input-dose"]') as HTMLInputElement;
      if (!doseInput) return {found: false};

      doseInput.focus();
      doseInput.select();
      const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      nativeSetter.call(doseInput, '20000');
      doseInput.dispatchEvent(new Event('input', {bubbles: true}));
      doseInput.dispatchEvent(new Event('change', {bubbles: true}));
      doseInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      await new Promise(r => setTimeout(r, 3000));

      const canvases = document.querySelectorAll('canvas').length;
      return {found: true, doseValue: doseInput.value, canvases};
    });

    expect(result.found).toBe(true);
    expect(result.canvases).toBeGreaterThanOrEqual(2);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
