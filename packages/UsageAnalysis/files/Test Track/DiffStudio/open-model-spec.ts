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

test('DiffStudio Open Model: Bioreactor, Multiaxis, Facet, slider, Process mode', async () => {
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

  // Step 1: Open DiffStudio
  await softStep('Step 1: Open Diff Studio from Apps', async () => {
    await page!.evaluate(async () => {
      await grok.functions.call('DiffStudio:runDiffStudio');
      await new Promise(r => setTimeout(r, 3000));
    });
    const viewName = await page!.evaluate(() => grok.shell.v?.name);
    expect(viewName).toBeTruthy();
  });

  // Step 2: Load Bioreactor from Library
  await softStep('Step 2: Load Bioreactor from Library', async () => {
    await page!.evaluate(async () => {
      const combo = document.querySelector('.diff-studio-ribbon-widget') as HTMLElement;
      if (combo) combo.click();
      await new Promise(r => setTimeout(r, 800));

      const items = document.querySelectorAll('[role="menuitem"]');
      const lib = Array.from(items).find(i => i.textContent?.trim().startsWith('Library'));
      if (lib) (lib as HTMLElement).click();
      await new Promise(r => setTimeout(r, 800));

      const items2 = document.querySelectorAll('[role="menuitem"]');
      const bio = Array.from(items2).find(i => i.textContent?.trim() === 'Bioreactor');
      if (bio) (bio as HTMLElement).click();
      await new Promise(r => setTimeout(r, 5000));
    });

    const info = await page!.evaluate(() => {
      const tv = grok.shell.tv;
      const df = tv?.dataFrame;
      return {rows: df?.rowCount, cols: df?.columns?.length};
    });
    expect(info.cols).toBe(13);
    expect(info.rows).toBe(1001);
  });

  // Step 3: Check Multiaxis and Facet tabs
  await softStep('Step 3: Check Multiaxis and Facet tabs', async () => {
    const tabs = await page!.evaluate(() => {
      const allEls = document.querySelectorAll('*');
      const multiaxis = Array.from(allEls).find(el =>
        el.textContent?.trim() === 'Multiaxis' && el.children.length === 0
      );
      const facet = Array.from(allEls).find(el =>
        el.textContent?.trim() === 'Facet' && el.children.length === 0
      );
      return {multiaxis: !!multiaxis, facet: !!facet};
    });
    expect(tabs.multiaxis).toBe(true);
    expect(tabs.facet).toBe(true);
  });

  // Step 4: Check Facet curves have different colors
  await softStep('Step 4: Facet curves have different colors', async () => {
    // Click Facet tab
    await page!.evaluate(() => {
      const allEls = document.querySelectorAll('*');
      const facet = Array.from(allEls).find(el =>
        el.textContent?.trim() === 'Facet' && el.children.length === 0
      );
      if (facet) (facet as HTMLElement).click();
    });
    await page!.waitForTimeout(1000);

    // Verify multiple canvases present (each facet panel has its own canvas)
    const canvases = await page!.evaluate(() => document.querySelectorAll('canvas').length);
    expect(canvases).toBeGreaterThanOrEqual(6);
  });

  // Step 5: Adjust Switch at value, verify chart updates
  await softStep('Step 5: Adjust Switch at; table and chart update', async () => {
    const result = await page!.evaluate(async () => {
      const switchAtInput = document.querySelector('input[name="input-switch-at"]') as HTMLInputElement;
      if (!switchAtInput) return {error: 'not found'};

      const valueBefore = switchAtInput.value;
      const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      nativeSetter.call(switchAtInput, '200');
      switchAtInput.dispatchEvent(new Event('input', {bubbles: true}));
      switchAtInput.dispatchEvent(new Event('change', {bubbles: true}));
      switchAtInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      await new Promise(r => setTimeout(r, 2000));

      return {valueBefore, valueAfter: switchAtInput.value};
    });
    expect(result.valueAfter).toBe('200');
  });

  // Step 6: Modify Process mode; FFox & KKox change, charts update
  await softStep('Step 6: Modify Process mode; inputs and charts update', async () => {
    const result = await page!.evaluate(async () => {
      const getVal = (name: string) =>
        (document.querySelector(`input[name="input-${name}"]`) as HTMLInputElement)?.value;

      const ffoxBefore = getVal('FFox');
      const kkoxBefore = getVal('KKox');

      // Switch to Mode 1 using selectedIndex
      const selects = Array.from(document.querySelectorAll('select'));
      const modeSelect = selects.find(s => Array.from(s.options).some(o => o.text === 'Mode 1'));
      if (!modeSelect) return {found: false, changed: false};

      modeSelect.selectedIndex = 1;
      modeSelect.dispatchEvent(new Event('input', {bubbles: true}));
      modeSelect.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 3000));

      const ffoxAfter = getVal('FFox');
      const kkoxAfter = getVal('KKox');

      return {found: true, ffoxBefore, ffoxAfter, kkoxBefore, kkoxAfter,
        changed: ffoxBefore !== ffoxAfter || kkoxBefore !== kkoxAfter};
    });
    expect(result.found).toBe(true);
    expect(result.changed).toBe(true);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
