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

test('DiffStudio Sensitivity Analysis: Bioreactor, select params, run Monte Carlo', async () => {
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

  // Step 1: Open DiffStudio, load Bioreactor
  await softStep('Step 1: Open DiffStudio and load Bioreactor', async () => {
    await page!.evaluate(async () => {
      await grok.functions.call('DiffStudio:runDiffStudio');
      await new Promise(r => setTimeout(r, 3000));

      const combo = document.querySelector('.diff-studio-ribbon-widget') as HTMLElement;
      if (combo) combo.click();
      await new Promise(r => setTimeout(r, 800));
      let items = document.querySelectorAll('[role="menuitem"]');
      const lib = Array.from(items).find(i => i.textContent?.trim().startsWith('Library'));
      if (lib) (lib as HTMLElement).click();
      await new Promise(r => setTimeout(r, 800));
      items = document.querySelectorAll('[role="menuitem"]');
      const bio = Array.from(items).find(i => i.textContent?.trim() === 'Bioreactor');
      if (bio) (bio as HTMLElement).click();
      await new Promise(r => setTimeout(r, 5000));
    });

    const info = await page!.evaluate(() => ({
      rows: grok.shell.tv?.dataFrame?.rowCount,
      cols: grok.shell.tv?.dataFrame?.columns?.length,
    }));
    expect(info.cols).toBe(13);
    expect(info.rows).toBe(1001);
  });

  // Step 2: Click Sensitivity icon
  await softStep('Step 2: Open Sensitivity Analysis view', async () => {
    await page!.evaluate(async () => {
      const spans = document.querySelectorAll('span.diff-studio-ribbon-text');
      const sensSpan = Array.from(spans).find(s => s.textContent?.trim() === 'Sensitivity');
      if (sensSpan) (sensSpan as HTMLElement).click();
      await new Promise(r => setTimeout(r, 5000));
    });

    const viewName = await page!.evaluate(() => grok.shell.v?.name);
    expect(viewName).toContain('comparison');
  });

  // Step 3: Modify Process mode, check FFox & KKox change; enable params
  await softStep('Step 3: Modify Process mode and select parameters', async () => {
    // Verify Process mode changes FFox
    const modeResult = await page!.evaluate(async () => {
      const getVal = (label: string) => {
        const spans = document.querySelectorAll('label.ui-input-label span');
        for (const span of spans) {
          if (span.textContent?.trim() === label) {
            const input = span.closest('label')?.nextElementSibling as HTMLInputElement;
            if (input?.tagName === 'INPUT') return input.value;
          }
        }
        return null;
      };

      const ffoxBefore = getVal('FFox');
      const selects = Array.from(document.querySelectorAll('select'));
      const modeSelect = selects.find(s => Array.from(s.options).some(o => o.text === 'Mode 1'));
      if (modeSelect) {
        modeSelect.selectedIndex = 1;
        modeSelect.dispatchEvent(new Event('input', {bubbles: true}));
        modeSelect.dispatchEvent(new Event('change', {bubbles: true}));
        await new Promise(r => setTimeout(r, 2000));
      }
      const ffoxAfter = getVal('FFox');

      // Reset
      if (modeSelect) {
        modeSelect.selectedIndex = 0;
        modeSelect.dispatchEvent(new Event('input', {bubbles: true}));
        modeSelect.dispatchEvent(new Event('change', {bubbles: true}));
        await new Promise(r => setTimeout(r, 1000));
      }
      return {ffoxBefore, ffoxAfter, changed: ffoxBefore !== ffoxAfter};
    });
    expect(modeResult.changed).toBe(true);

    // Enable switches for FFox, FFred, FKox
    const switchResult = await page!.evaluate(async () => {
      const inputMap: {label: string; top: number}[] = [];
      const spans = document.querySelectorAll('label.ui-input-label span');
      for (const span of spans) {
        const label = span.closest('label');
        const input = label?.nextElementSibling as HTMLInputElement;
        if (input?.tagName === 'INPUT') {
          const rect = input.getBoundingClientRect();
          if (rect.width > 0)
            inputMap.push({label: span.textContent?.trim() ?? '', top: rect.top});
        }
      }

      const switches = Array.from(document.querySelectorAll('div.ui-input-switch'));
      const sortedSwitches = switches
        .map(s => ({el: s, top: s.getBoundingClientRect().top, on: s.classList.contains('ui-input-switch-on')}))
        .filter(s => s.top > 0)
        .sort((a, b) => a.top - b.top);

      const targetParams = ['FFox', 'FKox', 'FFred'];
      const clicked: string[] = [];
      for (const sw of sortedSwitches) {
        const nearest = inputMap.reduce((best, inp) =>
          Math.abs(inp.top - sw.top) < Math.abs(best.top - sw.top) ? inp : best
        , inputMap[0]);
        if (nearest && Math.abs(nearest.top - sw.top) < 30) {
          if (targetParams.includes(nearest.label) && !sw.on) {
            (sw.el as HTMLElement).click();
            clicked.push(nearest.label);
          }
        }
      }
      await new Promise(r => setTimeout(r, 1000));
      return {clicked};
    });
    expect(switchResult.clicked.length).toBeGreaterThanOrEqual(1);
  });

  // Step 4: Run sensitivity analysis, verify 4 viewers
  await softStep('Step 4: Run analysis and verify viewers', async () => {
    // Make play button accessible and click via real mouse
    await page!.evaluate(() => {
      const playIcon = document.querySelector('.d4-ribbon-item i.fa-play');
      if (playIcon) {
        const parent = playIcon.closest('.d4-ribbon-item')!;
        parent.setAttribute('role', 'button');
        parent.setAttribute('aria-label', 'Run analysis');
      }
    });

    const playPos = await page!.evaluate(() => {
      const playIcon = document.querySelector('.d4-ribbon-item i.fa-play');
      if (!playIcon) return null;
      const rect = playIcon.getBoundingClientRect();
      return {x: rect.left + rect.width / 2, y: rect.top + rect.height / 2};
    });
    expect(playPos).not.toBeNull();
    await page!.mouse.click(playPos!.x, playPos!.y);

    // Wait for results
    let canvases = 0;
    let rows = 0;
    for (let i = 0; i < 25; i++) {
      await page!.waitForTimeout(2000);
      canvases = await page!.evaluate(() => document.querySelectorAll('canvas').length);
      rows = await page!.evaluate(() => grok.shell.tv?.dataFrame?.rowCount ?? 0);
      if (canvases >= 4 && rows > 0) break;
    }
    expect(rows).toBeGreaterThan(0);
    expect(canvases).toBeGreaterThanOrEqual(4);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
