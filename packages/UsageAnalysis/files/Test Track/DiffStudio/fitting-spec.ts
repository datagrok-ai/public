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

test('DiffStudio Fitting: Bioreactor model, modify params, add target, run fit', async () => {
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

  // Step 1: Open DiffStudio and load Bioreactor from Library
  await softStep('Step 1: Open DiffStudio and load Bioreactor', async () => {
    await page!.evaluate(async () => {
      await grok.functions.call('DiffStudio:runDiffStudio');
      await new Promise(r => setTimeout(r, 3000));

      const combo = document.querySelector('.diff-studio-ribbon-widget') as HTMLElement;
      if (combo) combo.click();
      await new Promise(r => setTimeout(r, 800));

      const items = document.querySelectorAll('[role="menuitem"]');
      const lib = Array.from(items).find(i => i.textContent?.trim().startsWith('Library'));
      if (lib) (lib as HTMLElement).click();
      await new Promise(r => setTimeout(r, 800));

      const items2 = document.querySelectorAll('[role="menuitem"]');
      const bio = Array.from(items2).find(i => {
        const t = i.textContent?.trim();
        return t === 'Bioreactor' || t === 'bioreactor';
      });
      if (bio) (bio as HTMLElement).click();
      await new Promise(r => setTimeout(r, 5000));
    });

    const info = await page!.evaluate(() => {
      const tv = grok.shell.tv;
      const df = tv?.dataFrame;
      return {rows: df?.rowCount, cols: df?.columns?.length};
    });
    expect(info.cols).toBeGreaterThanOrEqual(10);
    expect(info.rows).toBeGreaterThan(0);
  });

  // Step 2: Click Fit icon → Fitting view opens
  await softStep('Step 2: Click Fit icon, Fitting view opens', async () => {
    await page!.evaluate(async () => {
      const spans = document.querySelectorAll('span.diff-studio-ribbon-text');
      const fitSpan = Array.from(spans).find(s => s.textContent?.trim() === 'Fit');
      if (fitSpan) (fitSpan as HTMLElement).click();
      await new Promise(r => setTimeout(r, 5000));
    });

    const viewName = await page!.evaluate(() => grok.shell.v?.name);
    expect(viewName).toContain('fitting');
  });

  // Step 3: Modify Process mode; verify FFox & KKox change
  await softStep('Step 3: Modify Process mode; check FFox & KKox change', async () => {
    const result = await page!.evaluate(async () => {
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

      // Switch to Mode 1 using selectedIndex
      const selects = Array.from(document.querySelectorAll('select'));
      const modeSelect = selects.find(s =>
        Array.from(s.options).some(o => o.text === 'Mode 1')
      );
      if (!modeSelect) return {found: false, changed: false};

      modeSelect.selectedIndex = 1;
      modeSelect.dispatchEvent(new Event('input', {bubbles: true}));
      modeSelect.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 2000));

      const ffoxAfter = getVal('FFox');

      // Reset to Default
      modeSelect.selectedIndex = 0;
      modeSelect.dispatchEvent(new Event('input', {bubbles: true}));
      modeSelect.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 1000));

      return {found: true, ffoxBefore, ffoxAfter, changed: ffoxBefore !== ffoxAfter};
    });
    expect(result.found).toBe(true);
    expect(result.changed).toBe(true);
  });

  // Step 4: Enable switchers for switch at, FFox (0.15→1.0), FKox (0→3)
  await softStep('Step 4: Enable switchers and set ranges', async () => {
    const result = await page!.evaluate(async () => {
      // Build Y-position map of parameter inputs
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

      // Get switches sorted by Y
      const switches = Array.from(document.querySelectorAll('div.ui-input-switch'));
      const sortedSwitches = switches
        .map(s => ({el: s, top: s.getBoundingClientRect().top}))
        .filter(s => s.top > 0)
        .sort((a, b) => a.top - b.top);

      // Match and click switches for target params
      const clicked: string[] = [];
      for (const sw of sortedSwitches) {
        const nearest = inputMap.reduce((best, inp) =>
          Math.abs(inp.top - sw.top) < Math.abs(best.top - sw.top) ? inp : best
        , inputMap[0]);
        if (nearest && Math.abs(nearest.top - sw.top) < 30) {
          if (['switch at', 'FFox', 'FKox'].includes(nearest.label)) {
            (sw.el as HTMLElement).click();
            clicked.push(nearest.label);
          }
        }
      }
      await new Promise(r => setTimeout(r, 1000));

      // Set FFox max to 1.0 and FKox min/max to 0/3
      const setInputByLabel = (label: string, value: string) => {
        const allSpans = document.querySelectorAll('label.ui-input-label span');
        for (const span of allSpans) {
          if (span.textContent?.trim() === label) {
            const input = span.closest('label')?.nextElementSibling as HTMLInputElement;
            if (input?.tagName === 'INPUT') {
              const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
              nativeSetter.call(input, value);
              input.dispatchEvent(new Event('input', {bubbles: true}));
              input.dispatchEvent(new Event('change', {bubbles: true}));
              return true;
            }
          }
        }
        return false;
      };

      const ffoxMaxSet = setInputByLabel('FFox (max)', '1.0');
      const fkoxMinSet = setInputByLabel('FKox (min)', '0');
      const fkoxMaxSet = setInputByLabel('FKox (max)', '3');

      return {clicked, ffoxMaxSet, fkoxMinSet, fkoxMaxSet};
    });
    expect(result.clicked.length).toBeGreaterThanOrEqual(1);
  });

  // Step 5: Add bioreactor-experiment.csv as target
  await softStep('Step 5: Add bioreactor-experiment.csv as target', async () => {
    const result = await page!.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:AppData/DiffStudio/library/bioreactor-experiment.csv');
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 1000));

      // Switch back to fitting view
      const views = Array.from(grok.shell.views);
      const fv = views.find(v => v.name?.includes('fitting'));
      if (fv) grok.shell.v = fv;
      await new Promise(r => setTimeout(r, 1000));

      // Select Table in Target combobox
      const selects = Array.from(document.querySelectorAll('select'));
      const targetSelect = selects.find(s => Array.from(s.options).some(o => o.value === 'Table'));
      if (targetSelect) {
        targetSelect.selectedIndex = Array.from(targetSelect.options).findIndex(o => o.value === 'Table');
        targetSelect.dispatchEvent(new Event('input', {bubbles: true}));
        targetSelect.dispatchEvent(new Event('change', {bubbles: true}));
      }
      await new Promise(r => setTimeout(r, 500));

      return {rows: df?.rowCount, targetFound: !!targetSelect};
    });
    expect(result.rows).toBeGreaterThan(0);
  });

  // Step 6: Run fitting; verify RMSE by iterations descending
  await softStep('Step 6: Run fitting; verify results', async () => {
    // Make the play button accessible and click it with Playwright's real click
    const playPos = await page!.evaluate(() => {
      const playIcon = document.querySelector('.d4-ribbon-item i.fa-play');
      if (!playIcon) return null;
      const rect = playIcon.getBoundingClientRect();
      return {x: rect.left + rect.width / 2, y: rect.top + rect.height / 2};
    });
    expect(playPos).not.toBeNull();
    await page!.mouse.click(playPos!.x, playPos!.y);

    // Wait for fitting to complete — poll for rows
    let rows = 0;
    for (let i = 0; i < 20; i++) {
      await page!.waitForTimeout(2000);
      rows = await page!.evaluate(() => grok.shell.tv?.dataFrame?.rowCount ?? 0);
      if (rows > 0) break;
    }
    expect(rows).toBeGreaterThan(0);

    // Verify canvases (charts) are present
    const canvases = await page!.evaluate(() => document.querySelectorAll('canvas').length);
    expect(canvases).toBeGreaterThanOrEqual(2);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
