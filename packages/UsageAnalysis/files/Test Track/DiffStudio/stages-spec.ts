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

test('DiffStudio Stages: Acid Production, Multiaxis, Facet, inputs, tooltips', async () => {
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

  // Step 1: Open DiffStudio, load Acid Production
  await softStep('Step 1: Open DiffStudio and load Acid Production', async () => {
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
      const acid = Array.from(items).find(i => i.textContent?.trim() === 'Acid production');
      if (acid) (acid as HTMLElement).click();
      await new Promise(r => setTimeout(r, 5000));
    });

    const info = await page!.evaluate(() => ({
      name: grok.shell.tv?.name,
      rows: grok.shell.tv?.dataFrame?.rowCount,
      cols: grok.shell.tv?.dataFrame?.columns?.length,
    }));
    expect(info.cols).toBe(6);
    expect(info.rows).toBeGreaterThan(0);
  });

  // Step 2: Check Multiaxis and Facet plots
  await softStep('Step 2: Check Multiaxis and Facet plots', async () => {
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

  // Step 3: Modify inputs, observe real-time changes
  await softStep('Step 3: Modify inputs and observe changes', async () => {
    const result = await page!.evaluate(async () => {
      const stageInput = document.querySelector('input[name="input-1-st-stage"]') as HTMLInputElement;
      if (!stageInput) return {error: 'input not found'};

      const valueBefore = stageInput.value;
      const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      nativeSetter.call(stageInput, '40');
      stageInput.dispatchEvent(new Event('input', {bubbles: true}));
      stageInput.dispatchEvent(new Event('change', {bubbles: true}));
      stageInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      await new Promise(r => setTimeout(r, 3000));

      return {valueBefore, valueAfter: stageInput.value};
    });
    expect(result.valueAfter).toBe('40');
  });

  // Step 4: Check tooltips on input hover
  await softStep('Step 4: Check tooltips on input hover', async () => {
    const tooltips = await page!.evaluate(async () => {
      const labels = document.querySelectorAll('label.ui-label');
      const results: {label: string; tooltip: string}[] = [];

      for (const label of labels) {
        const rect = label.getBoundingClientRect();
        if (rect.width === 0) continue;
        const text = label.textContent?.trim();
        if (!text || text.length > 30) continue;

        label.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: rect.left + 5, clientY: rect.top + 5}));
        label.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise(r => setTimeout(r, 500));

        const tooltip = document.querySelector('.d4-tooltip');
        const tooltipText = tooltip?.textContent?.trim().substring(0, 80) || '';
        if (tooltipText && tooltipText !== 'No inputs selected')
          results.push({label: text.substring(0, 20), tooltip: tooltipText});

        label.dispatchEvent(new MouseEvent('mouseout', {bubbles: true}));
        label.dispatchEvent(new MouseEvent('mouseleave', {bubbles: true}));
        await new Promise(r => setTimeout(r, 200));

        if (results.length >= 3) break;
      }
      return results;
    });
    expect(tooltips.length).toBeGreaterThanOrEqual(1);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
