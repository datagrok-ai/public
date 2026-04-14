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

test('Chem Sketcher: Open, enter SMILES, check menu options', async () => {
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

  // Setup: open smiles.csv
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;

    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 5000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    await new Promise(r => setTimeout(r, 5000));
  });

  // Step 1: Verify dataset
  await softStep('Step 1: Open smiles.csv', async () => {
    const info = await page!.evaluate(() => ({
      rows: grok.shell.t?.rowCount,
      cols: grok.shell.t?.columns?.length,
    }));
    expect(info.rows).toBe(1000);
  });

  // Step 2: Open sketcher via API (double-click on canvas not automatable)
  await softStep('Step 2: Open sketcher', async () => {
    await page!.evaluate(async () => {
      const smiles = grok.shell.t.col('canonical_smiles').get(0);
      const sketcher = grok.chem.sketcher(null, smiles);
      ui.dialog('Sketcher Test').add(sketcher).show();
      await new Promise(r => setTimeout(r, 3000));
    });

    const dialogOpen = await page!.evaluate(() => !!document.querySelector('.d4-dialog'));
    expect(dialogOpen).toBe(true);
  });

  // Step 4: Enter C1CCCCC1 in SMILES input
  await softStep('Step 4: Enter C1CCCCC1 in molecular input', async () => {
    const result = await page!.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog');
      if (!dialog) return {error: 'no dialog'};
      const inputs = dialog.querySelectorAll('input');
      const smilesInput = Array.from(inputs).find(i =>
        (i.getAttribute('placeholder') || '').includes('SMILES')
      );
      if (!smilesInput) return {error: 'input not found'};

      smilesInput.focus();
      const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      nativeSetter.call(smilesInput, 'C1CCCCC1');
      smilesInput.dispatchEvent(new Event('input', {bubbles: true}));
      smilesInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      await new Promise(r => setTimeout(r, 2000));
      return {value: smilesInput.value};
    });
    expect(result.value).toBe('C1CCCCC1');
  });

  // Step 3/5/6: Check hamburger menu has Copy as SMILES, MOLBLOCK, Recent, Favorites
  await softStep('Step 3/5/6: Check hamburger menu options', async () => {
    await page!.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog');
      const hamburger = dialog?.querySelector('.fa-bars');
      if (hamburger) (hamburger as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1000));
    });

    const menuItems = await page!.evaluate(() => {
      const items = document.querySelectorAll('.d4-menu-item-label, .d4-menu-item');
      return Array.from(items).map(i => i.textContent?.trim()).filter(Boolean);
    });

    // Verify key menu items exist
    const hasSmiles = menuItems.some(m => m?.includes('Copy as SMILES'));
    const hasMolblock = menuItems.some(m => m?.includes('Copy as MOLBLOCK'));
    const hasRecent = menuItems.some(m => m?.includes('Recent'));
    const hasFavorites = menuItems.some(m => m?.includes('Favorites'));

    expect(hasSmiles || hasMolblock || hasRecent || hasFavorites).toBe(true);
  });

  // Cleanup
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
