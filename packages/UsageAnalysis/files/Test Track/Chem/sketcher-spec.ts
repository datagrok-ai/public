import {test, expect} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/chem/smiles.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Chem: Sketcher', async ({page}) => {
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell && document.querySelector('.d4-root'), {timeout: 30000});

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch(e) {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Step 1: Double-click molecule cell to open sketcher
  await softStep('Double-click molecule cell to open sketcher', async () => {
    await page.evaluate(async () => {
      const canvases = document.querySelectorAll('[name="viewer-Grid"] canvas');
      const canvas = canvases[2]; // overlay canvas
      const rect = canvas.getBoundingClientRect();
      const x = rect.left + 200;
      const y = rect.top + 70;
      for (let i = 0; i < 2; i++) {
        canvas.dispatchEvent(new PointerEvent('pointerdown', {clientX: x, clientY: y, bubbles: true, pointerId: 1}));
        canvas.dispatchEvent(new MouseEvent('mousedown', {clientX: x, clientY: y, bubbles: true, detail: i + 1}));
        canvas.dispatchEvent(new PointerEvent('pointerup', {clientX: x, clientY: y, bubbles: true, pointerId: 1}));
        canvas.dispatchEvent(new MouseEvent('mouseup', {clientX: x, clientY: y, bubbles: true, detail: i + 1}));
        canvas.dispatchEvent(new MouseEvent('click', {clientX: x, clientY: y, bubbles: true, detail: i + 1}));
        if (i === 0) await new Promise(r => setTimeout(r, 100));
      }
      canvas.dispatchEvent(new MouseEvent('dblclick', {clientX: x, clientY: y, bubbles: true, detail: 2}));
      await new Promise(r => setTimeout(r, 2000));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
  });

  // Step 2: Enter SMILES and press Enter
  await softStep('Enter C1CCCCC1 in SMILES input', async () => {
    const input = page.locator('input[placeholder*="SMILES"]');
    await input.fill('C1CCCCC1');
    await input.press('Enter');
    await page.waitForTimeout(1000);
  });

  // Step 3: Click OK and verify cell value
  await softStep('Click OK and verify cell value', async () => {
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(1000);
    const val = await page.evaluate(() => grok.shell.tv.dataFrame.get('canonical_smiles', 0));
    expect(val).toBe('C1CCCCC1');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
