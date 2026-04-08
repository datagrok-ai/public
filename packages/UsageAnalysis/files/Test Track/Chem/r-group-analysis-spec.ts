import {test, expect} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/chem/sar_small.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Chem: R-Group Analysis', async ({page}) => {
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell && document.querySelector('.d4-root'), {timeout: 30000});

  // Phase 2: Open dataset
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

  // Step 1: Open R-Groups Analysis dialog
  await softStep('Open R-Groups Analysis dialog', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]');
      chemMenu!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const menuItems = document.querySelectorAll('.d4-menu-item-label');
      const rGroup = Array.from(menuItems).find(m => m.textContent!.trim() === 'R-Groups Analysis...');
      rGroup!.closest('.d4-menu-item')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
  });

  // Step 2: Click MCS
  await softStep('Click MCS button', async () => {
    await page.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog');
      const links = dialog!.querySelectorAll('*');
      for (const l of links) {
        if (l.textContent!.trim() === 'MCS' && l.children.length === 0) {
          (l as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 5000));
    });
  });

  // Step 3: Check Visual analysis and click OK
  await softStep('Check Visual analysis and click OK', async () => {
    await page.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog');
      const inputs = dialog!.querySelectorAll('.ui-input-root');
      for (const inp of inputs) {
        if (inp.innerText.trim().startsWith('Visual analysis')) {
          const cb = inp.querySelector('input[type="checkbox"]') as HTMLInputElement;
          if (cb && !cb.checked) cb.click();
          break;
        }
      }
    });
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(10000);
    const colCount = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    expect(colCount).toBeGreaterThan(27);
  });

  // Step 4: Run again with Replace latest unchecked
  await softStep('Run again with Replace latest unchecked', async () => {
    const colsBefore = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]');
      chemMenu!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const menuItems = document.querySelectorAll('.d4-menu-item-label');
      const rGroup = Array.from(menuItems).find(m => m.textContent!.trim() === 'R-Groups Analysis...');
      rGroup!.closest('.d4-menu-item')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 2000));
      // Click MCS
      const dialog = document.querySelector('.d4-dialog');
      const links = dialog!.querySelectorAll('*');
      for (const l of links) {
        if (l.textContent!.trim() === 'MCS' && l.children.length === 0) {
          (l as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 5000));
      // Uncheck Replace latest
      const inputs = dialog!.querySelectorAll('.ui-input-root');
      for (const inp of inputs) {
        if (inp.innerText.trim().startsWith('Replace latest')) {
          const cb = inp.querySelector('input[type="checkbox"]') as HTMLInputElement;
          if (cb && cb.checked) cb.click();
          break;
        }
      }
    });
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(10000);
    const colsAfter = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    expect(colsAfter).toBeGreaterThan(colsBefore);
  });

  // Step 5: Run again with Replace latest checked
  await softStep('Run again with Replace latest checked', async () => {
    const colsBefore = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]');
      chemMenu!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const menuItems = document.querySelectorAll('.d4-menu-item-label');
      const rGroup = Array.from(menuItems).find(m => m.textContent!.trim() === 'R-Groups Analysis...');
      rGroup!.closest('.d4-menu-item')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 2000));
      // Click MCS
      const dialog = document.querySelector('.d4-dialog');
      const links = dialog!.querySelectorAll('*');
      for (const l of links) {
        if (l.textContent!.trim() === 'MCS' && l.children.length === 0) {
          (l as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 5000));
      // Check Replace latest
      const inputs = dialog!.querySelectorAll('.ui-input-root');
      for (const inp of inputs) {
        if (inp.innerText.trim().startsWith('Replace latest')) {
          const cb = inp.querySelector('input[type="checkbox"]') as HTMLInputElement;
          if (cb && !cb.checked) cb.click();
          break;
        }
      }
    });
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(10000);
    const colsAfter = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    expect(colsAfter).toBeLessThanOrEqual(colsBefore);
  });

  // Step 6: Run without MCS — expect error
  await softStep('Run without MCS — expect No core error', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]');
      chemMenu!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const menuItems = document.querySelectorAll('.d4-menu-item-label');
      const rGroup = Array.from(menuItems).find(m => m.textContent!.trim() === 'R-Groups Analysis...');
      rGroup!.closest('.d4-menu-item')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 2000));
    });
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(3000);
    const balloon = await page.evaluate(() => {
      const balloons = document.querySelectorAll('.d4-balloon, .d4-toast');
      return Array.from(balloons).map(b => b.innerText.trim());
    });
    expect(balloon).toContain('No core was provided');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
