import {test, expect} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/chem/smiles.csv';
const scaffoldsDatasetPath = 'System:AppData/Chem/chembl-scaffolds.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Chem: Info Panels', async ({page}) => {
  // Phase 1: Navigate
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

  // Step 1: Click canonical_smiles column header
  await softStep('Click canonical_smiles column header', async () => {
    await page.evaluate(() => {
      grok.shell.tv.grid.col('canonical_smiles').selected = true;
    });
    await page.waitForTimeout(1000);
  });

  // Step 2: Check Chemistry section in Context Pane
  await softStep('Check Chemistry section for column', async () => {
    const headers = await page.evaluate(() => {
      const headers = document.querySelectorAll('.d4-accordion-pane-header');
      return Array.from(headers).map(h => h.textContent!.trim());
    });
    expect(headers).toContain('Chemistry');
  });

  // Step 3: Expand Chemistry > Rendering
  await softStep('Check Chemistry > Rendering panel', async () => {
    await page.evaluate(() => {
      const headers = document.querySelectorAll('.d4-accordion-pane-header');
      const chemHeader = Array.from(headers).find(h => h.textContent!.trim() === 'Chemistry');
      chemHeader?.click();
      const renderingHeader = Array.from(headers).find(h => h.textContent!.trim() === 'Rendering');
      renderingHeader?.scrollIntoView({behavior: 'instant', block: 'start'});
    });
    await page.waitForTimeout(500);
    const hasScaffoldDropdown = await page.evaluate(() => {
      return !!Array.from(document.querySelectorAll('select'))
        .find(s => Array.from(s.options).some(o => o.text === 'None'));
    });
    expect(hasScaffoldDropdown).toBe(true);
  });

  // Step 4: Click first molecule cell and check info panels
  await softStep('Click first molecule cell and check info panels', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.dataFrame.currentRowIdx = 0;
      grok.shell.tv.dataFrame.currentCol = grok.shell.tv.dataFrame.col('canonical_smiles');
      await new Promise(r => setTimeout(r, 2000));
    });
    const panels = await page.evaluate(() => {
      const headers = document.querySelectorAll('.d4-accordion-pane-header');
      return Array.from(headers).map(h => h.textContent!.trim());
    });
    expect(panels.some(p => /^\d+$/.test(p) || p === 'Links' || p === 'MolregnoInfo')).toBe(true);
  });

  // Step 5: Expand Links and MolregnoInfo tabs
  await softStep('Expand Links and MolregnoInfo tabs', async () => {
    const content = await page.evaluate(async () => {
      const headers = document.querySelectorAll('.d4-accordion-pane-header');
      for (const h of headers) {
        const text = h.textContent!.trim();
        if (text === 'Links' || text === 'MolregnoInfo')
          h.click();
      }
      await new Promise(r => setTimeout(r, 2000));
      const panes = document.querySelectorAll('.d4-accordion-pane');
      const result: Record<string, string> = {};
      for (const p of panes) {
        const header = p.querySelector('.d4-accordion-pane-header');
        const name = header?.textContent?.trim();
        if (name === 'Links' || name === 'MolregnoInfo')
          result[name] = p.innerText.substring(0, 200);
      }
      return result;
    });
    expect(Object.keys(content).length).toBeGreaterThan(0);
  });

  // Step 6: Open chembl-scaffolds.csv and test Rendering with Scaffold column
  await softStep('Open chembl-scaffolds.csv and set Scaffold column', async () => {
    await page.evaluate(async (path) => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
      // Select Smiles column
      grok.shell.tv.grid.col('Smiles').selected = true;
    }, scaffoldsDatasetPath);
    await page.waitForTimeout(1000);

    // Set Scaffold column in Rendering
    await page.evaluate(() => {
      const selects = document.querySelectorAll('select');
      for (const s of selects) {
        if (Array.from(s.options).some(o => o.text === 'Scaffold')) {
          s.value = 'Scaffold';
          s.dispatchEvent(new Event('change', {bubbles: true}));
          break;
        }
      }
    });
    await page.waitForTimeout(1000);
  });

  // Step 7: Check Highlight scaffold
  await softStep('Check Highlight scaffold checkbox', async () => {
    await page.evaluate(() => {
      const checkboxes = document.querySelectorAll('input[type="checkbox"]');
      for (const cb of checkboxes) {
        const text = cb.closest('div')?.parentElement?.innerText || '';
        if (text.includes('Highlight scaffold') && !(cb as HTMLInputElement).checked)
          (cb as HTMLInputElement).click();
      }
    });
    await page.waitForTimeout(1000);
  });

  // Step 8: Add benzene to Highlight
  await softStep('Add benzene highlight structure', async () => {
    // Click Sketch in Highlight section
    await page.evaluate(() => {
      const headers = document.querySelectorAll('.d4-accordion-pane-header');
      const highlightHeader = Array.from(headers).find(h => h.textContent!.trim() === 'Highlight');
      const pane = highlightHeader?.closest('.d4-accordion-pane');
      const sketchLink = pane?.querySelector('.sketch-link, [class*="sketch"]');
      if (sketchLink) (sketchLink as HTMLElement).click();
      else {
        // Find "Sketch" text inside Highlight pane
        const texts = pane?.querySelectorAll('*');
        if (texts) {
          for (const t of texts) {
            if (t.textContent?.trim() === 'Sketch' && t.children.length === 0) {
              (t as HTMLElement).click();
              break;
            }
          }
        }
      }
    });
    await page.waitForTimeout(1000);

    // Type SMILES in the sketcher dialog
    const smilesInput = page.locator('input[placeholder*="SMILES"]');
    await smilesInput.fill('C1CCCCC1');
    await smilesInput.press('Enter');
    await page.waitForTimeout(1000);

    // Click OK
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(2000);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
