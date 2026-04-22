import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e?.message ?? String(e)}); }
}

test('Filter panel — Chem and Bio', async ({page}) => {
  test.setTimeout(600_000);

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  // === Part 1: Refresh with chem filter ===

  await softStep('Open spgi-100 and filter panel', async () => {
    // Phase 2: Open dataset — do NOT open filters here
    await page.evaluate(async () => {
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = true;
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 5000);
      });
      // Wait for Bio/Chem cell rendering + package filter registration
      const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
        .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
      if (hasBioChem) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise(r => setTimeout(r, 200));
        }
        await new Promise(r => setTimeout(r, 5000));
      }
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

    // Phase 3: Open filters after Grid is confirmed
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});
  });

  await softStep('Draw CCC(N(C)C)=O in Structure filter', async () => {
    await page.locator('[name="viewer-Filters"] .sketch-link').first().click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    const smilesInput = page.locator('.d4-dialog input[placeholder*="SMILES"]');
    await smilesInput.fill('CCC(N(C)C)=O');
    await smilesInput.press('Enter');
    await page.waitForTimeout(1000);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForTimeout(2000);

    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(14);
  });

  await softStep('Switch search types and verify counts', async () => {
    await page.evaluate(async () => {
      const filterViewer = document.querySelector('[name="viewer-Filters"]');
      const gearIcon = filterViewer?.querySelector('.chem-search-options-icon') as HTMLElement;
      gearIcon?.click();
      await new Promise(r => setTimeout(r, 500));
    });

    const results = await page.evaluate(async () => {
      const filterViewer = document.querySelector('[name="viewer-Filters"]');
      const cards = filterViewer?.querySelectorAll('.d4-filter');
      let structureCard: Element | null = null;
      cards?.forEach(card => {
        const header = card.querySelector('.d4-filter-header');
        if (header?.textContent?.trim() === 'Structure') structureCard = card;
      });
      const select = structureCard?.querySelector('select') as HTMLSelectElement;
      if (!select) return null;

      const df = grok.shell.tv.dataFrame;
      const results: Record<string, number> = {};
      const types = ['Included in', 'Exact', 'Similar', 'Not contains', 'Not included in', 'Contains'];
      for (const type of types) {
        select.value = type;
        select.dispatchEvent(new Event('input', {bubbles: true}));
        select.dispatchEvent(new Event('change', {bubbles: true}));
        await new Promise(r => setTimeout(r, 3000));
        results[type] = df.filter.trueCount;
      }
      return results;
    });

    expect(results).not.toBeNull();
    expect(results!['Included in']).toBe(0);
    expect(results!['Exact']).toBe(0);
    expect(results!['Similar']).toBe(0);
    expect(results!['Not contains']).toBe(86);
    expect(results!['Not included in']).toBe(100);
    expect(results!['Contains']).toBe(14);
  });

  await softStep('Close All after chem test', async () => {
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(500);
  });

  // === Part 2: Bio ===

  await softStep('Open peptides.csv and open filter panel', async () => {
    // Phase 2: Open dataset — do NOT open filters here
    await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:DemoFiles/bio/peptides.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 5000);
      });
      // Wait for Bio cell rendering + package filter registration
      const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
        .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
      if (hasBioChem) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise(r => setTimeout(r, 200));
        }
        await new Promise(r => setTimeout(r, 5000));
      }
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

    // Phase 3: Open filters — close/reopen if Bio filters missing
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 1000));
      if (document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length === 0) {
        tv.getFiltersGroup().close();
        await new Promise(r => setTimeout(r, 1000));
        tv.getFiltersGroup();
        await new Promise(r => setTimeout(r, 2000));
      }
    });
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});
  });

  await softStep('Enter T-T-Y-K-N-Y-V substructure', async () => {
    const count = await page.evaluate(async () => {
      const filterViewer = document.querySelector('[name="viewer-Filters"]');
      const cards = filterViewer?.querySelectorAll('.d4-filter');
      let seqCard: Element | null = null;
      cards?.forEach(card => {
        const header = card.querySelector('.d4-filter-header');
        if (header?.textContent?.trim() === 'AlignedSequence') seqCard = card;
      });
      const subInput = seqCard?.querySelector('input[placeholder="Substructure"]') as HTMLInputElement;
      if (subInput) {
        subInput.focus();
        subInput.value = 'T-T-Y-K-N-Y-V';
        subInput.dispatchEvent(new Event('input', {bubbles: true}));
        subInput.dispatchEvent(new Event('change', {bubbles: true}));
        subInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
        await new Promise(r => setTimeout(r, 3000));
      }
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(count).toBe(28);
  });

  await softStep('Close All after bio test', async () => {
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(500);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
