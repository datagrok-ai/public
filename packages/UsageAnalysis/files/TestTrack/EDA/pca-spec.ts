import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('PCA scenario', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/cars.csv');
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Step 1 — Open cars.csv from Demo files', async () => {
    const info = await page.evaluate(() => {
      const g: any = (window as any).grok;
      return {rows: g.shell.tv.dataFrame.rowCount, cols: g.shell.tv.dataFrame.columns.length};
    });
    expect(info.rows).toBe(30);
    expect(info.cols).toBe(17);
  });

  await softStep('Step 2 — Top Menu > ML > Analyze > PCA...', async () => {
    await page.evaluate(() => {
      const el = document.querySelector('[name="div-ML---Analyze---PCA..."]') as HTMLElement | null;
      if (el) el.click();
    });
    await page.waitForTimeout(1200);
    await page.locator('.d4-dialog').waitFor({timeout: 10000});
    await expect(page.locator('.d4-dialog')).toContainText('PCA');
    // Expected inputs present
    await expect(page.locator('[name="input-host-Table"]')).toBeVisible();
    await expect(page.locator('[name="input-host-Features"]')).toBeVisible();
    await expect(page.locator('[name="input-host-Components"]')).toBeVisible();
    await expect(page.locator('[name="input-host-Center"]')).toBeVisible();
    await expect(page.locator('[name="input-host-Scale"]')).toBeVisible();
  });

  await softStep('Step 3 — Select all Features, set Components = 3', async () => {
    // Open "Select columns..." dialog from Features input
    await page.locator('[name="input-host-Features"] .ui-input-editor').first().click();
    await page.waitForTimeout(700);
    // Click "All" label inside Select columns dialog
    await page.locator('[name="label-All"]').click();
    await page.waitForTimeout(200);
    // Confirm selection — click OK in the Select columns dialog
    // There may be two dialogs stacked; click the last (top-most) OK.
    await page.locator('[name="button-OK"]').last().click();
    await page.waitForTimeout(600);

    // Set Components = 3 on the PCA dialog
    const compInput = page.locator('[name="input-Components"]').first();
    await compInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('3');
    const compVal = await page.evaluate(() => {
      const el = document.querySelector('[name="input-Components"]') as HTMLInputElement | null;
      if (!el) return null;
      if ('value' in el) return (el as HTMLInputElement).value;
      const inner = el.querySelector('input') as HTMLInputElement | null;
      return inner?.value ?? null;
    });
    expect(compVal).toBe('3');
  });

  await softStep('Step 4 — Click OK; expect PC1/PC2/PC3 columns added', async () => {
    // Click PCA dialog OK
    await page.locator('.d4-dialog [name="button-OK"]').last().click();

    // Poll up to 120s for PC1/PC2/PC3 columns to appear
    const deadline = Date.now() + 120_000;
    let found: string[] = [];
    let totalCols = 0;
    while (Date.now() < deadline) {
      const snap = await page.evaluate(() => {
        const g: any = (window as any).grok;
        const df = g.shell.tv?.dataFrame;
        if (!df) return {cols: [], total: 0};
        const names: string[] = [];
        for (let i = 0; i < df.columns.length; i++)
          names.push(df.columns.byIndex(i).name);
        return {cols: names, total: df.columns.length};
      });
      totalCols = snap.total;
      found = snap.cols.filter((n: string) => /^PC\d+$/.test(n));
      if (found.includes('PC1') && found.includes('PC2') && found.includes('PC3'))
        break;
      await page.waitForTimeout(1000);
    }
    if (!(found.includes('PC1') && found.includes('PC2') && found.includes('PC3')))
      throw new Error(`PCA produced no PC1/PC2/PC3 columns within 120s. Total columns: ${totalCols}. PC-like columns found: ${JSON.stringify(found)}.`);
    expect(found).toEqual(expect.arrayContaining(['PC1', 'PC2', 'PC3']));
  });

  await softStep('Step 5 — Repeat with Center and Scale checked', async () => {
    // Re-open PCA dialog via menu
    await page.evaluate(() => {
      const el = document.querySelector('[name="div-ML---Analyze---PCA..."]') as HTMLElement | null;
      if (el) el.click();
    });
    await page.waitForTimeout(1200);
    await page.locator('.d4-dialog').waitFor({timeout: 10000});
    await expect(page.locator('.d4-dialog')).toContainText('PCA');

    // Features: select All
    await page.locator('[name="input-host-Features"] .ui-input-editor').first().click();
    await page.waitForTimeout(700);
    await page.locator('[name="label-All"]').click();
    await page.waitForTimeout(200);
    await page.locator('[name="button-OK"]').last().click();
    await page.waitForTimeout(600);

    // Components = 3
    const compInput = page.locator('[name="input-Components"]');
    await compInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('3');

    // Check Center and Scale checkboxes. Try direct input first; fall back to input
    // nested inside the host container. Dispatch change to ensure Dart listener fires.
    await page.evaluate(() => {
      for (const name of ['Center', 'Scale']) {
        const host = document.querySelector(`[name="input-host-${name}"]`) as HTMLElement | null;
        if (!host) continue;
        const cb = (host.querySelector('input[type="checkbox"]') as HTMLInputElement | null)
          ?? (document.querySelector(`[name="input-${name}"]`) as HTMLInputElement | null);
        if (cb && !cb.checked) {
          cb.click();
        }
      }
    });
    await page.waitForTimeout(200);

    // Click OK
    await page.locator('.d4-dialog [name="button-OK"]').last().click();

    // Poll up to 120s for a *new* PC batch to be added (PC1 (2) etc. or additional PC columns)
    const deadline = Date.now() + 120_000;
    let pcCols: string[] = [];
    let totalCols = 0;
    while (Date.now() < deadline) {
      const snap = await page.evaluate(() => {
        const g: any = (window as any).grok;
        const df = g.shell.tv?.dataFrame;
        if (!df) return {cols: [], total: 0};
        const names: string[] = [];
        for (let i = 0; i < df.columns.length; i++)
          names.push(df.columns.byIndex(i).name);
        return {cols: names, total: df.columns.length};
      });
      totalCols = snap.total;
      pcCols = snap.cols.filter((n: string) => /^PC\d+( \(\d+\))?$/.test(n));
      if (pcCols.length >= 6) break;
      await page.waitForTimeout(1000);
    }
    if (pcCols.length < 6)
      throw new Error(`Repeated PCA with Center/Scale produced insufficient PC columns. Total columns: ${totalCols}. PC columns found: ${JSON.stringify(pcCols)}.`);
    expect(pcCols.length).toBeGreaterThanOrEqual(6);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
