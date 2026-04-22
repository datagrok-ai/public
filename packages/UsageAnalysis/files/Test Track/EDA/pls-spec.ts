import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e.message ?? String(e)}); }
}

test('PLS scenario', async ({page}) => {
  test.setTimeout(300_000);

  const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
  const login = process.env.DATAGROK_LOGIN ?? 'admin';
  const password = process.env.DATAGROK_PASSWORD ?? 'admin';

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 60000});

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

  await softStep('Step 2 — Top Menu > ML > Analyze > PLS...', async () => {
    await page.evaluate(() => {
      const el = document.querySelector('[name="div-ML---Analyze---PLS..."]') as HTMLElement | null;
      if (el) el.click();
    });
    await page.waitForTimeout(1200);
    await page.locator('.d4-dialog').waitFor({timeout: 10000});
    await expect(page.locator('.d4-dialog')).toContainText('PLS');
    await expect(page.locator('[name="input-host-Predict"]')).toBeVisible();
    await expect(page.locator('[name="input-host-Using"]')).toBeVisible();
    await expect(page.locator('[name="input-host-Components"]')).toBeVisible();
    await expect(page.locator('[name="input-host-Quadratic"]')).toBeVisible();
  });

  await softStep('Step 3 — Select all features in Using; set Components = 3', async () => {
    // Open "Select columns..." dialog from Using input
    await page.locator('[name="input-host-Using"] .ui-input-editor').first().click();
    await page.waitForTimeout(700);
    await page.locator('[name="label-All"]').click();
    await page.waitForTimeout(200);
    // Confirm selection — click the top-most OK (Select columns sub-dialog)
    await page.locator('[name="button-OK"]').last().click();
    await page.waitForTimeout(600);

    // Set Components = 3 on the PLS dialog
    const compInput = page.locator('[name="input-Components"]').first();
    await compInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('3');
    await page.keyboard.press('Tab');
    const compVal = await page.evaluate(() => {
      const el = document.querySelector('[name="input-Components"]') as HTMLInputElement | null;
      if (!el) return null;
      if ('value' in el) return (el as HTMLInputElement).value;
      const inner = el.querySelector('input') as HTMLInputElement | null;
      return inner?.value ?? null;
    });
    expect(compVal).toBe('3');
  });

  await softStep('Step 4 — Click RUN; expect PLS1/PLS2/PLS3 columns added', async () => {
    // Assert RUN button is enabled. If it's disabled, fail loudly with the reason.
    const runBtn = page.locator('.d4-dialog [name="button-RUN"]').last();
    await runBtn.waitFor({timeout: 10000});
    const disabledAttr = await runBtn.getAttribute('disabled');
    const ariaDisabled = await runBtn.getAttribute('aria-disabled');
    const isDisabled = disabledAttr !== null || ariaDisabled === 'true';
    if (isDisabled) {
      // Capture validation styles on Predict / Using hosts for diagnostics
      const diag = await page.evaluate(() => {
        const names = ['Predict', 'Using', 'Components'];
        const out: Record<string, string> = {};
        for (const n of names) {
          const host = document.querySelector(`[name="input-host-${n}"]`) as HTMLElement | null;
          if (!host) { out[n] = '<missing>'; continue; }
          const sel = host.querySelector('.d4-column-selector, .ui-input-editor') as HTMLElement | null;
          out[n] = sel ? (sel.getAttribute('style') ?? '') : '';
        }
        return out;
      });
      throw new Error(
        'RUN button is disabled — cannot execute PLS. Scenario wording "select all available columns" conflicts '
        + 'with PLS validation (Predict must be excluded from Using). No tooltip explains the disabled state. '
        + `Input styles: ${JSON.stringify(diag)}`);
    }

    await runBtn.click();

    // Poll up to 120s for PLS1/PLS2/PLS3 columns to appear
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
      found = snap.cols.filter((n: string) => /^PLS\d+$/.test(n));
      if (found.includes('PLS1') && found.includes('PLS2') && found.includes('PLS3'))
        break;
      await page.waitForTimeout(1000);
    }
    if (!(found.includes('PLS1') && found.includes('PLS2') && found.includes('PLS3')))
      throw new Error(`PLS produced no PLS1/PLS2/PLS3 columns within 120s. Total columns: ${totalCols}. PLS-like columns found: ${JSON.stringify(found)}.`);
    expect(found).toEqual(expect.arrayContaining(['PLS1', 'PLS2', 'PLS3']));
  });

  await softStep('Step 4b (fallback) — Deselect `price` from Using and retry RUN', async () => {
    // Re-open PLS dialog if it was closed by a prior RUN click.
    const dialogVisible = await page.locator('.d4-dialog').isVisible().catch(() => false);
    if (!dialogVisible) {
      await page.evaluate(() => {
        const el = document.querySelector('[name="div-ML---Analyze---PLS..."]') as HTMLElement | null;
        if (el) el.click();
      });
      await page.waitForTimeout(1200);
      await page.locator('.d4-dialog').waitFor({timeout: 10000});
      // Re-select All before toggling off price
      await page.locator('[name="input-host-Using"] .ui-input-editor').first().click();
      await page.waitForTimeout(700);
      await page.locator('[name="label-All"]').click();
      await page.waitForTimeout(200);
      await page.locator('[name="button-OK"]').last().click();
      await page.waitForTimeout(600);
      const compInput = page.locator('[name="input-Components"]').first();
      await compInput.click();
      await page.keyboard.press('Control+A');
      await page.keyboard.type('3');
      await page.keyboard.press('Tab');
    }

    // Re-open Select columns... sub-dialog for Using
    await page.locator('[name="input-host-Using"] .ui-input-editor').first().click();
    await page.waitForTimeout(700);

    // Try to type 'price' into a search field inside the sub-dialog to narrow down the list.
    const typedInSearch = await page.evaluate(() => {
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      const top = dialogs[dialogs.length - 1] as HTMLElement | undefined;
      if (!top) return false;
      const search = top.querySelector('input[type="search"], input[type="text"]') as HTMLInputElement | null;
      if (!search) return false;
      search.focus();
      search.value = 'price';
      search.dispatchEvent(new Event('input', {bubbles: true}));
      search.dispatchEvent(new Event('change', {bubbles: true}));
      return true;
    });
    if (typedInSearch)
      await page.waitForTimeout(300);

    // Toggle 'price' off — click its label in the top-most dialog.
    const toggled = await page.evaluate(() => {
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      const top = dialogs[dialogs.length - 1] as HTMLElement | undefined;
      if (!top) return false;
      const candidates = Array.from(top.querySelectorAll('label, .ui-label, .d4-link-label, div, span')) as HTMLElement[];
      const target = candidates.find((el) => (el.textContent || '').trim() === 'price');
      if (target) { target.click(); return true; }
      return false;
    });
    if (!toggled) {
      const priceLabel = page.locator('[name="label-price"]').first();
      if (await priceLabel.isVisible({timeout: 2000}).catch(() => false))
        await priceLabel.click();
    }
    await page.waitForTimeout(200);
    await page.locator('[name="button-OK"]').last().click();
    await page.waitForTimeout(600);

    // Assert RUN now enabled
    const runBtn = page.locator('.d4-dialog [name="button-RUN"]').last();
    await runBtn.waitFor({timeout: 10000});
    const disabledAttr = await runBtn.getAttribute('disabled');
    const ariaDisabled = await runBtn.getAttribute('aria-disabled');
    if (disabledAttr !== null || ariaDisabled === 'true')
      throw new Error('RUN button still disabled after deselecting `price` from Using.');

    await runBtn.click();

    // Poll up to 120s for PLS1/PLS2/PLS3
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
      found = snap.cols.filter((n: string) => /^PLS\d+$/.test(n));
      if (found.includes('PLS1') && found.includes('PLS2') && found.includes('PLS3'))
        break;
      await page.waitForTimeout(1000);
    }
    if (!(found.includes('PLS1') && found.includes('PLS2') && found.includes('PLS3')))
      throw new Error(`Fallback PLS run produced no PLS1/PLS2/PLS3 columns within 120s. Total columns: ${totalCols}. PLS-like columns found: ${JSON.stringify(found)}.`);
    expect(found).toEqual(expect.arrayContaining(['PLS1', 'PLS2', 'PLS3']));
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
