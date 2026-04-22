import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e?.message ?? String(e)}); }
}

test('Add New Column — formula, rename propagation, save/reopen', async ({page}) => {
  test.setTimeout(300_000);

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

  // Setup phase + open dataset (step 1: Home-dir file).
  await page.evaluate(async () => {
    const g: any = (window as any);
    document.body.classList.add('selenium');
    g.grok.shell.settings.showFiltersIconsConstantly = true;
    g.grok.shell.windows.simpleMode = true;
    g.grok.shell.closeAll();
    const df = await g.grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    const tv = g.grok.shell.addTableView(df);
    await new Promise<void>(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    (window as any).__df = df;
    (window as any).__tv = tv;
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Step 1: Open Home-dir table (demog.csv)', async () => {
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      return {rows: df.rowCount, hasHeight: df.columns.contains('HEIGHT'), hasWeight: df.columns.contains('WEIGHT')};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.hasHeight).toBe(true);
    expect(info.hasWeight).toBe(true);
  });

  // Step 1 variants (Query/GetTop100/GetAll from Northwind) are not available on localhost — skipped.
  test.skip(false, 'Northwind query variants require a Northwind connection; skipped on localhost');

  async function addCalcColumn(name: string, formula: string) {
    await page.locator('[name="icon-add-new-column"]').click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    // Name input: first text input in dialog
    const nameInput = page.locator('.d4-dialog input[type="text"]').first();
    await nameInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type(name);
    // Formula editor is a CodeMirror contenteditable
    const editor = page.locator('.d4-dialog .cm-content');
    await editor.click();
    await page.keyboard.type(formula);
    await page.waitForTimeout(500);
    await page.locator('.d4-dialog [name="button-OK"], .d4-dialog button:has-text("OK")').first().click();
    await page.waitForFunction((n) => {
      const g: any = (window as any);
      return !document.querySelector('.d4-dialog') && g.grok.shell.tv?.dataFrame?.columns?.contains(n);
    }, name, {timeout: 15000});
  }

  await softStep('Step 2: Click Add new column icon', async () => {
    // Already exercised by addCalcColumn in step 3; verify the icon is present.
    await expect(page.locator('[name="icon-add-new-column"]')).toBeVisible();
  });

  await softStep('Step 3: Add formula column CalcCol1 = ${HEIGHT} + ${WEIGHT}', async () => {
    await addCalcColumn('CalcCol1', '${HEIGHT} + ${WEIGHT}');
    const check = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      return {exists: df.columns.contains('CalcCol1'), v0: df.col('CalcCol1').get(0), h0: df.col('HEIGHT').get(0), w0: df.col('WEIGHT').get(0)};
    });
    expect(check.exists).toBe(true);
    expect(Math.abs(check.v0 - (check.h0 + check.w0))).toBeLessThan(0.01);
  });

  await softStep('Step 4: Add dependent CalcCol2 = ${CalcCol1} * 2', async () => {
    await addCalcColumn('CalcCol2', '${CalcCol1} * 2');
    const check = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv.dataFrame;
      return {exists: df.columns.contains('CalcCol2'), c1: df.col('CalcCol1').get(0), c2: df.col('CalcCol2').get(0)};
    });
    expect(check.exists).toBe(true);
    expect(Math.abs(check.c2 - check.c1 * 2)).toBeLessThan(0.01);
  });

  await softStep('Step 5: Rename column and change value — formulas/values update', async () => {
    const result = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      df.col('HEIGHT').name = 'HT';
      await new Promise(r => setTimeout(r, 500));
      const formulaAfterRename = df.col('CalcCol1').tags.get('formula');
      df.col('HT').set(0, 200.0);
      df.fireValuesChanged();
      await new Promise(r => setTimeout(r, 1500));
      return {
        formulaAfterRename,
        ht0: df.col('HT').get(0),
        w0: df.col('WEIGHT').get(0),
        c1: df.col('CalcCol1').get(0),
        c2: df.col('CalcCol2').get(0),
      };
    });
    expect(result.formulaAfterRename).toBe('${HT} + ${WEIGHT}');
    expect(Math.abs(result.c1 - (result.ht0 + result.w0))).toBeLessThan(0.01);
    expect(Math.abs(result.c2 - result.c1 * 2)).toBeLessThan(0.01);
  });

  await softStep('Step 6: Save project (datasync not available — regular save)', async () => {
    const projName = 'TestAddNewColumn_' + Date.now();
    await page.evaluate((n) => { (window as any).__projName = n; }, projName);
    await page.locator('[name="button-Save"]').click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    const nameInput = page.locator('.d4-dialog input[type="text"]').first();
    await nameInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type(projName);
    await page.locator('.d4-dialog [name="button-OK"], .d4-dialog button:has-text("OK")').first().click();
    await page.waitForFunction(() => !document.querySelector('.d4-dialog'), {timeout: 30000});
    await page.waitForTimeout(2000);
    const found = await page.evaluate(async (n) => {
      const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
      return {id: p?.id ?? null, name: p?.name ?? null};
    }, projName);
    expect(found.id).not.toBeNull();
  });

  await softStep('Step 7: Close All', async () => {
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(1500);
    const tables = await page.evaluate(() => (window as any).grok.shell.tables.map((t: any) => t.name));
    expect(tables.length).toBe(0);
  });

  await softStep('Step 8: Reopen the saved project', async () => {
    const result = await page.evaluate(async () => {
      const n = (window as any).__projName;
      const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
      await p.open();
      for (let i = 0; i < 100; i++) {
        if ((window as any).grok.shell.tables.length > 0) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 1000));
      const df = (window as any).grok.shell.tv?.dataFrame;
      return {
        hasHT: df?.columns?.contains('HT'),
        hasCalc1: df?.columns?.contains('CalcCol1'),
        hasCalc2: df?.columns?.contains('CalcCol2'),
        formula: df?.col('CalcCol1')?.tags?.get('formula') ?? null,
      };
    });
    expect(result.hasHT).toBe(true);
    expect(result.hasCalc1).toBe(true);
    expect(result.hasCalc2).toBe(true);
    expect(result.formula).toBe('${HT} + ${WEIGHT}');
  });

  await softStep('Step 9: Rename column in reopened project — formula updates', async () => {
    const result = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      df.col('HT').name = 'HEIGHT_2';
      await new Promise(r => setTimeout(r, 500));
      return {
        hasHEIGHT_2: df.columns.contains('HEIGHT_2'),
        hasHT: df.columns.contains('HT'),
        formula: df.col('CalcCol1').tags.get('formula'),
      };
    });
    expect(result.hasHEIGHT_2).toBe(true);
    expect(result.hasHT).toBe(false);
    expect(result.formula).toBe('${HEIGHT_2} + ${WEIGHT}');
  });

  await softStep('Step 10: Change values — calc columns recalculate', async () => {
    const result = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      df.col('HEIGHT_2').set(0, 180.0);
      df.fireValuesChanged();
      await new Promise(r => setTimeout(r, 1500));
      return {
        h0: df.col('HEIGHT_2').get(0),
        w0: df.col('WEIGHT').get(0),
        c1: df.col('CalcCol1').get(0),
        c2: df.col('CalcCol2').get(0),
      };
    });
    expect(Math.abs(result.c1 - (result.h0 + result.w0))).toBeLessThan(0.01);
    expect(Math.abs(result.c2 - result.c1 * 2)).toBeLessThan(0.01);
  });

  // Cleanup: delete the test project
  await page.evaluate(async () => {
    try {
      const n = (window as any).__projName;
      const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
      if (p) await (window as any).grok.dapi.projects.delete(p);
    } catch (e) { /* ignore cleanup errors */ }
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `- [${e.step}] ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
