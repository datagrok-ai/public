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

test('Collaborative Filtering for Linked Tables', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

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

  // Phase 2: Run linking script — open 3 tables and link them
  await softStep('1. Run JS script to link 3 tables', async () => {
    const result = await page.evaluate(async () => {
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = true;
      grok.shell.closeAll();

      const df1 = await grok.data.files.openTable('System:DemoFiles/SPGI.csv');
      const df2 = await grok.data.files.openTable('System:DemoFiles/SPGI-linked1.csv');
      const df3 = await grok.data.files.openTable('System:DemoFiles/SPGI-linked2.csv');

      grok.data.linkTables(df3, df2,
        ['Sample Name', 'link column 1', 'link column 2', 'link column 3'],
        ['Sample Name', 'link column 1', 'link column 2', 'link column 3'],
        [DG.SYNC_TYPE.FILTER_TO_FILTER]);
      grok.data.linkTables(df1, df2, ['Id'], ['Concept Id'], [DG.SYNC_TYPE.SELECTION_TO_FILTER]);

      grok.shell.addTableView(df1);
      grok.shell.addTableView(df2);
      grok.shell.addTableView(df3);

      for (const df of [df1, df2, df3]) {
        await new Promise(resolve => {
          const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(resolve, 5000);
        });
      }

      const hasBioChem = [df1, df2, df3].some(df =>
        Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
          .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule'));
      if (hasBioChem) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise(r => setTimeout(r, 200));
        }
        await new Promise(r => setTimeout(r, 5000));
      }

      return {df1: df1.rowCount, df2: df2.rowCount, df3: df3.rowCount};
    });
    expect(result.df1).toBe(3624);
    expect(result.df2).toBe(3624);
    expect(result.df3).toBe(224);
  });

  // Step 2: Go to SPGI view, select 5 rows on top
  await softStep('2. Select 5 rows in SPGI', async () => {
    const result = await page.evaluate(async () => {
      for (const v of grok.shell.tableViews)
        if (v.dataFrame.name === 'SPGI') { grok.shell.v = v; break; }
      await new Promise(r => setTimeout(r, 500));
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      for (let i = 0; i < 5; i++) df.selection.set(i, true);
      return {view: df.name, selected: df.selection.trueCount};
    });
    expect(result.view).toBe('SPGI');
    expect(result.selected).toBe(5);
  });

  // Step 3: Switch to SPGI-linked1 — should contain 9 filtered rows
  await softStep('3. SPGI-linked1 should have 9 filtered rows', async () => {
    const result = await page.evaluate(async () => {
      for (const v of grok.shell.tableViews)
        if (v.dataFrame.name === 'SPGI-linked1') { grok.shell.v = v; break; }
      await new Promise(r => setTimeout(r, 1000));
      return {view: grok.shell.tv.dataFrame.name, filtered: grok.shell.tv.dataFrame.filter.trueCount};
    });
    expect(result.view).toBe('SPGI-linked1');
    expect(result.filtered).toBe(9);
  });

  // Step 4: Switch to SPGI-linked2 view
  await softStep('4. Switch to SPGI-linked2', async () => {
    const view = await page.evaluate(async () => {
      for (const v of grok.shell.tableViews)
        if (v.dataFrame.name === 'SPGI-linked2') { grok.shell.v = v; break; }
      await new Promise(r => setTimeout(r, 500));
      return grok.shell.tv.dataFrame.name;
    });
    expect(view).toBe('SPGI-linked2');
  });

  // Step 5: Open the Filter Panel
  await softStep('5. Open Filter Panel on SPGI-linked2', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      return {open: !!document.querySelector('[name="viewer-Filters"]')};
    });
    expect(result.open).toBe(true);
  });

  // Step 6: For link column 3, select 'v ii'
  await softStep('6. Filter link column 3 to v ii', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'link column 3', selected: ['v ii']});
      await new Promise(r => setTimeout(r, 1000));
      return {filtered: grok.shell.tv.dataFrame.filter.trueCount};
    });
    expect(result.filtered).toBeGreaterThan(0);
  });

  // Step 7: Switch to SPGI-linked1 — should have 5 filtered rows
  await softStep('7. SPGI-linked1 should have 5 filtered rows', async () => {
    const result = await page.evaluate(async () => {
      for (const v of grok.shell.tableViews)
        if (v.dataFrame.name === 'SPGI-linked1') { grok.shell.v = v; break; }
      await new Promise(r => setTimeout(r, 1000));
      return {filtered: grok.shell.tv.dataFrame.filter.trueCount};
    });
    expect(result.filtered).toBe(5);
  });

  // Step 8: Open Filter Panel on SPGI-linked1
  await softStep('8. Open Filter Panel on SPGI-linked1', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      return {open: !!document.querySelector('[name="viewer-Filters"]')};
    });
    expect(result.open).toBe(true);
  });

  // Step 9: Filter PAMPA Classification to 'inconclusive' — 2 rows
  await softStep('9. PAMPA Classification = inconclusive — 2 rows', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'PAMPA Classification', selected: ['inconclusive']});
      await new Promise(r => setTimeout(r, 1000));
      return {filtered: grok.shell.tv.dataFrame.filter.trueCount};
    });
    expect(result.filtered).toBe(2);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
