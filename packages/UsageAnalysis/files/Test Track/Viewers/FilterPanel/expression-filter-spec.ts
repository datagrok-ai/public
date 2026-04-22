import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';
const datasetPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e?.message ?? String(e)}); }
}

test('Filter panel — Expression filter', async ({page}) => {
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

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Open filters
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  await softStep('Add expression filter with 5 rules', async () => {
    const count = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'expression',
        gridNames: [
          '${WEIGHT} > 50', '${HEIGHT} < 160', '${SEX} equals F',
          '${RACE} contains an', '${STARTED} after 01/01/1991',
        ],
        gridValues: [true, true, true, true, true],
        mode: 'AND',
      });
      await new Promise(r => setTimeout(r, 1000));
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(count).toBe(288);
  });

  await softStep('Switch to OR then back to AND', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const rules = [
        '${WEIGHT} > 50', '${HEIGHT} < 160', '${SEX} equals F',
        '${RACE} contains an', '${STARTED} after 01/01/1991',
      ];
      fg.updateOrAdd({type: 'expression', gridNames: rules, gridValues: [true, true, true, true, true], mode: 'OR'});
      await new Promise(r => setTimeout(r, 1000));
      const orCount = grok.shell.tv.dataFrame.filter.trueCount;
      fg.updateOrAdd({type: 'expression', gridNames: rules, gridValues: [true, true, true, true, true], mode: 'AND'});
      await new Promise(r => setTimeout(r, 1000));
      const andCount = grok.shell.tv.dataFrame.filter.trueCount;
      return {orCount, andCount};
    });
    expect(result.orCount).toBe(5850);
    expect(result.andCount).toBe(288);
  });

  await softStep('Remove first rule (WEIGHT > 50)', async () => {
    const count = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'expression',
        gridNames: ['${HEIGHT} < 160', '${SEX} equals F', '${RACE} contains an', '${STARTED} after 01/01/1991'],
        gridValues: [true, true, true, true],
        mode: 'AND',
      });
      await new Promise(r => setTimeout(r, 1000));
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(count).toBe(330);
  });

  await softStep('Switch to free-text, add 2 rules', async () => {
    const count = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'expression',
        gridNames: [
          '${HEIGHT} < 160', '${SEX} equals F', '${RACE} contains an', '${STARTED} after 01/01/1991',
          'AGE > 30 and SEX = M', 'AGE < 60 and HEIGHT > 190',
        ],
        gridValues: [true, true, true, true, true, true],
        mode: 'AND',
        expressionMode: 'free-text',
      });
      await new Promise(r => setTimeout(r, 1000));
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(count).toBe(0); // Contradictory: SEX=F AND SEX=M
  });

  await softStep('Uncheck first 4 rules', async () => {
    const count = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'expression',
        gridNames: [
          '${HEIGHT} < 160', '${SEX} equals F', '${RACE} contains an', '${STARTED} after 01/01/1991',
          'AGE > 30 and SEX = M', 'AGE < 60 and HEIGHT > 190',
        ],
        gridValues: [false, false, false, false, true, true],
        mode: 'AND',
        expressionMode: 'free-text',
      });
      await new Promise(r => setTimeout(r, 1000));
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(count).toBe(73);
  });

  await softStep('Save layout, close, apply — verify restored state', async () => {
    const result = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));
      const fg = grok.shell.tv.getFiltersGroup();
      fg.close();
      await new Promise(r => setTimeout(r, 1000));
      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
      const fg2 = grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 1000));
      const states = fg2.getStates(null, 'expression');
      try { await grok.dapi.layouts.delete(saved); } catch (e) {}
      return {
        count: grok.shell.tv.dataFrame.filter.trueCount,
        ruleCount: states?.[0]?.gridNames?.length,
        gridValues: states?.[0]?.gridValues,
      };
    });
    expect(result.count).toBe(73);
    expect(result.ruleCount).toBe(6);
    expect(result.gridValues).toEqual([false, false, false, false, true, true]);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
