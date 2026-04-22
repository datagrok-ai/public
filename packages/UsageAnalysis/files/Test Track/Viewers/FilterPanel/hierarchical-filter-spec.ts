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

test('Hierarchical filter', async ({page}) => {
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
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (grok as any).shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Test case — open filter panel
  await page.evaluate(() => (grok.shell as any).tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  await softStep('Steps 3-4: Add hierarchical filter with SEX, RACE, SEVERITY', async () => {
    await page.evaluate(async () => {
      const fg = (grok.shell as any).tv.getFiltersGroup();
      fg.add({type: 'hierarchical'});
      await new Promise(r => setTimeout(r, 1000));
      for (const f of fg.filters) {
        if (f.filterType === 'hierarchical') {
          const dart = f.dart || f;
          (window as any).grok_GridFilterBase_ApplyState(dart, {
            type: 'hierarchical', active: true,
            colNames: ['SEX', 'RACE', 'SEVERITY'],
            allEnabled: true
          });
          break;
        }
      }
    });
    await page.waitForTimeout(1000);
    const headerText = await page.evaluate(() => {
      const headers = document.querySelectorAll('.d4-filter-header');
      for (const h of headers)
        if (h.textContent!.includes('SEVERITY')) return h.textContent!.trim();
      return '';
    });
    expect(headerText).toContain('SEX');
    expect(headerText).toContain('RACE');
    expect(headerText).toContain('SEVERITY');
  });

  await softStep('Step 5: Expand F > Caucasian, see SEVERITY values', async () => {
    await page.locator('[name="tree-expander-F"]').click();
    await page.waitForTimeout(500);
    await page.evaluate(() => {
      const nodes = document.querySelectorAll('.d4-tree-view-node');
      for (const n of nodes) {
        const v = n.querySelector('.d4-hierarchical-filter-caption-value');
        if (v?.textContent?.trim() === 'Caucasian') {
          const exp = n.querySelector('.d4-tree-view-tri') as HTMLElement;
          if (exp) exp.click();
          break;
        }
      }
    });
    await page.waitForTimeout(500);
    const nodes = await page.evaluate(() => {
      const result: string[] = [];
      for (const n of document.querySelectorAll('.d4-tree-view-node')) {
        const v = n.querySelector('.d4-hierarchical-filter-caption-value');
        if (v) result.push(v.textContent!.trim());
      }
      return result;
    });
    expect(nodes).toContain('None');
    expect(nodes).toContain('High');
    expect(nodes).toContain('Low');
    expect(nodes).toContain('Medium');
    expect(nodes).toContain('Critical');
  });

  await softStep('Step 6: Select Caucasian - all children checked', async () => {
    await page.evaluate(() => {
      for (const n of document.querySelectorAll('.d4-tree-view-node')) {
        const v = n.querySelector('.d4-hierarchical-filter-caption-value');
        if (v?.textContent?.trim() === 'Caucasian') {
          (n.querySelector('.d4-hierarchical-filter-checkbox') as HTMLElement)?.click();
          break;
        }
      }
    });
    await page.waitForTimeout(500);
    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBe(2823);
  });

  await softStep('Step 7: Uncheck Low and Medium - indeterminate parent', async () => {
    await page.evaluate(async () => {
      for (const name of ['Low', 'Medium']) {
        for (const n of document.querySelectorAll('.d4-tree-view-node')) {
          const v = n.querySelector('.d4-hierarchical-filter-caption-value');
          if (v?.textContent?.trim() === name) {
            (n.querySelector('.d4-hierarchical-filter-checkbox') as HTMLElement)?.click();
            break;
          }
        }
        await new Promise(r => setTimeout(r, 200));
      }
    });
    await page.waitForTimeout(500);
    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBe(1917);
  });

  await softStep('Step 8: Reorder to SEX / SEVERITY / RACE', async () => {
    await page.evaluate(() => {
      const fg = (grok.shell as any).tv.getFiltersGroup();
      for (const f of fg.filters) {
        if (f.filterType === 'hierarchical') {
          const dart = f.dart || f;
          (window as any).grok_GridFilterBase_ApplyState(dart, {
            type: 'hierarchical', active: true,
            colNames: ['SEX', 'SEVERITY', 'RACE'],
            allEnabled: true
          });
          grok.shell.tv.dataFrame.rows.requestFilter();
          break;
        }
      }
    });
    await page.waitForTimeout(1000);
    const headerText = await page.evaluate(() => {
      const headers = document.querySelectorAll('.d4-filter-header');
      for (const h of headers)
        if (h.textContent!.includes('SEVERITY')) return h.textContent!.trim();
      return '';
    });
    expect(headerText).toBe('SEX / SEVERITY / RACE');
  });

  await softStep('Step 9: Collaborative filtering - None + DIS_POP', async () => {
    await page.locator('[name="tree-expander-F"]').click();
    await page.waitForTimeout(500);
    await page.evaluate(() => {
      for (const n of document.querySelectorAll('.d4-tree-view-node')) {
        const v = n.querySelector('.d4-hierarchical-filter-caption-value');
        if (v?.textContent?.trim() === 'None') {
          (v as HTMLElement).click();
          break;
        }
      }
    });
    await page.waitForTimeout(500);
    let filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBe(1815);

    await page.evaluate(() => {
      const fg = (grok.shell as any).tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'DIS_POP', selected: ['AS', 'Indigestion']});
    });
    await page.waitForTimeout(500);
    filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBe(339);
  });

  await softStep('Steps 10-12: Save layout, close, restore', async () => {
    const result = await page.evaluate(async () => {
      const layout = (grok.shell as any).tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      const fg = (grok.shell as any).tv.getFiltersGroup();
      fg.close();
      await new Promise(r => setTimeout(r, 1000));
      const afterClose = grok.shell.tv.dataFrame.filter.trueCount;

      const saved = await grok.dapi.layouts.find(layoutId);
      (grok.shell as any).tv.loadLayout(saved);
      (grok.shell as any).tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 3000));
      const afterRestore = grok.shell.tv.dataFrame.filter.trueCount;

      await grok.dapi.layouts.delete(saved);
      return {afterClose, afterRestore};
    });
    expect(result.afterClose).toBe(5850);
    expect(result.afterRestore).toBe(339);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
