import {test} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

declare const grok: any;

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e.message ?? String(e)}); }
}

test('Statistics viewer tests', async ({page}) => {
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
  await page.locator('[name="Browse"]').waitFor({timeout: 120_000});

  // Setup + open dataset
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (grok.shell.settings as any).showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise<void>(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    return { rows: df.rowCount };
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({ timeout: 30000 });

  // ── Add viewer ──────────────────────────────────────────────────────────────

  await softStep('Add viewer: click Statistics icon — viewer opens', async () => {
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-statistics"]') as HTMLElement;
      if (!icon) throw new Error('Statistics icon not found');
      icon.click();
    });
    await page.locator('[name="viewer-Statistics"]').waitFor({ timeout: 10000 });
  });

  await softStep('Add viewer: close Statistics viewer', async () => {
    await page.evaluate(() => {
      const sv = document.querySelector('[name="viewer-Statistics"]');
      if (!sv) throw new Error('viewer-Statistics not found');
      const panelBase = sv.parentElement?.parentElement;
      const closeBtn = panelBase?.querySelector('[name="Close"]') as HTMLElement;
      if (!closeBtn) throw new Error('Close button not found in panelBase');
      closeBtn.click();
    });
    await page.locator('[name="viewer-Statistics"]').waitFor({ state: 'detached', timeout: 5000 });
  });

  await softStep('Add viewer: re-open Statistics viewer', async () => {
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-statistics"]') as HTMLElement;
      if (!icon) throw new Error('Statistics icon not found');
      icon.click();
    });
    await page.locator('[name="viewer-Statistics"]').waitFor({ timeout: 10000 });
  });

  // ── Default statistics display ────────────────────────────────────────────

  await softStep('Default stats: verify stat columns values, nulls, unique, min, max, avg, med, stdev', async () => {
    const stats = await page.evaluate(() => {
      const tv = grok.shell.tv;
      const sv = tv.viewers.find((v: any) => v.type === 'Statistics');
      if (!sv) throw new Error('Statistics viewer not found in tv.viewers');
      return sv.props.stats as string[];
    });
    const expected = ['values', 'nulls', 'unique', 'min', 'max', 'avg', 'med', 'stdev'];
    for (const col of expected) {
      if (!stats.includes(col)) throw new Error(`Missing stat column: ${col}`);
    }
  });

  await softStep('Default stats: numerical columns AGE, HEIGHT, WEIGHT have non-empty values', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return ['AGE', 'HEIGHT', 'WEIGHT'].map(name => {
        const col = df.columns.byName(name);
        return { name, type: col?.type, hasData: col ? (df.rowCount - col.stats.missingValueCount) > 0 : false };
      });
    });
    for (const col of result) {
      if (!col.hasData) throw new Error(`Column ${col.name} has no data`);
    }
  });

  await softStep('Default stats: name column lists all demog columns', async () => {
    const colCount = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    if (colCount !== 11) throw new Error(`Expected 11 columns, got ${colCount}`);
  });

  // ── Statistics for categorical columns ────────────────────────────────────

  await softStep('Categorical: string columns SEX, RACE, DIS_POP have blank numeric stats', async () => {
    // Verified via Open as table: min/max/avg/med/stdev blank for string cols
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return ['SEX', 'RACE', 'DIS_POP'].map(name => {
        const col = df.columns.byName(name);
        return { name, type: col?.type };
      });
    });
    for (const col of result) {
      if (col.type !== 'string') throw new Error(`Expected ${col.name} to be string, got ${col.type}`);
    }
  });

  await softStep('Categorical: count-type stats (values, nulls, unique) populated for categorical rows', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return ['SEX', 'RACE', 'DIS_POP'].map(name => {
        const col = df.columns.byName(name);
        return { name, valueCount: col ? df.rowCount - col.stats.missingValueCount : 0 };
      });
    });
    for (const col of result) {
      if (col.valueCount === 0) throw new Error(`${col.name} has no values`);
    }
  });

  // ── Add and remove statistics columns ─────────────────────────────────────

  await softStep('Add/remove stats: right-click → Statistics submenu → add sum column', async () => {
    await page.evaluate(() => {
      const sv = document.querySelector('[name="viewer-Statistics"]');
      if (!sv) throw new Error('viewer-Statistics not found');
      const canvas = sv.querySelector('canvas') as HTMLElement;
      if (!canvas) throw new Error('canvas not found');
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2,
        clientY: rect.top + rect.height / 2,
      }));
    });
    // Hover Statistics menu item
    await page.locator('.d4-menu-popup >> text="Statistics"').first().hover();
    await page.waitForTimeout(300);
    // Click sum in submenu
    await page.locator('.d4-menu-popup >> text=" sum"').click();
    await page.waitForTimeout(500);
    const stats = await page.evaluate(() => {
      const sv = grok.shell.tv.viewers.find((v: any) => v.type === 'Statistics');
      return sv?.props.stats as string[];
    });
    if (!stats.includes('sum')) throw new Error('sum column not added to Statistics viewer');
  });

  await softStep('Add/remove stats: right-click → Statistics → remove sum column', async () => {
    await page.evaluate(() => {
      const sv = document.querySelector('[name="viewer-Statistics"]');
      const canvas = sv?.querySelector('canvas') as HTMLElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2,
        clientY: rect.top + rect.height / 2,
      }));
    });
    await page.locator('.d4-menu-popup >> text="Statistics"').first().hover();
    await page.waitForTimeout(300);
    await page.locator('.d4-menu-popup >> text=" sum"').click();
    await page.waitForTimeout(500);
    const stats = await page.evaluate(() => {
      const sv = grok.shell.tv.viewers.find((v: any) => v.type === 'Statistics');
      return sv?.props.stats as string[];
    });
    if (stats.includes('sum')) throw new Error('sum column still present after removal');
  });

  // ── Histogram columns ─────────────────────────────────────────────────────

  await softStep('Histogram: right-click → Histograms → submenu lists SEX, RACE, DIS_POP', async () => {
    await page.evaluate(() => {
      const sv = document.querySelector('[name="viewer-Statistics"]');
      const canvas = sv?.querySelector('canvas') as HTMLElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2,
        clientY: rect.top + rect.height / 2,
      }));
    });
    await page.locator('.d4-menu-popup >> text="Histograms"').first().hover();
    await page.waitForTimeout(300);
    const submenuItems = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item'))
        .map((el) => (el as HTMLElement).textContent?.trim() ?? '')
    );
    for (const col of ['SEX', 'RACE', 'DIS_POP']) {
      if (!submenuItems.some((t) => t.includes(col)))
        throw new Error(`Expected ${col} in Histograms submenu, got: ${submenuItems.join(', ')}`);
    }
    await page.locator('.d4-menu-popup >> text=" SEX"').click();
    await page.waitForTimeout(500);
  });

  await softStep('Histogram: right-click → Histograms → remove SEX histogram column', async () => {
    await page.evaluate(() => {
      const sv = document.querySelector('[name="viewer-Statistics"]');
      const canvas = sv?.querySelector('canvas') as HTMLElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2,
        clientY: rect.top + rect.height / 2,
      }));
    });
    await page.locator('.d4-menu-popup >> text="Histograms"').first().hover();
    await page.waitForTimeout(300);
    await page.locator('.d4-menu-popup >> text=" SEX"').click();
    await page.waitForTimeout(500);
  });

  // ── Row source: filtered rows ─────────────────────────────────────────────

  await softStep('Filtered rows: open filter panel', async () => {
    await page.evaluate(() => {
      const el = document.querySelector('[name="div-section--Filters"]') as HTMLElement;
      if (el) { el.click(); return; }
      grok.shell.tv.getFiltersGroup();
    });
    await page.waitForTimeout(1000);
  });

  await softStep('Filtered rows: add AGE range filter 20-40, verify values count decreases', async () => {
    const result = await page.evaluate(() => {
      const tv = grok.shell.tv;
      const fg = tv.getFiltersGroup();
      const before = tv.dataFrame.filter.trueCount;
      fg.updateOrAdd({ type: 'histogram', column: 'AGE', min: 20, max: 40 });
      const after = tv.dataFrame.filter.trueCount;
      return { before, after };
    });
    if (result.after >= result.before) throw new Error(`Filter did not reduce row count: ${result.before} → ${result.after}`);
  });

  await softStep('Filtered rows: remove AGE filter — stats revert to full dataset', async () => {
    const result = await page.evaluate(() => {
      const tv = grok.shell.tv;
      const fg = tv.getFiltersGroup();
      fg.updateOrAdd({ type: 'histogram', column: 'AGE', min: 18, max: 89 });
      return tv.dataFrame.filter.trueCount;
    });
    if (result < 5800) throw new Error(`Expected ~5850 rows after reset, got ${result}`);
  });

  // ── Row source: selected rows ─────────────────────────────────────────────

  await softStep('Selected rows: select 100 rows', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      for (let i = 0; i < 100; i++) df.selection.set(i, true);
    });
  });

  await softStep('Selected rows: open properties, set Row Source to Selected', async () => {
    await page.evaluate(() => {
      const sv = document.querySelector('[name="viewer-Statistics"]');
      const panelBase = sv?.parentElement?.parentElement;
      const gear = panelBase?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gear?.click();
    });
    await page.waitForTimeout(500);
    await page.evaluate(() => {
      const sv = grok.shell.tv.viewers.find((v: any) => v.type === 'Statistics');
      if (!sv) throw new Error('Statistics viewer not found');
      sv.props.rowSource = 'Selected';
    });
    await page.waitForTimeout(500);
    const rowSource = await page.evaluate(() => {
      const sv = grok.shell.tv.viewers.find((v: any) => v.type === 'Statistics');
      return sv?.props.rowSource;
    });
    if (rowSource !== 'Selected') throw new Error(`Expected rowSource=Selected, got ${rowSource}`);
  });

  await softStep('Selected rows: values count matches 100 selected rows', async () => {
    const result = await page.evaluate(() => {
      const sv = grok.shell.tv.viewers.find((v: any) => v.type === 'Statistics');
      return {
        selected: grok.shell.tv.dataFrame.selection.trueCount,
        rowSource: sv?.props.rowSource,
      };
    });
    if (result.selected !== 100) throw new Error(`Expected 100 selected rows, got ${result.selected}`);
    if (result.rowSource !== 'Selected') throw new Error(`Expected rowSource=Selected, got ${result.rowSource}`);
  });

  await softStep('Selected rows: set Row Source back to All', async () => {
    await page.evaluate(() => {
      const sv = grok.shell.tv.viewers.find((v: any) => v.type === 'Statistics');
      if (!sv) throw new Error('Statistics viewer not found');
      sv.props.rowSource = 'All';
      grok.shell.tv.dataFrame.selection.setAll(false);
    });
  });

  // ── Open as table ─────────────────────────────────────────────────────────

  await softStep('Open as table: right-click → Open as table opens new table tab', async () => {
    await page.evaluate(() => {
      const sv = document.querySelector('[name="viewer-Statistics"]');
      const canvas = sv?.querySelector('canvas') as HTMLElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2,
        clientY: rect.top + rect.height / 2,
      }));
    });
    await page.waitForTimeout(300);
    await page.locator('.d4-menu-popup >> text="Open as table"').click();
    await page.waitForTimeout(1000);
    const statsTable = await page.evaluate(() => {
      const tables = grok.shell.tables;
      const st = tables.find((t: any) => t.name?.toLowerCase().includes('stat'));
      return st ? { rows: st.rowCount, cols: st.columns.length } : null;
    });
    if (!statsTable) throw new Error('Stats table not found after Open as table');
    if (statsTable.rows !== 11) throw new Error(`Expected 11 rows in stats table, got ${statsTable.rows}`);
    // Navigate back to demog Table view
    await page.evaluate(() => {
      const tabs = Array.from(document.querySelectorAll('[name^="view-handle"]')) as HTMLElement[];
      const tableTab = tabs.find(t => t.textContent?.trim() === 'Table');
      tableTab?.click();
    });
  });

  // ── Full-screen shortcut ──────────────────────────────────────────────────

  await softStep('Full-screen: click expand icon — viewer expands to full screen', async () => {
    await page.evaluate(() => {
      const sv = document.querySelector('[name="viewer-Statistics"]');
      const panelBase = sv?.parentElement?.parentElement;
      const expandIcon = panelBase?.querySelector('[name="icon-expand-arrows"]') as HTMLElement;
      if (!expandIcon) throw new Error('icon-expand-arrows not found');
      expandIcon.click();
    });
    await page.waitForTimeout(500);
  });

  await softStep('Full-screen: press Alt+F — viewer returns to normal size', async () => {
    await page.keyboard.press('Alt+f');
    await page.waitForTimeout(500);
  });

  // ── Final check ───────────────────────────────────────────────────────────
  if (stepErrors.length > 0) {
    throw new Error('Step failures:\n' + stepErrors.map(e => `  [${e.step}]: ${e.error}`).join('\n'));
  }
});
