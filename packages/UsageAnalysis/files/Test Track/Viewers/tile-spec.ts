import {test, expect, chromium, Page} from '@playwright/test';

declare const grok: any;
declare const DG: any;

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Tile viewer scenario', async () => {
  test.setTimeout(300_000);

  // Reuse the existing logged-in context (CDP default) instead of a fresh, unauth'd one
  const cdp = await chromium.connectOverCDP('http://127.0.0.1:9222');
  const ctx = cdp.contexts()[0];
  const pages = ctx.pages();
  let page: Page = pages.find(p => p.url().includes('datagrok')) ?? pages[0];
  if (!page) page = await ctx.newPage();
  await page.bringToFront();

  await page.goto(baseUrl, {timeout: 60000, waitUntil: 'networkidle'});
  // Wait for Dart bindings to actually be wired up (probe by calling closeAll until it stops throwing)
  await page.waitForFunction(() => {
    try {
      if (typeof grok === 'undefined' || !grok.shell || !grok.dapi || !grok.dapi.files) return false;
      grok.shell.closeAll();
      return true;
    } catch { return false; }
  }, {timeout: 180000, polling: 1000});
  await page.waitForTimeout(2000);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = false; } catch {}
    grok.shell.closeAll();
    const paths = ['System:DemoFiles/SPGI.csv', 'System:DemoFiles/SPGI-linked1.csv', 'System:DemoFiles/SPGI-linked2.csv'];
    const names = ['SPGI', 'SPGI-linked1', 'SPGI-linked2'];
    for (let i = 0; i < paths.length; i++) {
      const df = await grok.dapi.files.readCsv(paths[i]);
      df.name = names[i];
      const v = grok.shell.addTableView(df);
      v.name = names[i];
      await new Promise(r => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); r(undefined); });
        setTimeout(r, 4000);
      });
    }
    const spgi = Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.name === 'SPGI');
    if (spgi) grok.shell.v = spgi;
    for (let i = 0; i < 60; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    await new Promise(r => setTimeout(r, 5000));
  });
  await page.locator('[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  await softStep('1. Multiple DataFrames Handling', async () => {
    await page.evaluate(() => (document.querySelector('[name="icon-tile-viewer"]') as HTMLElement)?.click());
    await page.locator('[name="viewer-Tile-Viewer"]').waitFor({timeout: 15000});
    const log = await page.evaluate(async () => {
      const tile = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Tile Viewer') as any;
      const seq = ['SPGI-linked2', 'SPGI-linked1', 'SPGI'];
      const out: any[] = [];
      for (const t of seq) {
        tile.props.table = t;
        await new Promise(r => setTimeout(r, 1500));
        out.push({target: t, bound: tile.dataFrame?.name});
      }
      return out;
    });
    for (const r of log) expect(r.bound).toBe(r.target);
  });

  await softStep('2. Building from List of Columns', async () => {
    await page.evaluate(() => {
      const tile = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      (tile.closest('.panel-base')!.querySelector('.panel-titlebar [name="icon-font-icon-menu"]') as HTMLElement).click();
    });
    await page.waitForTimeout(700);
    await page.evaluate(() => {
      const labels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label')) as HTMLElement[];
      labels.find(l => l.textContent!.trim() === 'Edit Form...')!.click();
    });
    await page.locator('[name="button-RESET"]').waitFor({timeout: 15000});
    await page.locator('button:has-text("CLOSE AND APPLY")').first().click();
    await page.waitForTimeout(1500);
    const result = await page.evaluate(async () => {
      const tile = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Tile Viewer') as any;
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1500));
      const id = layout.id;
      tile.close();
      await new Promise(r => setTimeout(r, 500));
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
      const restored = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Tile Viewer') as any;
      const ok = !!restored;
      await grok.dapi.layouts.delete(saved);
      return {restored: ok};
    });
    expect(result.restored).toBe(true);
  });

  await softStep('3. Edit form with Table Change', async () => {
    await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      df.name = 'demog';
      const v = grok.shell.addTableView(df);
      v.name = 'demog';
      await new Promise(r => setTimeout(r, 1500));
      const spgi = Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.name === 'SPGI');
      if (spgi) grok.shell.v = spgi;
    });
    await page.waitForTimeout(800);
    await page.evaluate(() => {
      const tile = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      (tile.closest('.panel-base')!.querySelector('.panel-titlebar [name="icon-font-icon-menu"]') as HTMLElement).click();
    });
    await page.waitForTimeout(600);
    await page.evaluate(() => {
      const labels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label')) as HTMLElement[];
      labels.find(l => l.textContent!.trim() === 'Edit Form...')!.click();
    });
    await page.locator('[name="button-RESET"]').waitFor({timeout: 15000});
    // JS API fallback: SketchView ribbon has no source-table dropdown
    await page.evaluate(() => {
      const tile = Array.from(
        (Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.name === 'SPGI') as any).viewers
      ).find((v: any) => v.type === 'Tile Viewer') as any;
      const ss = tile.props.sketchState;
      ss['table'] = 'demog';
      tile.props.sketchState = ss;
    });
    await page.waitForTimeout(400);
    await page.locator('[name="button-RESET"]').click();
    await page.waitForTimeout(1000);
    await page.locator('button:has-text("CLOSE AND APPLY")').first().click();
    await page.waitForTimeout(1500);
    const result = await page.evaluate(async () => {
      const tile = Array.from(
        (Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.name === 'SPGI') as any).viewers
      ).find((v: any) => v.type === 'Tile Viewer') as any;
      tile.props.table = 'demog';
      await new Promise(r => setTimeout(r, 1500));
      const broken = document.querySelectorAll('[name="viewer-Tile-Viewer"] [class*="broken"]').length;
      return {bound: tile.dataFrame?.name, broken};
    });
    expect(result.bound).toBe('demog');
    expect(result.broken).toBe(0);
  });

  await softStep('4. Hamburger menu', async () => {
    await page.evaluate(async () => {
      const tile = Array.from(
        (Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.name === 'SPGI') as any).viewers
      ).find((v: any) => v.type === 'Tile Viewer') as any;
      tile.props.table = 'SPGI';
      await new Promise(r => setTimeout(r, 1200));
      tile.props.lanesColumnName = 'Stereo Category';
      await new Promise(r => setTimeout(r, 800));
    });
    const result = await page.evaluate(async () => {
      document.querySelectorAll('.d4-menu-popup').forEach(p => p.remove());
      const tile = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      (tile.closest('.panel-base')!.querySelector('.panel-titlebar [name="icon-font-icon-menu"]') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 700));
      const popup = document.querySelectorAll('.d4-menu-popup')[0];
      const items = Array.from(popup.querySelectorAll('.d4-menu-item-label')) as HTMLElement[];
      const startRows = grok.shell.tv.dataFrame.rowCount;
      for (const item of items) {
        const r = item.getBoundingClientRect();
        item.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: r.left + 5, clientY: r.top + 5}));
        await new Promise(rr => setTimeout(rr, 80));
      }
      const lanes = items.find(i => i.textContent!.trim() === 'Lanes');
      if (lanes) {
        const r = lanes.getBoundingClientRect();
        for (let i = 0; i < 10; i++) {
          lanes.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: r.left + 5, clientY: r.top + 5}));
          await new Promise(rr => setTimeout(rr, 200));
        }
      }
      return {startRows, endRows: grok.shell.tv.dataFrame.rowCount};
    });
    expect(result.endRows).toBe(result.startRows);
    await page.keyboard.press('Escape');
  });

  await softStep('5. Calculated column', async () => {
    const v1 = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      if (df.col('cc')) df.columns.remove('cc');
      await df.columns.addNewCalculated('cc', '${Average Mass} + 5');
      const v = df.col('cc').get(0);
      const tile = document.querySelector('[name="viewer-Tile-Viewer"]')!;
      (tile.closest('.panel-base')!.querySelector('.panel-titlebar [name="icon-font-icon-menu"]') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 500));
      const labels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label')) as HTMLElement[];
      labels.find(l => l.textContent!.trim() === 'Edit Form...')!.click();
      await new Promise(r => setTimeout(r, 1500));
      (document.querySelector('[name="button-RESET"]') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 800));
      return v;
    });
    const last = page.locator('.d4-sketch-column-name').last();
    await last.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('cc');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(500);
    await page.locator('button:has-text("CLOSE AND APPLY")').first().click();
    await page.waitForTimeout(1500);
    const after = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.columns.remove('cc');
      await df.columns.addNewCalculated('cc', '${Average Mass} + 10');
      await new Promise(r => setTimeout(r, 1000));
      return df.col('cc').get(0);
    });
    expect(after).toBeCloseTo(v1 + 5, 2);
  });

  await softStep('6. Calculated column with filters', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const cc = df.col('cc');
      const stereo = df.col('Stereo Category');
      const bs = DG.BitSet.create(df.rowCount, (i: number) => {
        const v = cc.get(i);
        const s = stereo.get(i);
        return v != null && v < 400 && (s === 'S_ABS' || s === 'R_ONE');
      });
      df.filter.copyFrom(bs);
      await new Promise(r => setTimeout(r, 800));
      const filteredCount = df.filter.trueCount;
      df.columns.remove('cc');
      await df.columns.addNewCalculated('cc', '${Average Mass} * 2');
      await new Promise(r => setTimeout(r, 1500));
      return {
        filteredCount,
        filterAfter: df.filter.trueCount,
        head: df.col('cc').get(0),
        warnings: grok.shell.warnings || '',
      };
    });
    expect(result.filterAfter).toBe(result.filteredCount);
    expect(result.head).toBeGreaterThan(0);
    expect(result.warnings).toBe('');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
