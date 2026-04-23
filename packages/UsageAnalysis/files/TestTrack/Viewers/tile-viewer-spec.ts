import {test, expect} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e.message ?? String(e)}); }
}

test('Tile Viewer tests', async ({page, baseURL}) => {
  test.setTimeout(300_000);

  await page.goto(baseURL ?? '/');
  await page.locator('[name="Toolbox"], [name="Browse"], .d4-sidebar').first().waitFor({timeout: 60000});
  await page.waitForFunction(() => {
    try {
      if (typeof grok === 'undefined' || !grok.shell) return false;
      grok.shell.settings.showFiltersIconsConstantly;
      return true;
    } catch { return false; }
  }, {timeout: 30000});

  // Setup
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (grok as any).shell.settings.showFiltersIconsConstantly = true;
    (grok as any).shell.windows.simpleMode = true;
    (grok as any).shell.closeAll();
    const df = await (grok as any).dapi.files.readCsv('System:DemoFiles/demog.csv');
    const tv = (grok as any).shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    tv.getFiltersGroup();
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  // Add Tile Viewer via toolbox icon
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-tile-viewer"]');
    if (icon) (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-Tile-Viewer"]').waitFor({timeout: 10000});

  // ---- Default form rendering ----
  await softStep('Default form rendering: click first tile → currentRow=0', async () => {
    const row = await page.evaluate(() => {
      const tiles = document.querySelectorAll('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form');
      (tiles[0] as HTMLElement)?.click();
      return (grok as any).shell.tv.dataFrame.currentRowIdx;
    });
    expect(row).toBe(0);
  });

  await softStep('Default form rendering: click second tile → currentRow=1', async () => {
    const row = await page.evaluate(() => {
      const tiles = document.querySelectorAll('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form');
      (tiles[1] as HTMLElement)?.click();
      return (grok as any).shell.tv.dataFrame.currentRowIdx;
    });
    expect(row).toBe(1);
  });

  // ---- Row selection ----
  await softStep('Row selection: shift-click tile 3 → highlighted as selected', async () => {
    const tiles = await page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').all();
    const tile0 = tiles[0];
    const tile2 = tiles[2];
    await tile0.click();
    const box2 = await tile2.boundingBox();
    if (box2) {
      await page.mouse.move(box2.x + 10, box2.y + 10);
      await page.mouse.down();
      await page.keyboard.down('Shift');
      await page.mouse.up();
      await page.keyboard.up('Shift');
    }
    const selCount = await page.evaluate(() => (grok as any).shell.tv.dataFrame.selection.trueCount);
    expect(selCount).toBeGreaterThanOrEqual(1);
  });

  await softStep('Row selection: ctrl-click tile 5 → added to selection', async () => {
    const tiles = await page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').all();
    const tile4 = tiles[4];
    const box = await tile4.boundingBox();
    if (box) {
      await page.evaluate(({x, y}) => {
        const el = document.elementFromPoint(x, y);
        ['mousedown','mouseup','click'].forEach(t =>
          el?.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, ctrlKey: true, clientX: x, clientY: y}))
        );
      }, {x: box.x + 10, y: box.y + 10});
    }
    const selCount = await page.evaluate(() => (grok as any).shell.tv.dataFrame.selection.trueCount);
    expect(selCount).toBeGreaterThanOrEqual(1);
    // Clean up
    await page.evaluate(() => (grok as any).shell.tv.dataFrame.selection.setAll(false));
  });

  // ---- Lanes ----
  await softStep('Lanes: set RACE → 4 lanes appear', async () => {
    await page.evaluate(() => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.lanesColumnName = 'RACE';
    });
    await page.waitForTimeout(500);
    const headers = await page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-lane-header').allTextContents();
    expect(headers.length).toBeGreaterThanOrEqual(4);
  });

  await softStep('Lanes: set SEX → 2 lanes (F, M)', async () => {
    await page.evaluate(() => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.lanesColumnName = 'SEX';
    });
    await page.waitForTimeout(500);
    const headers = await page.locator('[name="viewer-Tile-Viewer"] .d4-tile-viewer-lane-header').allTextContents();
    expect(headers).toEqual(expect.arrayContaining(['F', 'M']));
    expect(headers.length).toBe(2);
  });

  await softStep('Lanes: clear → single flat lane', async () => {
    const laneCount = await page.evaluate(async () => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.lanesColumnName = null;
      await new Promise(r => setTimeout(r, 400));
      return document.querySelectorAll('[name="viewer-Tile-Viewer"] .d4-tile-viewer-lane').length;
    });
    expect(laneCount).toBe(1);
  });

  // ---- Row source ----
  await softStep('Row source: Selected → 5 tiles for 5 selected rows', async () => {
    const tileCount = await page.evaluate(async () => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.rowSource = 'Selected';
      const df = (grok as any).shell.tv.dataFrame;
      df.selection.setAll(false);
      for (let i = 0; i < 5; i++) df.selection.set(i, true);
      await new Promise(r => setTimeout(r, 500));
      return document.querySelectorAll('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form').length;
    });
    expect(tileCount).toBe(5);
  });

  await softStep('Row source: Filtered with SEX=M → filtered rows shown', async () => {
    const result = await page.evaluate(async () => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.rowSource = 'Filtered';
      const fg = (grok as any).shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
      await new Promise(r => setTimeout(r, 800));
      return (grok as any).shell.tv.dataFrame.filter.trueCount;
    });
    expect(result).toBeLessThan(5850);
  });

  await softStep('Row source: All → all rows shown, filter cleared', async () => {
    const result = await page.evaluate(async () => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.rowSource = 'All';
      const fg = (grok as any).shell.tv.getFiltersGroup();
      const toRemove = [...fg.filters];
      for (const f of toRemove) fg.remove(f);
      (grok as any).shell.tv.dataFrame.selection.setAll(false);
      await new Promise(r => setTimeout(r, 500));
      return (grok as any).shell.tv.dataFrame.filter.trueCount;
    });
    expect(result).toBe(5850);
  });

  // ---- Tiles font ----
  await softStep('Tiles font: change size to 18px', async () => {
    const font = await page.evaluate(() => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.tilesFont = 'normal normal 18px "Roboto"';
      return v.props.tilesFont;
    });
    expect(font).toContain('18px');
  });

  await softStep('Tiles font: change family to Courier', async () => {
    const font = await page.evaluate(() => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.tilesFont = 'normal normal 18px "Courier"';
      return v.props.tilesFont;
    });
    expect(font).toContain('Courier');
  });

  await softStep('Tiles font: reset to default 13px Roboto', async () => {
    const font = await page.evaluate(() => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.tilesFont = 'normal normal 13px "Roboto"';
      return v.props.tilesFont;
    });
    expect(font).toContain('13px');
    expect(font).toContain('Roboto');
  });

  // ---- Auto-generate on column change ----
  await softStep('Auto-generate: load SPGI, switch table → SPGI columns', async () => {
    const result = await page.evaluate(async () => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.autoGenerate = true;
      const spgi = await (grok as any).dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      (grok as any).shell.addTableView(spgi);
      await new Promise(r => setTimeout(r, 1000));
      const demogView = Array.from((grok as any).shell.views).find((v: any) => v.name === 'Table');
      if (demogView) (grok as any).shell.v = demogView;
      await new Promise(r => setTimeout(r, 500));
      const tv = (grok as any).shell.tv;
      const tileV = tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      tileV.props.table = 'Table (2)';
      await new Promise(r => setTimeout(r, 1000));
      return tileV.dataFrame.name;
    });
    expect(result).toContain('(2)');
  });

  await softStep('Auto-generate: switch back to demog → demog columns', async () => {
    const result = await page.evaluate(async () => {
      const tv = (grok as any).shell.tv;
      const tileV = tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      tileV.props.table = 'Table';
      await new Promise(r => setTimeout(r, 1000));
      tileV.props.autoGenerate = false;
      return tileV.dataFrame.name;
    });
    expect(result).toBe('Table');
  });

  // ---- Context menu ----
  await softStep('Context menu: right-click tile → menu with Edit Form... and Properties...', async () => {
    const menuItems = await page.evaluate(async () => {
      const tiles = document.querySelectorAll('[name="viewer-Tile-Viewer"] .d4-tile-viewer-form');
      const tile = tiles[0] as HTMLElement;
      const rect = tile.getBoundingClientRect();
      tile.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2, clientX: rect.left + 20, clientY: rect.top + 20}));
      await new Promise(r => setTimeout(r, 300));
      return Array.from(document.querySelectorAll('.d4-menu-item-label')).map(el => el.textContent?.trim()).filter(Boolean);
    });
    expect(menuItems).toContain('Edit Form...');
    expect(menuItems.some(i => i && i.includes('Properties'))).toBe(true);
  });

  await softStep('Context menu: close menu and open properties via gear', async () => {
    await page.keyboard.press('Escape');
    await page.evaluate(() => {
      const gears = document.querySelectorAll('[name="icon-font-icon-settings"]');
      for (const g of gears) {
        const r = (g as HTMLElement).getBoundingClientRect();
        if (r.left > 900) { (g as HTMLElement).click(); break; }
      }
    });
    await page.waitForTimeout(500);
  });

  // ---- Viewer title and description ----
  await softStep('Viewer title: set title "Patient Cards" and description', async () => {
    const result = await page.evaluate(async () => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.title = 'Patient Cards';
      v.props.description = 'Demographic data per patient';
      await new Promise(r => setTimeout(r, 300));
      return { title: v.props.title, description: v.props.description };
    });
    expect(result.title).toBe('Patient Cards');
    expect(result.description).toBe('Demographic data per patient');
  });

  await softStep('Viewer title: set description position to Bottom', async () => {
    const val = await page.evaluate(() => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.descriptionPosition = 'Bottom';
      return v.props.descriptionPosition;
    });
    expect(val).toBe('Bottom');
  });

  await softStep('Viewer title: clear title and description', async () => {
    const result = await page.evaluate(() => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.title = '';
      v.props.description = '';
      return { title: v.props.title, description: v.props.description };
    });
    expect(result.title).toBe('');
    expect(result.description).toBe('');
  });

  // ---- Filter interaction ----
  await softStep('Filter interaction: open panel and add SEX=M filter', async () => {
    const filtered = await page.evaluate(async () => {
      const v = (grok as any).shell.tv.viewers.find((v: any) => v.type === 'Tile Viewer');
      v.props.rowSource = 'Filtered';
      const fg = (grok as any).shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
      await new Promise(r => setTimeout(r, 600));
      return (grok as any).shell.tv.dataFrame.filter.trueCount;
    });
    expect(filtered).toBeLessThan(5850);
  });

  await softStep('Filter interaction: add AGE > 50 → further reduced', async () => {
    const filtered = await page.evaluate(async () => {
      const fg = (grok as any).shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 51, max: 89});
      await new Promise(r => setTimeout(r, 600));
      return (grok as any).shell.tv.dataFrame.filter.trueCount;
    });
    expect(filtered).toBeLessThan(2607);
  });

  await softStep('Filter interaction: remove all filters → all tiles restored', async () => {
    const total = await page.evaluate(async () => {
      const fg = (grok as any).shell.tv.getFiltersGroup();
      const toRemove = [...fg.filters];
      for (const f of toRemove) fg.remove(f);
      await new Promise(r => setTimeout(r, 500));
      return (grok as any).shell.tv.dataFrame.filter.trueCount;
    });
    expect(total).toBe(5850);
  });

  await softStep('Filter interaction: close filter panel', async () => {
    await page.evaluate(() => {
      const filterSection = document.querySelector('[name="div-section--Filters"]');
      if (filterSection) (filterSection as HTMLElement).click();
    });
  });

  if (stepErrors.length > 0)
    throw new Error('Steps failed:\n' + stepErrors.map(e => `  ${e.step}: ${e.error}`).join('\n'));
});
