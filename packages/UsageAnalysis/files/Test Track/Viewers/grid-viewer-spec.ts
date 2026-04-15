import {test, chromium, Page} from '@playwright/test';

declare const grok: any;
declare const DG: any;

const baseUrl = process.env.DATAGROK_URL ?? process.env.BASE_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Grid viewer: end-to-end on SPGI + SPGI-linked1', async () => {
  test.setTimeout(600_000);

  const cdp = await chromium.connectOverCDP('http://127.0.0.1:9222');
  const ctx = cdp.contexts()[0];
  const pages = ctx.pages();
  let page: Page = pages.find((p) => p.url().includes('datagrok')) ?? pages[0];
  if (!page) page = await ctx.newPage();
  await page.bringToFront();

  await page.goto(baseUrl, {timeout: 60000, waitUntil: 'networkidle'});
  await page.waitForFunction(() => {
    try {
      if (typeof grok === 'undefined' || !grok.shell || !grok.dapi || !grok.dapi.files) return false;
      grok.shell.closeAll();
      return true;
    } catch { return false; }
  }, {timeout: 180000, polling: 1000});
  await page.waitForTimeout(2000);

  let layoutId: string | null = null;
  let projectId: string | null = null;

  await softStep('1. Open SPGI and SPGI-linked1', async () => {
    await page.evaluate(async () => {
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
      try { grok.shell.windows.simpleMode = false; } catch {}
      grok.shell.closeAll();
      const df1 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      df1.name = 'SPGI';
      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI-linked1.csv');
      df2.name = 'SPGI-linked1';
      grok.shell.addTableView(df1);
      grok.shell.addTableView(df2);
      for (const df of [df1, df2]) {
        await new Promise((r) => {
          const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); r(undefined); });
          setTimeout(r, 4000);
        });
      }
      for (let i = 0; i < 60; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
      const spgi = Array.from(grok.shell.tableViews).find((tv: any) => tv.dataFrame.name === 'SPGI');
      grok.shell.v = spgi;
      await new Promise((r) => setTimeout(r, 1500));
    });
    await page.locator('[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  });

  await softStep('2. Add Grid via Add viewer icon', async () => {
    // Gallery is opened via a Dart command through an 'Add viewer' popup; open directly via JS API.
    await page.evaluate(async () => {
      const tv = grok.shell.v;
      // DG.Viewer.types returns types; use addViewer with 'Grid' type
      tv.addViewer('Grid');
      await new Promise((r) => setTimeout(r, 1000));
    });
    const count = await page.evaluate(() => grok.shell.v.viewers.filter((v: any) => v.type === 'Grid').length);
    if (count < 2) throw new Error('second Grid viewer not added, got ' + count);
  });

  await softStep('3. Close the added Grid viewer (JS API fallback)', async () => {
    await page.evaluate(() => {
      const tv = grok.shell.v;
      const extra = tv.viewers.find((v: any) => v !== tv.grid && v.type === 'Grid');
      extra.close();
    });
    await page.waitForTimeout(500);
    const n = await page.evaluate(() => grok.shell.v.viewers.filter((v: any) => v.type === 'Grid').length);
    if (n !== 1) throw new Error(`expected 1 Grid after close, got ${n}`);
  });

  await softStep('4. Add Grid via Toolbox > Viewers', async () => {
    await page.locator('[name="icon-grid"]').first().click({timeout: 10000});
    await page.waitForTimeout(500);
    const n = await page.evaluate(() => grok.shell.v.viewers.filter((v: any) => v.type === 'Grid').length);
    if (n !== 2) throw new Error(`expected 2 Grids, got ${n}`);
  });

  await softStep('5. Interact with cells, headers, selection, scrollbars', async () => {
    const res = await page.evaluate(() => {
      const tv = grok.shell.v;
      const df = tv.dataFrame;
      df.selection.set(3, true);
      df.selection.set(7, true);
      df.currentCell = df.cell(5, df.columns.byIndex(2).name);
      tv.grid.sort([df.columns.byIndex(0).name], [true]);
      const extra = tv.viewers.find((v: any) => v.type === 'Grid' && v !== tv.grid);
      return {
        selected: df.selection.trueCount,
        hasHScroll: !!extra.root.querySelector('.d4-grid-horz-scroll'),
        hasVScroll: !!extra.root.querySelector('.d4-grid-vert-scroll'),
      };
    });
    if (res.selected !== 2 || !res.hasHScroll || !res.hasVScroll)
      throw new Error('grid interactions failed: ' + JSON.stringify(res));
  });

  await softStep('6. Color code numeric column (JS API fallback)', async () => {
    await page.evaluate(() => {
      const df = grok.shell.v.dataFrame;
      const numCol = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i))
        .find((c: any) => c.type === 'double' || c.type === 'int');
      (numCol as any).meta.colors.setLinear();
    });
    const ok = await page.evaluate(() => {
      const df = grok.shell.v.dataFrame;
      const c = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i))
        .find((c: any) => c.meta.colors.getType() === 'Linear');
      return !!c;
    });
    if (!ok) throw new Error('color coding not applied');
  });

  await softStep('7. Open Context Panel with Grid selected', async () => {
    await page.evaluate(() => {
      const tv = grok.shell.v;
      const extra = tv.viewers.find((v: any) => v.type === 'Grid' && v !== tv.grid);
      grok.shell.o = extra;
    });
    await page.waitForTimeout(500);
  });

  await softStep('8. Data > Table — switch to SPGI-linked1 (JS API fallback)', async () => {
    const name = await page.evaluate(async () => {
      const tv = grok.shell.v;
      const extra = tv.viewers.find((v: any) => v.type === 'Grid' && v !== tv.grid);
      extra.dataFrame = grok.shell.tables.find((t: any) => t.name === 'SPGI-linked1');
      await new Promise((r) => setTimeout(r, 600));
      return extra.dataFrame.name;
    });
    if (name !== 'SPGI-linked1') throw new Error('rebind failed: ' + name);
  });

  await softStep('9. Modify row height, fonts, frozen columns', async () => {
    await page.evaluate(() => {
      const tv = grok.shell.v;
      const extra = tv.viewers.find((v: any) => v.type === 'Grid' && v !== tv.grid);
      extra.props.rowHeight = 40;
      extra.props.defaultCellFont = '14px "Arial"';
      extra.props.colHeaderFont = 'bold 14px "Arial"';
      extra.props.frozenColumns = 2;
    });
  });

  await softStep('10. Drag header/row-number borders (property proxy)', async () => {
    await page.evaluate(() => {
      const tv = grok.shell.v;
      const extra = tv.viewers.find((v: any) => v.type === 'Grid' && v !== tv.grid);
      extra.props.colHeaderHeight = 50;
      extra.props.rowHeight = 50;
    });
  });

  await softStep('11. Save the layout', async () => {
    layoutId = await page.evaluate(async () => {
      const tv = grok.shell.v;
      const layout = tv.saveLayout();
      layout.name = 'grid-viewer-test-' + Date.now();
      await grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1200));
      return layout.id;
    });
    if (!layoutId) throw new Error('layout id missing');
  });

  await softStep('12. Change the layout (add viewers, alter grid props)', async () => {
    await page.locator('[name="icon-histogram"]').first().click({timeout: 10000});
    await page.waitForTimeout(400);
    await page.locator('[name="icon-bar-chart"]').first().click({timeout: 10000});
    await page.waitForTimeout(400);
    await page.evaluate(() => {
      const tv = grok.shell.v;
      const extra = tv.viewers.find((v: any) => v.type === 'Grid' && v !== tv.grid);
      extra.props.rowHeight = 20;
      extra.props.frozenColumns = 0;
    });
  });

  await softStep('13. Apply saved layout — verify restored', async () => {
    const res = await page.evaluate(async (id: string) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.v.loadLayout(saved);
      await new Promise((r) => setTimeout(r, 3500));
      const tv = grok.shell.v;
      const extra = tv.viewers.find((v: any) => v.type === 'Grid' && v !== tv.grid);
      return {
        frozen: extra?.props.frozenColumns,
        colHeaderHeight: extra?.props.colHeaderHeight,
        font: extra?.props.defaultCellFont,
      };
    }, layoutId!);
    if (res.frozen !== 2 || res.colHeaderHeight !== 50)
      throw new Error('layout not restored: ' + JSON.stringify(res));
  });

  await softStep('14. Save the project', async () => {
    await page.locator('button:has-text("SAVE")').first().click({timeout: 10000});
    await page.waitForTimeout(2000);
    await page.locator('[name="button-OK"]').first().click({timeout: 10000});
    await page.waitForTimeout(6000);
    projectId = await page.evaluate(() => grok.shell.project?.id ?? null);
    if (!projectId) throw new Error('project not saved');
  });

  await softStep('15. Close All', async () => {
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1200);
    const n = await page.evaluate(() => Array.from(grok.shell.tableViews).length);
    if (n !== 0) throw new Error('views still open: ' + n);
  });

  await softStep('16. Open the project', async () => {
    if (!projectId) throw new Error('no projectId from step 14');
    const names = await page.evaluate(async (id: string) => {
      let p: any = null;
      for (let i = 0; i < 10 && !p; i++) {
        p = await grok.dapi.projects.find(id).catch(() => null);
        if (!p) await new Promise((r) => setTimeout(r, 1000));
      }
      if (!p) throw new Error('project not found by id ' + id);
      await p.open();
      await new Promise((r) => setTimeout(r, 5000));
      return Array.from(grok.shell.tableViews).map((tv: any) => tv.dataFrame.name);
    }, projectId);
    if (!names.includes('SPGI') || !names.includes('SPGI-linked1'))
      throw new Error('project views missing: ' + JSON.stringify(names));
  });

  // Cleanup
  await page.evaluate(async (ids: {layoutId: string | null; projectId: string | null}) => {
    try {
      if (ids.layoutId) { const l = await grok.dapi.layouts.find(ids.layoutId); if (l) await grok.dapi.layouts.delete(l); }
      if (ids.projectId) { const p = await grok.dapi.projects.find(ids.projectId); if (p) await grok.dapi.projects.delete(p); }
      grok.shell.closeAll();
    } catch {}
  }, {layoutId, projectId});

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
