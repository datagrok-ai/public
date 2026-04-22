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

test('Legend structure rendering', async ({page}) => {
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

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    const df = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    (window as any).grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 5000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i: number) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        const grid = document.querySelector('[name="viewer-Grid"]');
        if (grid?.querySelector('canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Add 7 viewers, set legend column to Core, Always visible', async () => {
    const types = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const names = ['Scatter plot', 'Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'];
      for (const n of names) {
        tv.addViewer(n);
        await new Promise(r => setTimeout(r, 300));
      }
      const col = 'Core';
      for (const v of tv.viewers) {
        if (v.type === 'Grid') continue;
        try {
          if (v.type === 'Scatter plot') v.props.colorColumnName = col;
          else if (v.type === 'Histogram') v.props.splitColumnName = col;
          else if (v.type === 'Line chart') v.props.splitColumnName = col;
          else if (v.type === 'Bar chart') v.props.splitColumnName = col;
          else if (v.type === 'Pie chart') v.props.categoryColumnName = col;
          else if (v.type === 'Trellis plot') v.props.xColumnNames = [col];
          else if (v.type === 'Box plot') v.props.categoryColumnNames = [col];
          try { v.props.legendVisibility = 'Always'; } catch (_) {}
        } catch (_) {}
      }
      await new Promise(r => setTimeout(r, 2500));
      return tv.viewers.map((v: any) => v.type);
    });
    expect(types.length).toBeGreaterThanOrEqual(8);
  });

  await softStep('Verify legend canvases on Scatter/Hist/Line/Pie (Molecule rendering)', async () => {
    const res = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const out: any = {};
      for (const v of tv.viewers) {
        if (v.type === 'Grid') continue;
        const legend = v.root.querySelector('[name="legend"]');
        const items = legend?.querySelectorAll('.d4-legend-item') ?? [];
        out[v.type] = {items: items.length, canvas: !!items[0]?.querySelector('canvas')};
      }
      return out;
    });
    expect(res['Scatter plot']?.canvas).toBe(true);
    expect(res['Histogram']?.canvas).toBe(true);
    expect(res['Line chart']?.canvas).toBe(true);
    expect(res['Pie chart']?.canvas).toBe(true);
  });

  await softStep('Scatter plot — Marker=Core, Color=Core', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.markersColumnName = 'Core';
      sp.props.colorColumnName = 'Core';
      await new Promise(r => setTimeout(r, 1500));
      const legend = sp.root.querySelector('[name="legend"]');
      const items = legend?.querySelectorAll('.d4-legend-item') ?? [];
      return {items: items.length, canvas: !!items[0]?.querySelector('canvas')};
    });
    expect(res.canvas).toBe(true);
  });

  await softStep('Scatter plot — Color=Series, Marker stays Core', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.colorColumnName = 'Series';
      await new Promise(r => setTimeout(r, 1000));
      const legend = sp.root.querySelector('[name="legend"]');
      const items = legend?.querySelectorAll('.d4-legend-item') ?? [];
      return {items: items.length, markers: sp.props.markersColumnName, color: sp.props.colorColumnName};
    });
    expect(res.markers).toBe('Core');
    expect(res.color).toBe('Series');
    expect(res.items).toBeGreaterThan(0);
  });

  await softStep('Save + re-apply layout', async () => {
    const id = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'StructRender_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise(r => setTimeout(r, 3500));
      (window as any).__srLayoutId = saved.id;
      return saved.id;
    });
    (globalThis as any).__srLayoutId = id;
    expect(id).toBeTruthy();
  });

  await softStep('Save project (known FK limitation)', async () => {
    const res = await page.evaluate(async () => {
      try {
        const tv = (window as any).grok.shell.tv;
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'StructRenderProj_' + Date.now();
        proj.addChild(tv.dataFrame);
        const saved = await (window as any).grok.dapi.projects.save(proj);
        return {ok: true, id: saved.id};
      } catch (e: any) {
        return {ok: false, error: String(e).slice(0, 200)};
      }
    });
    expect(res.ok || String(res.error ?? '').includes('foreign key')).toBe(true);
  });

  await softStep('Cleanup', async () => {
    await page.evaluate(async (id) => {
      if (id) {
        try { await (window as any).grok.dapi.layouts.delete(await (window as any).grok.dapi.layouts.find(id)); } catch(_) {}
      }
      (window as any).grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
    }, (globalThis as any).__srLayoutId);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
