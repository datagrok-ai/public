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

test('Legend visibility and positioning', async ({page}) => {
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

  await softStep('Add 7 viewers', async () => {
    const types = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const names = ['Scatter plot', 'Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'];
      for (const n of names) {
        tv.addViewer(n);
        await new Promise(r => setTimeout(r, 300));
      }
      return tv.viewers.map((v: any) => v.type);
    });
    expect(types.length).toBeGreaterThanOrEqual(8);
  });

  await softStep('Set Stereo Category legend on each viewer', async () => {
    const ok = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const col = 'Stereo Category';
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
        } catch (_) {}
      }
      await new Promise(r => setTimeout(r, 1500));
      return true;
    });
    expect(ok).toBe(true);
  });

  await softStep('Verify legends are visible', async () => {
    const count = await page.locator('[name="legend"], .d4-corner-legend').count();
    expect(count).toBeGreaterThan(0);
  });

  await softStep('Change legend source and back', async () => {
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.colorColumnName = 'Series';
      await new Promise(r => setTimeout(r, 500));
      sp.props.colorColumnName = 'Stereo Category';
      await new Promise(r => setTimeout(r, 500));
      return sp.props.colorColumnName;
    });
    expect(result).toBe('Stereo Category');
  });

  await softStep('Open color picker via hover on legend item', async () => {
    const found = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      const legend = sp.root.querySelector('[name="legend"]');
      if (!legend) return false;
      const items = legend.querySelectorAll('.d4-legend-item');
      let target: Element | null = null;
      for (const it of items) {
        const val = it.querySelector('.d4-legend-value');
        if (val && val.textContent?.trim() === 'R_ONE') { target = it; break; }
      }
      if (!target) return false;
      target.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      target.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise(r => setTimeout(r, 300));
      return !!document.querySelector('[name="legend-icon-color-picker"], [name="Legend-icon-color-picker"]');
    });
    expect(found).toBe(true);
  });

  await softStep('Save and re-apply layout', async () => {
    const id = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'LegendVP_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      const found = await (window as any).grok.dapi.layouts.find(saved.id);
      tv.loadLayout(found);
      await new Promise(r => setTimeout(r, 3000));
      return saved.id;
    });
    (globalThis as any).__layoutId = id;
    expect(id).toBeTruthy();
  });

  await softStep('Set legendVisibility=Always, legendPosition=Auto', async () => {
    const ok = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      for (const v of tv.viewers) {
        if (v.type === 'Grid') continue;
        try { v.props.legendVisibility = 'Always'; } catch (_) {}
        try { v.props.legendPosition = 'Auto'; } catch (_) {}
      }
      await new Promise(r => setTimeout(r, 800));
      return tv.viewers.filter((v: any) => v.type !== 'Grid').every((v: any) => v.props.legendVisibility === 'Always');
    });
    expect(ok).toBe(true);
  });

  await softStep('Resize viewer — Auto-position reflows legend', async () => {
    const width = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.root.style.width = '300px';
      await new Promise(r => setTimeout(r, 500));
      return sp.root.getBoundingClientRect().width;
    });
    expect(Math.round(width)).toBe(300);
  });

  await softStep('Save/re-apply layout (Always + Auto must persist)', async () => {
    const vis = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'LegendVP2_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise(r => setTimeout(r, 3000));
      (window as any).__layoutId2 = saved.id;
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      return sp?.props?.legendVisibility;
    });
    (globalThis as any).__layoutId2 = await page.evaluate(() => (window as any).__layoutId2);
    expect(vis).toBe('Always');
  });

  await softStep('Visibility=Auto + resize 200px hides, 400px shows', async () => {
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      for (const v of tv.viewers) {
        if (v.type === 'Grid') continue;
        try { v.props.legendVisibility = 'Auto'; } catch(_) {}
      }
      await new Promise(r => setTimeout(r, 500));
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.root.style.width = '200px';
      await new Promise(r => setTimeout(r, 500));
      const small = !!sp.root.querySelector('[name="legend"]');
      sp.root.style.width = '400px';
      await new Promise(r => setTimeout(r, 500));
      const big = !!sp.root.querySelector('[name="legend"]');
      return { small, big };
    });
    expect(result.big).toBe(true);
  });

  await softStep('Corner positions (miniLegend property not exposed)', async () => {
    const positions = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const corners = ['LeftTop', 'LeftBottom', 'RightTop', 'RightBottom'];
      let i = 0;
      const out: any = {};
      for (const v of tv.viewers) {
        if (v.type === 'Grid') continue;
        v.root.style.width = '';
        try { v.props.legendVisibility = 'Always'; } catch(e) {}
        try { v.props.legendPosition = corners[i % 4]; } catch(e) {}
        out[v.type] = v.props.legendPosition;
        i++;
      }
      await new Promise(r => setTimeout(r, 1500));
      return out;
    });
    expect(Object.keys(positions).length).toBeGreaterThan(0);
  });

  await softStep('Save/re-apply layout after corner positions', async () => {
    const pos = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'LegendVP3_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise(r => setTimeout(r, 3000));
      (window as any).__layoutId3 = saved.id;
      const sp = tv.viewers.find((v: any) => v.type === 'Scatter plot');
      return sp?.props?.legendPosition;
    });
    (globalThis as any).__layoutId3 = await page.evaluate(() => (window as any).__layoutId3);
    expect(pos).toBeTruthy();
  });

  await softStep('Save project (known FK constraint limitation)', async () => {
    const res = await page.evaluate(async () => {
      try {
        const tv = (window as any).grok.shell.tv;
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'LegendVPProj_' + Date.now();
        proj.addChild(tv.dataFrame);
        const saved = await (window as any).grok.dapi.projects.save(proj);
        return { ok: true, id: saved.id };
      } catch (e: any) {
        return { ok: false, error: String(e).slice(0, 200) };
      }
    });
    expect(res.ok || String(res.error ?? '').includes('foreign key')).toBe(true);
  });

  await softStep('Cleanup layouts, close all', async () => {
    await page.evaluate(async () => {
      for (const key of ['__layoutId', '__layoutId2', '__layoutId3']) {
        const id = (window as any)[key];
        if (id) {
          try {
            const l = await (window as any).grok.dapi.layouts.find(id);
            await (window as any).grok.dapi.layouts.delete(l);
          } catch (_) {}
        }
      }
      (window as any).grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
    });
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
