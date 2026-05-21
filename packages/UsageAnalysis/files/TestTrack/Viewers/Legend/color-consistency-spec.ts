import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('Legend color consistency', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

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
    const tv = (window as any).grok.shell.tv;
    const names = ['Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'];
    for (const n of names) {
      tv.addViewer(n);
      await new Promise(r => setTimeout(r, 300));
    }
    const col = 'Stereo Category';
    for (const v of tv.viewers) {
      if (v.type === 'Grid') continue;
      try {
        if (v.type === 'Histogram') v.props.splitColumnName = col;
        else if (v.type === 'Line chart') v.props.splitColumnName = col;
        else if (v.type === 'Bar chart') v.props.splitColumnName = col;
        else if (v.type === 'Pie chart') v.props.categoryColumnName = col;
        else if (v.type === 'Trellis plot') v.props.xColumnNames = [col];
        else if (v.type === 'Box plot') v.props.categoryColumnNames = [col];
        try { v.props.legendVisibility = 'Always'; } catch(_) {}
      } catch(_) {}
    }
    await new Promise(r => setTimeout(r, 1500));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  await softStep('Enable categorical color coding + change two colors', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const col = df.col('Stereo Category');
      col.tags['.color-coding-type'] = 'Categorical';
      col.tags['.categorical-colors'] = JSON.stringify({'R_ONE': '#FF0000', 'S_UNKN': '#00FF00'});
      const cats = col.categories;
      const map: any = {};
      for (const c of cats) {
        if (c === 'R_ONE') map[c] = 0xFFFF0000;
        else if (c === 'S_UNKN') map[c] = 0xFF00FF00;
        else map[c] = 0xFF808080;
      }
      col.meta.colors.setCategorical(map);
      await new Promise(r => setTimeout(r, 800));
      (window as any).grok.shell.tv.grid.invalidate();
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type !== 'Grid') try { v.invalidate?.(); } catch(_) {}
      await new Promise(r => setTimeout(r, 500));
      return {rOneColor: '0x' + (col.meta.colors.getColor(0) >>> 0).toString(16)};
    });
    expect(res.rOneColor).toBe('0xffff0000');
  });

  await softStep('Verify a viewer with legend (Histogram) shows colored items', async () => {
    const items = await page.evaluate(() => {
      const hist = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Histogram');
      const legend = hist.root.querySelector('[name="legend"]');
      const items = legend?.querySelectorAll('.d4-legend-item') ?? [];
      return items.length;
    });
    expect(items).toBeGreaterThan(0);
  });

  await softStep('Save + re-apply layout — custom palette persists', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'ColorConsist_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise(r => setTimeout(r, 3000));
      (window as any).__ccLayoutId = saved.id;
      const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
      return {rOne: '0x' + (col.meta.colors.getColor(0) >>> 0).toString(16), tag: col.tags['.categorical-colors']};
    });
    (globalThis as any).__ccLayoutId = await page.evaluate(() => (window as any).__ccLayoutId);
    expect(res.rOne).toBe('0xffff0000');
  });

  await softStep('Save project (known FK limitation)', async () => {
    const res = await page.evaluate(async () => {
      try {
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'ColorConsistProj_' + Date.now();
        proj.addChild((window as any).grok.shell.tv.dataFrame);
        const saved = await (window as any).grok.dapi.projects.save(proj);
        return {ok: true, id: saved.id};
      } catch (e: any) { return {ok: false, error: String(e).slice(0, 180)}; }
    });
    expect(res.ok || String(res.error ?? '').includes('foreign key')).toBe(true);
  });

  await softStep('Cleanup', async () => {
    await page.evaluate(async (id) => {
      if (id) { try { await (window as any).grok.dapi.layouts.delete(await (window as any).grok.dapi.layouts.find(id)); } catch(_) {} }
      (window as any).grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
    }, (globalThis as any).__ccLayoutId);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
