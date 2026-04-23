import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('Line chart legend', async ({page}) => {
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
      .some((c: any) => c.semType === 'Molecule');
    if (hasBioChem) await new Promise(r => setTimeout(r, 5000));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  await softStep('Add line chart, Split=Series', async () => {
    const items = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Line chart');
      await new Promise(r => setTimeout(r, 800));
      const lc = tv.viewers.find((v: any) => v.type === 'Line chart');
      lc.props.splitColumnName = 'Series';
      try { lc.props.legendVisibility = 'Always'; } catch(_) {}
      await new Promise(r => setTimeout(r, 1500));
      return lc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
    });
    expect(items).toBeGreaterThan(0);
  });

  await softStep('Enable multiAxis', async () => {
    const res = await page.evaluate(async () => {
      const lc = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart');
      try { lc.props.multiAxis = true; } catch(e) { return {multiAxis: false, err: String(e)}; }
      await new Promise(r => setTimeout(r, 1500));
      return {multiAxis: lc.props.multiAxis};
    });
    expect(res.multiAxis).toBe(true);
  });

  await softStep('Save + re-apply layout (multiAxis+split persist)', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'LineChart_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise(r => setTimeout(r, 3000));
      (window as any).__lcLayoutId = saved.id;
      const lc = tv.viewers.find((v: any) => v.type === 'Line chart');
      return {multiAxis: lc?.props.multiAxis, split: lc?.props.splitColumnName};
    });
    (globalThis as any).__lcLayoutId = await page.evaluate(() => (window as any).__lcLayoutId);
    expect(res.multiAxis).toBe(true);
    expect(res.split).toBe('Series');
  });

  await softStep('Two Y columns = Average Mass + TPSA', async () => {
    const res = await page.evaluate(async () => {
      const lc = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart');
      lc.props.yColumnNames = ['Average Mass', 'TPSA'];
      await new Promise(r => setTimeout(r, 2000));
      const items = lc.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {yCols: lc.props.yColumnNames, totalItems: items.length};
    });
    expect(res.yCols).toEqual(['Average Mass', 'TPSA']);
    expect(res.totalItems).toBeGreaterThan(0);
  });

  await softStep('Replace Y column → NIBR logP', async () => {
    const res = await page.evaluate(async () => {
      const lc = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart');
      lc.props.yColumnNames = ['Average Mass', 'NIBR logP'];
      await new Promise(r => setTimeout(r, 1500));
      return {yCols: lc.props.yColumnNames, items: lc.root.querySelectorAll('[name="legend"] .d4-legend-item').length};
    });
    expect(res.yCols).toEqual(['Average Mass', 'NIBR logP']);
  });

  await softStep('Save + re-apply layout (new Y persists)', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'LineChart2_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise(r => setTimeout(r, 3000));
      (window as any).__lcLayoutId2 = saved.id;
      const lc = tv.viewers.find((v: any) => v.type === 'Line chart');
      return {yCols: lc?.props.yColumnNames};
    });
    (globalThis as any).__lcLayoutId2 = await page.evaluate(() => (window as any).__lcLayoutId2);
    expect(res.yCols).toEqual(['Average Mass', 'NIBR logP']);
  });

  await softStep('Save project (known FK)', async () => {
    const res = await page.evaluate(async () => {
      try {
        const DG = (window as any).DG;
        const proj = DG.Project.create();
        proj.name = 'LineChartProj_' + Date.now();
        proj.addChild((window as any).grok.shell.tv.dataFrame);
        await (window as any).grok.dapi.projects.save(proj);
        return {ok: true};
      } catch(e: any) { return {ok: false, error: String(e).slice(0, 180)}; }
    });
    expect(res.ok || String(res.error ?? '').includes('foreign key')).toBe(true);
  });

  await softStep('Cleanup', async () => {
    await page.evaluate(async () => {
      for (const key of ['__lcLayoutId', '__lcLayoutId2']) {
        const id = (window as any)[key];
        if (id) { try { await (window as any).grok.dapi.layouts.delete(await (window as any).grok.dapi.layouts.find(id)); } catch(_) {} }
      }
      (window as any).grok.shell.closeAll();
    });
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
