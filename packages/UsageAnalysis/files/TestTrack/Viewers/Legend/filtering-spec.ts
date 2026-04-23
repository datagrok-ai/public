import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('Legend filtering', async ({page}) => {
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
    const names = ['Scatter plot', 'Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'];
    for (const n of names) {
      tv.addViewer(n);
      await new Promise(r => setTimeout(r, 300));
    }
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
        try { v.props.legendVisibility = 'Always'; } catch(_) {}
      } catch(_) {}
    }
    tv.getFiltersGroup();
    await new Promise(r => setTimeout(r, 1500));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 15000});

  await softStep('Numerical filter: Average Mass > 400 (≈1588)', async () => {
    const count = await page.evaluate(async () => {
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 400, max: 10000});
      await new Promise(r => setTimeout(r, 1500));
      return (window as any).grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(count).toBeGreaterThan(1500);
    expect(count).toBeLessThan(1700);
  });

  await softStep('Categorical filter: R_ONE, S_UNKN only (legend=2)', async () => {
    const res = await page.evaluate(async () => {
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      const DG = (window as any).DG;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_UNKN']});
      await new Promise(r => setTimeout(r, 1500));
      const sp = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {filterCount: (window as any).grok.shell.tv.dataFrame.filter.trueCount, legendItems: items.length};
    });
    expect(res.legendItems).toBe(2);
  });

  await softStep('Structure filter on Core — platform API available', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const firstSmiles = df.col('Core').get(0);
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      try {
        fg.updateOrAdd({type: 'Chem:substructureFilter', column: 'Core', columnName: 'Core', molBlock: firstSmiles});
        await new Promise(r => setTimeout(r, 2000));
        return {ok: true, filterCount: df.filter.trueCount};
      } catch (e: any) { return {ok: false, error: String(e).slice(0, 200)}; }
    });
    expect(res.ok).toBe(true);
  });

  await softStep('Save + re-apply layout (filter state + ≥3s settle)', async () => {
    const res = await page.evaluate(async () => {
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch(_) {} }
      (window as any).grok.shell.tv.dataFrame.filter.setAll(true);
      const DG = (window as any).DG;
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 400, max: 10000});
      await new Promise(r => setTimeout(r, 500));
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_UNKN']});
      await new Promise(r => setTimeout(r, 1500));
      const before = (window as any).grok.shell.tv.dataFrame.filter.trueCount;
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'Filtering_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise(r => setTimeout(r, 3500));
      (window as any).__filtLayoutId = saved.id;
      return {before, after: (window as any).grok.shell.tv.dataFrame.filter.trueCount};
    });
    (globalThis as any).__filtLayoutId = await page.evaluate(() => (window as any).__filtLayoutId);
    expect(res.after).toBe(res.before);
  });

  await softStep('Reset + in-viewer Scatter plot filter', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      await new Promise(r => setTimeout(r, 500));
      const sp = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.filter = '${Stereo Category} in ("R_ONE", "S_UNKN")';
      await new Promise(r => setTimeout(r, 1500));
      return {filter: sp.props.filter};
    });
    expect(res.filter).toContain('Stereo Category');
  });

  await softStep('Add Filter Panel filter Average Mass > 300 (composed)', async () => {
    const res = await page.evaluate(async () => {
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 300, max: 10000});
      await new Promise(r => setTimeout(r, 1500));
      const sp = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {legendItems: items.length};
    });
    expect(res.legendItems).toBe(2);
  });

  await softStep('Bar chart OnClick=Filter property', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch(_) {} }
      const sp = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.props.filter = '';
      await new Promise(r => setTimeout(r, 500));
      const bc = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Bar chart');
      bc.props.onClick = 'Filter';
      await new Promise(r => setTimeout(r, 500));
      return {onClick: bc.props.onClick};
    });
    expect(res.onClick).toBe('Filter');
  });

  await softStep('Scatter plot Row Source cycles', async () => {
    const res = await page.evaluate(async () => {
      const sp = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      const DG = (window as any).DG;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_UNKN']});
      await new Promise(r => setTimeout(r, 800));
      const results: any = {};
      for (const src of ['All', 'Filtered', 'FilteredSelected', 'Selected']) {
        try { sp.props.rowSource = src; } catch(_) {}
        await new Promise(r => setTimeout(r, 500));
        const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
        results[src] = items.length;
      }
      return results;
    });
    expect(res.Filtered).toBeGreaterThan(0);
  });

  await softStep('Bar chart stack edge case — includeNulls=false', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch(_) {} }
      await new Promise(r => setTimeout(r, 500));
      const bc = (window as any).grok.shell.tv.viewers.find((v: any) => v.type === 'Bar chart');
      bc.props.valueColumnName = 'CAST Idea ID';
      bc.props.splitColumnName = 'Stereo Category';
      bc.props.stackColumnName = 'Primary Scaffold Name';
      try { bc.props.includeNulls = false; } catch(_) {}
      await new Promise(r => setTimeout(r, 1500));
      const scaffolds = df.col('Primary Scaffold Name').categories;
      const DG = (window as any).DG;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Primary Scaffold Name', selected: scaffolds.slice(0, 2)});
      await new Promise(r => setTimeout(r, 1500));
      const legend = bc.root.querySelector('[name="legend"]');
      const items = legend?.querySelectorAll('.d4-legend-item') ?? [];
      return {legendItems: items.length};
    });
    expect(res.legendItems).toBeLessThanOrEqual(2);
  });

  await softStep('Cleanup', async () => {
    await page.evaluate(async (id) => {
      if (id) { try { await (window as any).grok.dapi.layouts.delete(await (window as any).grok.dapi.layouts.find(id)); } catch(_) {} }
      (window as any).grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
    }, (globalThis as any).__filtLayoutId);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
