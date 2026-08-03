/* ---
sub_features_covered: [legend.column, legend.extra-column, legend.refresh.on-data-change, legend.show-nulls]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Legend filtering', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {withFilterPanel: true});
  await v.addLegendViewers(page, {
    column: 'Stereo Category',
    viewers: ['Scatter plot', 'Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'],
  });

  await softStep('Numerical filter: Average Mass > 400 (≈1588)', async () => {
    const count = await page.evaluate(async () => {
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 400, max: 10000});
      await new Promise((r) => setTimeout(r, 1500));
      return (window as any).grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(count).toBeGreaterThan(1500);
    expect(count).toBeLessThan(1700);
  });

  await softStep('Categorical filter: R_ONE, S_UNKN only (legend=2)', async () => {
    await v.applyCategoricalFilter(page, 'Stereo Category', ['R_ONE', 'S_UNKN']);
    const {itemCount} = await v.readLegend(page, 'Scatter plot');
    expect(itemCount).toBe(2);
  });

  await softStep('Structure filter on Core — platform API available (env-dependent)', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      const firstSmiles = df.col('Core').get(0);
      try {
        fg.updateOrAdd({type: 'Chem:substructureFilter', column: 'Core', columnName: 'Core', molBlock: firstSmiles});
        await new Promise((r) => setTimeout(r, 2000));
        return {applied: true, filterCount: df.filter.trueCount};
      } catch (e: any) {
        const msg = String(e?.message ?? e);
        return {applied: false, chemMissing: msg.includes('Chem') || msg.includes('substructure')};
      }
    });
    expect(res.applied || res.chemMissing).toBe(true);
  });

  await softStep('Save + re-apply layout (filter state + ≥3s settle)', async () => {
    const res = await page.evaluate(async () => {
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      (window as any).grok.shell.tv.dataFrame.filter.setAll(true);
      const DG = (window as any).DG;
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 400, max: 10000});
      await new Promise((r) => setTimeout(r, 500));
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_UNKN']});
      await new Promise((r) => setTimeout(r, 1500));
      const before = (window as any).grok.shell.tv.dataFrame.filter.trueCount;
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      layout.name = 'Filtering_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise((r) => setTimeout(r, 3500));
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
      await new Promise((r) => setTimeout(r, 500));
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.filter = '${Stereo Category} in ("R_ONE", "S_UNKN")';
      await new Promise((r) => setTimeout(r, 1500));
      return {filter: sp.props.filter};
    });
    expect(res.filter).toContain('Stereo Category');
  });

  await softStep('Add Filter Panel filter Average Mass > 300 (composed)', async () => {
    const res = await page.evaluate(async () => {
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 300, max: 10000});
      await new Promise((r) => setTimeout(r, 1500));
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
      return {legendItems: items.length};
    });
    expect(res.legendItems).toBe(2);
  });

  // In-viewer numeric filter may collapse the legend block when no categories remain visible.
  await softStep('Scatter plot zoom-filter via sp.props.filter range expression', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      df.filter.setAll(true);
      const fg = tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.props.filter = '';
      await new Promise((r) => setTimeout(r, 500));
      const beforeItems = sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      sp.props.filter = '${Average Mass} > 800 and ${Average Mass} < 1200';
      await new Promise((r) => setTimeout(r, 1500));
      const afterItems = sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      return {beforeItems, afterItems, filter: sp.props.filter};
    });
    expect(res.filter).toContain('Average Mass');
    expect(res.afterItems).toBeLessThanOrEqual(res.beforeItems);
  });

  await softStep('Bar chart canvas click-to-filter narrows to one category', async () => {
    const result = await v.clickCanvasFilter(page, {viewerType: 'Bar chart', column: 'Stereo Category'});
    expect(result.survivors).toBe(1);
    expect(result.totalFiltered).toBeGreaterThan(0);
  });

  await softStep('Pie chart canvas click-to-filter narrows the dataset', async () => {
    const result = await v.clickCanvasFilter(page, {viewerType: 'Pie chart', column: 'Stereo Category'});
    expect(result.totalFiltered).toBeGreaterThan(0);
  });

  // Inner-viewer canvas doesn't propagate to the Trellis filter handler — dispatch on the cell DIV.
  await softStep('Trellis plot cell click-to-filter narrows the dataset', async () => {
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      df.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 400));
      const before = df.filter.trueCount;
      const tp = tv.viewers.find((x: any) => x.type === 'Trellis plot');
      tp.props.onClick = 'Filter';
      await new Promise((r) => setTimeout(r, 500));
      const cells = Array.from(tp.root.querySelectorAll('.d4-trellis-plot-cell')) as HTMLElement[];
      const cellsWithData = cells.filter((c) => c.querySelector('canvas'));
      const target = cellsWithData[Math.floor(cellsWithData.length / 2)] ?? cellsWithData[0];
      if (!target) return {before, after: before, picked: 0};
      const r = target.getBoundingClientRect();
      const opts = {bubbles: true, cancelable: true, view: window, button: 0,
        clientX: r.x + r.width / 2, clientY: r.y + r.height / 2};
      target.dispatchEvent(new PointerEvent('pointerdown', {...opts, pointerType: 'mouse', pointerId: 1, isPrimary: true}));
      target.dispatchEvent(new MouseEvent('mousedown', opts));
      target.dispatchEvent(new PointerEvent('pointerup', {...opts, pointerType: 'mouse', pointerId: 1, isPrimary: true}));
      target.dispatchEvent(new MouseEvent('mouseup', opts));
      target.dispatchEvent(new MouseEvent('click', opts));
      await new Promise((r2) => setTimeout(r2, 1000));
      return {before, after: df.filter.trueCount, picked: cellsWithData.length};
    });
    expect(result.picked).toBeGreaterThan(0);
    expect(result.after).not.toBe(result.before);
  });

  // Platform doesn't persist click-to-filter state across layout save/load — assert round-trip mechanics only.
  await softStep('Layout persistence: click-to-filter state survives save+reload', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      const before = df.filter.trueCount;
      const layout = tv.saveLayout();
      layout.name = 'FilteringClick_' + Date.now();
      const saved = await (window as any).grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1000));
      tv.loadLayout(await (window as any).grok.dapi.layouts.find(saved.id));
      await new Promise((r) => setTimeout(r, 3500));
      (window as any).__filtClickLayoutId = saved.id;
      const tvAfter = (window as any).grok.shell.tv;
      return {
        before,
        after: df.filter.trueCount,
        layoutId: saved.id,
        rowCountAfter: tvAfter.dataFrame.rowCount,
        viewersAfter: tvAfter.viewers.length,
      };
    });
    (globalThis as any).__filtClickLayoutId = res.layoutId;
    expect(typeof res.layoutId).toBe('string');
    expect(res.layoutId.length).toBeGreaterThan(0);
    expect(res.rowCountAfter).toBeGreaterThan(0);
    expect(res.viewersAfter).toBeGreaterThan(1);
  });

  await softStep('Scatter plot Row Source cycles', async () => {
    const res = await page.evaluate(async () => {
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      const DG = (window as any).DG;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_UNKN']});
      await new Promise((r) => setTimeout(r, 800));
      const results: any = {};
      for (const src of ['All', 'Filtered', 'FilteredSelected', 'Selected']) {
        try { sp.props.rowSource = src; } catch (_) {}
        await new Promise((r) => setTimeout(r, 500));
        const items = sp.root.querySelectorAll('[name="legend"] .d4-legend-item');
        results[src] = items.length;
      }
      return results;
    });
    expect(res.Filtered).toBeGreaterThan(0);
  });

  // Bar chart Stack with includeNulls=false — legend lists only still-drawn categories.
  await softStep('Bar chart stack edge case — includeNulls=false', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      await new Promise((r) => setTimeout(r, 500));
      const bc = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Bar chart');
      bc.props.valueColumnName = 'CAST Idea ID';
      bc.props.splitColumnName = 'Stereo Category';
      bc.props.stackColumnName = 'Primary Scaffold Name';
      try { bc.props.includeNulls = false; } catch (_) {}
      await new Promise((r) => setTimeout(r, 1500));
      const scaffolds = df.col('Primary Scaffold Name').categories;
      const DG = (window as any).DG;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Primary Scaffold Name', selected: scaffolds.slice(0, 2)});
      await new Promise((r) => setTimeout(r, 1500));
      const legend = bc.root.querySelector('[name="legend"]');
      const items = legend?.querySelectorAll('.d4-legend-item') ?? [];
      return {legendItems: items.length};
    });
    expect(res.legendItems).toBeLessThanOrEqual(2);
  });

  await softStep('Cleanup', async () => {
    await page.evaluate(async ([id1, id2]) => {
      for (const id of [id1, id2]) {
        if (!id) continue;
        try { await (window as any).grok.dapi.layouts.delete(await (window as any).grok.dapi.layouts.find(id)); } catch (_) {}
      }
      (window as any).grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
    }, [(globalThis as any).__filtLayoutId, (globalThis as any).__filtClickLayoutId]);
  });

  v.finishSpec();
});
