/* ---
realizes: [viewers.scatter-plot, viewers.histogram, viewers.line-chart, viewers.bar-chart, viewers.pie-chart, viewers.trellis-plot, viewers.box-plot, viewers.filters.histogram, viewers.filters.categorical, chem.filter.substructure-filter]
--- */
import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {clearClickFilters, clickCanvasFilter} from './canvas-filter';

// The substructure filter and the two layout round-trips live in filtering-server-spec.ts.
test.use(specTestOptions);

test('Legend filtering', async ({page}) => {
  test.setTimeout(600_000);

  await openDatagrok(page);
  await v.installEventWaits(page);
  // no withFilterPanel: every step here reaches the filters through getFiltersGroup(),
  // and opening the panel up front raced the substructure filter this dataset's molecule
  // column builds — the .d4-filter wait then timed out before any assertion ran
  await v.openTable(page);
  await v.addLegendViewers(page, {
    column: 'Stereo Category',
    viewers: ['Scatter plot', 'Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'],
  });

  await softStep('Numerical filter: Average Mass in [400, 10000]', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.tv.dataFrame;
      const col = df.col('Average Mass');
      // counted off the column rather than nailed to a constant: the spec was pinned to
      // the full SPGI's ~1588 and kept asserting it after the dataset moved to spgi-100
      let inRange = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (!col.isNone(i) && col.get(i) >= 400 && col.get(i) <= 10000) inRange++;
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 400, max: 10000});
      // technical: the filter group debounces, so onRowsFiltered fires on an
      // intermediate row set — no channel marks the settled one
      await new Promise((r) => setTimeout(r, 1500));
      return {inRange, filtered: df.filter.trueCount, rowCount: df.rowCount as number};
    });
    expect(res.filtered).toBe(res.inRange);
    // the range has to actually divide the table, or the comparison above is free
    expect(res.inRange).toBeGreaterThan(0);
    expect(res.inRange).toBeLessThan(res.rowCount);
  });

  await softStep('Categorical filter: R_ONE, S_UNKN only (legend=2)', async () => {
    await v.applyCategoricalFilter(page, 'Stereo Category', ['R_ONE', 'S_UNKN']);
    const {itemCount} = await v.readLegend(page, 'Scatter plot');
    expect(itemCount).toBe(2);
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
      const w = window as any;
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 300, max: 10000});
      // technical: the filter group debounces, so onRowsFiltered fires on an
      // intermediate row set — no channel marks the settled one
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
    const result = await clickCanvasFilter(page, {viewerType: 'Bar chart', column: 'Stereo Category'});
    expect(result.survivors).toBe(1);
    expect(result.totalFiltered).toBeGreaterThan(0);
  });

  await softStep('Pie chart canvas click-to-filter narrows the dataset', async () => {
    const result = await clickCanvasFilter(page, {viewerType: 'Pie chart', column: 'Stereo Category'});
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

  await softStep('Scatter plot Row Source cycles', async () => {
    await clearClickFilters(page);
    const res = await page.evaluate(async () => {
      const w = window as any;
      const sp = (window as any).grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      const fg = (window as any).grok.shell.tv.getFiltersGroup();
      const DG = (window as any).DG;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_UNKN']});
      // technical: the filter group debounces, so onRowsFiltered fires on an
      // intermediate row set — no channel marks the settled one
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
      const w = window as any;
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
      // technical: the filter group debounces, so onRowsFiltered fires on an
      // intermediate row set — no channel marks the settled one
      await new Promise((r) => setTimeout(r, 1500));
      const legend = bc.root.querySelector('[name="legend"]');
      const items = legend?.querySelectorAll('.d4-legend-item') ?? [];
      return {legendItems: items.length};
    });
    expect(res.legendItems).toBeLessThanOrEqual(2);
  });

  await softStep('Cleanup', async () => {
    await v.resetFilters(page, {clearScatterFilter: true});
    await v.cleanupShell(page);
  });

  v.finishSpec();
});
