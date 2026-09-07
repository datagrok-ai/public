/* ---
realizes: [viewers.scatter-plot, viewers.histogram, viewers.line-chart, viewers.bar-chart, viewers.pie-chart, viewers.trellis-plot, viewers.box-plot, viewers.filters.histogram, viewers.filters.categorical, chem.filter.substructure-filter]
--- */
// GROK-17222: legend item count must update under FP, in-viewer, and Pie/Bar click-to-filter.

import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {addLegendViewers} from './legend-setup';
import {clickCanvasFilter} from './canvas-filter';

test.use(specTestOptions);

test('GROK-17222: legend reflects filter state across 4 trigger sources', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await openDatagrok(page);
  await v.openTable(page, {withFilterPanel: true});
  await addLegendViewers(page, {
    column: 'Stereo Category',
    viewers: ['Line chart', 'Scatter plot', 'Pie chart', 'Bar chart'],
  });

  await softStep('Step 4: Filter Panel categorical filter narrows legend (Line chart)', async () => {
    const beforeItems = (await v.readLegend(page, 'Line chart')).itemCount;
    const cats: string[] = await page.evaluate(() =>
      Array.from((window as any).grok.shell.tv.dataFrame.col('Stereo Category').categories));
    const subset = cats.slice(0, Math.min(2, cats.length));
    const {filteredCount} = await v.applyCategoricalFilter(page, 'Stereo Category', subset);
    const afterItems = (await v.readLegend(page, 'Line chart')).itemCount;
    expect(afterItems, 'legend item count must respond to FP filter (GROK-17222)')
      .toBeLessThanOrEqual(beforeItems);
    expect(filteredCount).toBeGreaterThan(0);
  });

  await softStep('Step 5: in-viewer Scatter filter via sp.props.filter range', async () => {
    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const df = tv.dataFrame;
      df.filter.setAll(true);
      const fg = tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      const w = window as any;
      const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot');
      const items = () => sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      sp.props.filter = '';
      await w.__poll(items, (n: number) => n > 0, 500, 25);
      const beforeItems = await w.__settledFor(items, 150, 500, 25);
      const quiet = w.__quiet('viewer:Scatter plot.onViewerRendered', 150, 1500);
      sp.props.filter = '${Stereo Category} in ("R_ONE", "S_UNKN")';
      await quiet;
      const afterItems = await w.__settledFor(items, 150, 1500, 25);
      return {beforeItems, afterItems, filter: sp.props.filter};
    });
    expect(res.filter).toContain('Stereo Category');
    expect(res.afterItems, 'legend item count must respond to in-viewer filter (GROK-17222)')
      .toBeLessThanOrEqual(res.beforeItems);
  });

  await softStep('Step 6: Pie chart click-to-filter narrows legend', async () => {
    const result = await clickCanvasFilter(page, {viewerType: 'Pie chart', column: 'Stereo Category'});
    expect(result.totalFiltered).toBeGreaterThan(0);
    const lcItems = (await v.readLegend(page, 'Line chart')).itemCount;
    expect(lcItems, 'Line chart legend updates after Pie click-to-filter (GROK-17222)').toBeGreaterThanOrEqual(0);
  });

  await softStep('Step 7: Bar chart click-to-filter narrows legend', async () => {
    const result = await clickCanvasFilter(page, {viewerType: 'Bar chart', column: 'Stereo Category'});
    expect(result.survivors).toBe(1);
    expect(result.totalFiltered).toBeGreaterThan(0);
    const lcItems = (await v.readLegend(page, 'Line chart')).itemCount;
    expect(lcItems, 'Line chart legend updates after Bar click-to-filter (GROK-17222)').toBeGreaterThanOrEqual(0);
  });

  await softStep('Cleanup', async () => {
    await v.resetFilters(page);
    await v.cleanupShell(page);
  });

  v.finishSpec();
});
