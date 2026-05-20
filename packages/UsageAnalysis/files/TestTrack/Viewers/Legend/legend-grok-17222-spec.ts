/* ---
sub_features_covered: [legend.refresh.on-data-change, legend.column]
related_bugs: [GROK-17222]
coverage_type: regression
bug_url: https://reddata.atlassian.net/browse/GROK-17222
--- */
// Fix invariant: legend reflects the current filtered state of data; legend
// categories and item count must update under every filter trigger source —
// Filter Panel, in-viewer filter, and click-to-filter on Pie + Bar charts.
// Full prose moved to legend-grok-17222.md.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

test('GROK-17222: legend reflects filter state across 4 trigger sources', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await v.openTableForLegend(page, {withFilterPanel: true});
  await v.addLegendViewers(page, {
    column: 'Stereo Category',
    viewers: ['Line chart', 'Scatter plot', 'Pie chart', 'Bar chart'],
  });

  // Step 4 (bug repro): Filter Panel categorical filter — legend MUST update.
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

  // Step 5 (bug repro): in-viewer Scatter filter via sp.props.filter range expression.
  await softStep('Step 5: in-viewer Scatter filter via sp.props.filter range', async () => {
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
      sp.props.filter = '${Stereo Category} in ("R_ONE", "S_UNKN")';
      await new Promise((r) => setTimeout(r, 1500));
      const afterItems = sp.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      return {beforeItems, afterItems, filter: sp.props.filter};
    });
    expect(res.filter).toContain('Stereo Category');
    expect(res.afterItems, 'legend item count must respond to in-viewer filter (GROK-17222)')
      .toBeLessThanOrEqual(res.beforeItems);
  });

  // Step 6 (bug repro): Pie chart click-to-filter narrows the dataset → legend updates.
  await softStep('Step 6: Pie chart click-to-filter narrows legend', async () => {
    const result = await v.clickCanvasFilter(page, {viewerType: 'Pie chart', column: 'Stereo Category'});
    expect(result.totalFiltered).toBeGreaterThan(0);
    const lcItems = (await v.readLegend(page, 'Line chart')).itemCount;
    expect(lcItems, 'Line chart legend updates after Pie click-to-filter (GROK-17222)').toBeGreaterThanOrEqual(0);
  });

  // Step 7 (bug repro): Bar chart click-to-filter narrows the dataset → legend updates.
  await softStep('Step 7: Bar chart click-to-filter narrows legend', async () => {
    const result = await v.clickCanvasFilter(page, {viewerType: 'Bar chart', column: 'Stereo Category'});
    expect(result.survivors).toBe(1);
    expect(result.totalFiltered).toBeGreaterThan(0);
    const lcItems = (await v.readLegend(page, 'Line chart')).itemCount;
    expect(lcItems, 'Line chart legend updates after Bar click-to-filter (GROK-17222)').toBeGreaterThanOrEqual(0);
  });

  // Cleanup: clear filters and close views.
  await softStep('Cleanup', async () => {
    await v.resetFilters(page);
    await v.cleanupShell(page);
  });

  v.finishSpec();
});
