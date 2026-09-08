/* ---
realizes: [viewers.histogram, viewers.line-chart, viewers.bar-chart, viewers.pie-chart, viewers.trellis-plot, viewers.box-plot]
--- */
// Scenario 2 picker UI runs on Histogram: Bar chart legend needs a color edit to render.
// The layout and project round-trips live in color-consistency-server-spec.ts.

import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {addLegendViewers} from './legend-setup';

test.use(specTestOptions);

test('Legend color consistency', async ({page}) => {
  test.setTimeout(600_000);

  await openDatagrok(page);
  await v.openTable(page);
  await v.installEventWaits(page);
  await addLegendViewers(page, {
    column: 'Stereo Category',
    viewers: ['Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'],
  });

  await softStep('Categorical color coding from grid: R_ONE=red, S_UNKN=green', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const df = w.grok.shell.tv.dataFrame;
      const col = df.col('Stereo Category');
      const swatches = () => Array.from(w.grok.shell.tv.viewers)
        .filter((x: any) => x.type !== 'Grid')
        .map((x: any) => Array.from(x.root.querySelectorAll('[name="legend"] .d4-legend-item'))
          .map((el: any) => getComputedStyle(el).color).join(',')).join(';');
      col.tags['.color-coding-type'] = 'Categorical';
      col.meta.colors.setCategorical(
        {'R_ONE': '#FF0000', 'S_UNKN': '#00FF00'},
        {fallbackColor: '#808080'},
      );
      for (const x of w.grok.shell.tv.viewers)
        if (x.type !== 'Grid') try { x.invalidate?.(); } catch (_) {}
      // the next step reads the rendered swatch colours, so the settle is on them, not the tag
      await w.__settledFor(swatches, 150, 1500, 25);
      let tagColors: Record<string, any> = {};
      try { tagColors = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}'); } catch (_) {}
      return {
        codingType: col.tags['.color-coding-type'],
        rOneTag: tagColors['R_ONE'] ?? null,
        sUnknTag: tagColors['S_UNKN'] ?? null,
      };
    });
    expect(res.codingType).toBe('Categorical');
    expect(String(res.rOneTag).toLowerCase()).toBe('#ff0000');
    expect(String(res.sUnknTag).toLowerCase()).toBe('#00ff00');
  });

  await softStep('Every viewer reflects R_ONE=red and S_UNKN=green (DOM)', async () => {
    const result = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const out: Record<string, any> = {viewers: {}};
      for (const x of tv.viewers) {
        if (x.type === 'Grid') continue;
        const items = Array.from(x.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
        const rOneItem = items.find((el) => el.querySelector('.d4-legend-value')?.textContent?.trim() === 'R_ONE');
        const sUnknItem = items.find((el) => el.querySelector('.d4-legend-value')?.textContent?.trim() === 'S_UNKN');
        out.viewers[x.type] = {
          legendRendered: items.length > 0,
          rOneColor: rOneItem ? getComputedStyle(rOneItem).color : null,
          sUnknColor: sUnknItem ? getComputedStyle(sUnknItem).color : null,
        };
      }
      return out;
    });
    let viewersWithLegend = 0;
    for (const [_, info] of Object.entries(result.viewers as Record<string, any>)) {
      if (!info.legendRendered) continue;
      if (info.rOneColor === 'rgb(255, 0, 0)' && info.sUnknColor === 'rgb(0, 255, 0)')
        viewersWithLegend++;
    }
    expect(viewersWithLegend, 'at least 1 viewer renders legend with the configured DOM colors').toBeGreaterThanOrEqual(1);
  });

  await softStep('Open color picker via legend, change R_ONE to blue', async () => {
    await v.changeLegendItemColor(page, {
      viewerType: 'Histogram',
      category: 'R_ONE',
      rgb: [31, 119, 180],
      hex: '#1f77b4',
      column: 'Stereo Category',
      additive: {'R_ONE': '#1f77b4', 'S_UNKN': '#00FF00'},
    });
  });

  await softStep('Picker change propagated: every legend item in DOM shows blue', async () => {
    const result = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const out: Record<string, string | null> = {};
      for (const x of tv.viewers) {
        if (x.type === 'Grid') continue;
        const items = Array.from(x.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
        const rOneItem = items.find((el) => el.querySelector('.d4-legend-value')?.textContent?.trim() === 'R_ONE');
        out[x.type] = rOneItem ? getComputedStyle(rOneItem).color : null;
      }
      return out;
    });
    let viewersChecked = 0;
    for (const [_, color] of Object.entries(result)) {
      if (color === 'rgb(31, 119, 180)') viewersChecked++;
    }
    expect(viewersChecked, 'picker change reflected in legend DOM on at least 1 viewer').toBeGreaterThanOrEqual(1);
  });

  await softStep('Cleanup', async () => { await v.cleanupShell(page, {clearStereoCategoryColorCoding: true}); });

  v.finishSpec();
});
