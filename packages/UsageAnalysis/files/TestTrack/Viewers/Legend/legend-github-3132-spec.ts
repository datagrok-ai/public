/* ---
sub_features_covered: [legend.item.color-picker, legend.allow-item-coloring]
related_bugs: [github-3132]
coverage_type: regression
bug_url: https://github.com/datagrok-ai/public/issues/3132
--- */
// Fix invariant: each color change made through the legend persists independently;
// changing color B must NOT reset previously-applied color A. Full prose moved to
// legend-github-3132.md.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

test('github-3132: sequential legend color changes persist independently', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await v.openTableForLegend(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Histogram', 'Scatter plot']});

  // Step 3 (bug repro): change FIRST colour through the legend (R_ONE → red).
  await softStep('Step 3: change R_ONE to red via legend picker', async () => {
    await v.changeLegendItemColor(page, {
      viewerType: 'Histogram',
      category: 'R_ONE',
      rgb: [255, 0, 0],
      hex: '#FF0000',
      column: 'Stereo Category',
    });
  });

  // Step 4 (bug repro): change SECOND colour (S_UNKN → green). Bug invariant
  // under test: this MUST NOT reset R_ONE to default. additive: R_ONE retained.
  await softStep('Step 4: change S_UNKN to green via legend picker', async () => {
    await v.changeLegendItemColor(page, {
      viewerType: 'Histogram',
      category: 'S_UNKN',
      rgb: [0, 255, 0],
      hex: '#00FF00',
      column: 'Stereo Category',
      additive: {'R_ONE': '#FF0000', 'S_UNKN': '#00FF00'},
    });
  });

  // Step 4 (bug invariant assertion): verify R_ONE STILL holds the previously-applied red.
  await softStep('Step 4 invariant: R_ONE retains red after S_UNKN change', async () => {
    const r = await page.evaluate(() => {
      const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
      const t = JSON.parse(col.tags['.color-coding-categorical'] ?? '{}');
      return {
        rOne: String(t['R_ONE'] ?? '').toLowerCase(),
        sUnkn: String(t['S_UNKN'] ?? '').toLowerCase(),
      };
    });
    expect(r.rOne, 'R_ONE color must persist after S_UNKN change (github-3132)').toBe('#ff0000');
    expect(r.sUnkn).toBe('#00ff00');
  });

  await softStep('Cleanup', async () => { await v.cleanupShell(page, {clearStereoCategoryColorCoding: true}); });

  v.finishSpec();
});
