/* ---
sub_features_covered: [legend.allow-item-coloring, legend.item.color-picker]
--- */
// github-3132: each legend color change persists independently; changing B must not reset A.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import * as v from '../helpers/viewers';

test.use(specTestOptions);

test('github-3132: sequential legend color changes persist independently', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Histogram', 'Scatter plot']});

  await softStep('Step 3: change R_ONE to red via legend picker', async () => {
    await v.changeLegendItemColor(page, {
      viewerType: 'Histogram',
      category: 'R_ONE',
      rgb: [255, 0, 0],
      hex: '#FF0000',
      column: 'Stereo Category',
    });
  });

  // Changing S_UNKN must not reset R_ONE — additive retains R_ONE.
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
