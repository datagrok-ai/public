/* ---
sub_features_covered: [legend.corner.late-geometry-reanchor]
--- */
// The CI-only stale-anchor bug: the bar chart's corner offsets derive from geometry its own
// render settles (the centered content band arrives after bars aggregate), and the changed
// placement extra key used to be observed only at the next frame — which nothing scheduled
// on an idle viewer, so on a slow machine the legend sat 120px off on the axis until a
// resize. legendCommit now re-reads the key post-layout and requests the re-resolving frame
// itself. Pinned under 6x CPU throttling — the condition that separated CI from dev.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('A throttled bar chart corner legend re-anchors itself after late geometry', async ({page}) => {
  test.setTimeout(300_000);
  await loginToDatagrok(page);
  await v.openTable(page);

  const cdp = await page.context().newCDPSession(page);
  await cdp.send('Emulation.setCPUThrottlingRate', {rate: 6});

  await page.evaluate(() => {
    const w = (window as any).grok.shell.windows;
    w.showContextPanel = false; w.showProperties = false; w.showHelp = false;
    (window as any).grok.shell.tv.addViewer('Bar chart',
      {splitColumnName: 'Stereo Category', stackColumnName: 'Series'});
  });
  await v.waitForLegendIdle(page, 'Bar chart');
  await v.resizeViewer(page, 'Bar chart', 620, 760);

  const read = () => page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Bar chart');
    const root = x.root.querySelector('[name="legend"]') as HTMLElement;
    return {state: `${root.getAttribute('data-legend-mode')}/${root.getAttribute('data-legend-slot')}`,
      bottom: parseFloat(root.style.bottom || '0'), right: parseFloat(root.style.right || '0')};
  });

  const first = await read();
  await v.resizeViewer(page, 'Bar chart', 622, 760);
  const settled = await read();

  expect(first.state).toBe(settled.state);
  expect(Math.abs(first.bottom - settled.bottom),
    `stale anchor: first ${JSON.stringify(first)} vs settled ${JSON.stringify(settled)}`)
    .toBeLessThanOrEqual(2);
  expect(Math.abs(first.right - settled.right)).toBeLessThanOrEqual(2);

  await cdp.send('Emulation.setCPUThrottlingRate', {rate: 1});
  await v.cleanupShell(page);
  v.finishSpec('Late-geometry re-anchor failures');
});
