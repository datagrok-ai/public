/* ---
sub_features_covered: [legend.line-chart.auto-markers-with-column]
--- */
// showMarkers=Auto hides markers above a drawn-row-count threshold, which used to apply even
// when the user explicitly picked a marker column: the marker legend popped in when a legend
// selection or filter shrank the row count, and vanished again on deselect. An explicit marker
// column now always shows markers in Auto — the heuristic only governs the no-column default.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Auto markers with an explicit marker column stay visible (demog)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
  });
  await page.waitForFunction(() =>
    (window as any).grok?.shell?.tv?.dataFrame?.rowCount > 0, null, {timeout: 30000});

  await page.evaluate(async () => {
    const tv = (window as any).grok.shell.tv;
    tv.addViewer('Line chart');
    await new Promise((r) => setTimeout(r, 1200));
    const lc = tv.viewers.find((x: any) => x.type === 'Line chart');
    lc.setOptions({yColumnNames: ['AGE'], splitColumnNames: ['RACE'], markersColumnName: 'SEX'});
  });
  await v.waitForLegendIdle(page, 'Line chart');
  await v.resizeViewer(page, 'Line chart', 900, 500);

  const legendLabels = () => page.evaluate(() => {
    const root = document.querySelector('[name="legend"]') as HTMLElement;
    return (Array.from(root?.querySelectorAll('[name="legend-item"]') ?? []) as HTMLElement[])
      .map((i) => (i.textContent || '').trim());
  });
  const hasMarkerItems = (labels: string[]) =>
    labels.some((l) => l.startsWith('F')) && labels.some((l) => l.startsWith('M'));

  await softStep('Marker legend is shown from the start despite 10k drawn rows', async () => {
    const labels = await legendLabels();
    expect(hasMarkerItems(labels), `marker items missing in: ${labels.join(', ')}`).toBe(true);
  });

  await softStep('Filtering below the density threshold does not change the legend set', async () => {
    const before = await legendLabels();
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.dataFrame.rows.filter((r: any) => r.idx < 300);
    });
    await v.waitForLegendIdle(page, 'Line chart');
    const filtered = await legendLabels();
    expect(hasMarkerItems(filtered)).toBe(true);
    await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      tv.dataFrame.filter.setAll(true);
    });
    await v.waitForLegendIdle(page, 'Line chart');
    const reset = await legendLabels();
    expect(hasMarkerItems(reset), 'marker legend must not vanish when the filter resets').toBe(true);
    expect(reset).toEqual(before);
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Auto-marker visibility failures');
});
