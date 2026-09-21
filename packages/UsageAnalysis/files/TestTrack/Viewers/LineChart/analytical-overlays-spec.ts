/* ---
realizes: [linechart.cp.analytical-overlays]
--- */
import {expect, type Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

// The client-side half of the scenario; the layout round-trip (S2 steps 6-8) lives in
// line-chart-server-spec.ts on the server lane.
test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

async function setProps(page: Page, props: Record<string, any>) {
  await v.setViewerProps(page, 'Line chart', [{set: props, wait: 500}]);
}

async function formulaLinesCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const lc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Line chart') as any;
    try {
      const parsed = JSON.parse(lc.props.formulaLines || '[]');
      return Array.isArray(parsed) ? parsed.length : -1;
    } catch (e) {
      return -1;
    }
  });
}

test('Line Chart — Analytical Overlays', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error' && !isLocalBootNoise(m.text())) consoleErrors.push(m.text()); });
  const errorCount = () => pageErrors.length + consoleErrors.length;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'line-chart', 'Line-chart', 15_000, 'Line chart');

  await setProps(page, {xColumnName: 'CAST Idea ID', yColumnNames: ['Chemical Space X']});

  await softStep('S1 steps 1-3: enable regression line, no-error floor', async () => {
    const before = errorCount();
    await setProps(page, {showRegressionLine: true});
    expect(errorCount()).toBe(before);
  });

  await softStep('S1 steps 4-5: enable rolling average overlay, no split', async () => {
    const before = errorCount();
    await setProps(page, {showMovingAverageLine: true});
    expect(errorCount()).toBe(before);
  });

  await softStep('S1 steps 6-7: enable standard-deviation overlay, no split', async () => {
    const before = errorCount();
    await setProps(page, {showMovingAverageDeviation: true});
    expect(errorCount()).toBe(before);
  });

  await softStep('S1 steps 8-9: apply split, overlays render per category', async () => {
    const before = errorCount();
    await setProps(page, {splitColumnNames: ['Stereo Category']});
    expect(errorCount()).toBe(before);
  });

  await softStep('S1 steps 10-11: revert overlays and split', async () => {
    const before = errorCount();
    await setProps(page, {
      showRegressionLine: false,
      showMovingAverageLine: false,
      showMovingAverageDeviation: false,
      splitColumnNames: [],
    });
    expect(errorCount()).toBe(before);
  });

  await softStep('S2 steps 1-5: add a formula line and a formula band', async () => {
    const before = errorCount();
    await setProps(page, {formulaLines: JSON.stringify([
      {type: 'line', formula: '${Chemical Space X} = 500', title: 'const-line', color: '#FF0000'},
      {type: 'band', formula: '${Chemical Space X} in(400, 600)', title: 'const-band', color: '#00FF00'},
    ])});
    expect(await formulaLinesCount(page)).toBe(2);
    expect(errorCount()).toBe(before);
    await setProps(page, {formulaLines: ''});
  });

  await v.closeAllAndWait(page);
  v.finishSpec();
});
