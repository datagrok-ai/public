/* ---
realizes: [linechart.cp.legend-color-and-persistence]
--- */
import {expect, type Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

// The client-side half of the scenario (S1: the legend click filter); the layout and project
// round-trips of S2 live in line-chart-server-spec.ts on the server lane.
test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitColumn = 'Stereo Category';
const baselineColors: Record<string, number> = {
  R_ONE: 0xFFFF0000,
  S_ABS: 0xFF00FF00,
  S_ACHIR: 0xFF0000FF,
  S_PART: 0xFFFFFF00,
  S_UNKN: 0xFFFF00FF,
};

async function readCategoryColors(page: Page): Promise<Record<string, number>> {
  return page.evaluate((col) => {
    const df = grok.shell.tv.dataFrame;
    const cat = df.col(col);
    const colors = cat.meta.colors;
    const out: Record<string, number> = {};
    const seen: Record<string, boolean> = {};
    for (let i = 0; i < df.rowCount; i++) {
      const v = cat.get(i);
      if (v && !seen[v]) { seen[v] = true; out[v] = colors.getColor(i, cat); }
    }
    return out;
  }, splitColumn);
}

async function lcFilterCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const lc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Line chart') as any;
    return lc.filter.trueCount;
  });
}

test('Line Chart — legend filter-color and layout persistence', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error' && !isLocalBootNoise(m.text())) consoleErrors.push(m.text()); });
  const errorCount = () => pageErrors.length + consoleErrors.length;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'line-chart', 'Line-chart', 15_000, 'Line chart');

  await page.evaluate((args) => {
    const lc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Line chart') as any;
    lc.props.xColumnName = 'Chemical Space X';
    lc.props.yColumnNames = ['Chemical Space Y'];
    lc.props.splitColumnNames = [args.col];
    grok.shell.tv.dataFrame.col(args.col).meta.colors.setCategorical(args.colors);
  }, {col: splitColumn, colors: baselineColors});
  await v.waitForViewerRendered(page, 'Line chart', 1200);

  const legend = page.locator('[name="viewer-Line-chart"] [name="legend"]');
  await legend.waitFor({timeout: 10000});
  await expect(legend.locator('.d4-legend-item.d4-legend-text-item')).toHaveCount(5);

  const baselineFilter = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  const colorBefore = await readCategoryColors(page);
  expect(colorBefore.R_ONE).toBe(0xFFFF0000);

  await softStep('S1: legend click filters + preserves remaining line colors', async () => {
    const before = errorCount();
    await legend.locator('.d4-legend-value', {hasText: /^S_ABS$/}).click();
    const filtered = await v.pollValue(() => lcFilterCount(page), (n) => n < baselineFilter, 700, 100);
    expect(filtered).toBeLessThan(baselineFilter);

    const colorAfter = await readCategoryColors(page);
    expect(colorAfter.R_ONE).toBe(colorBefore.R_ONE);
    expect(colorAfter.R_ONE).not.toBe(baselineColors.S_ABS);
    expect(errorCount()).toBe(before);
  });

  await softStep('S1: re-click legend category resets filter to full count', async () => {
    const before = errorCount();
    await legend.locator('.d4-legend-value', {hasText: /^S_ABS$/}).click();
    await v.pollValue(() => lcFilterCount(page), (n) => n === baselineFilter, 700, 100);
    expect(await lcFilterCount(page)).toBe(baselineFilter);
    expect(errorCount()).toBe(before);
  });

  await v.cleanupShell(page, {clearStereoCategoryColorCoding: true});
  v.finishSpec();
});
