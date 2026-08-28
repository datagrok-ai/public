/* ---
realizes: [barchart.cp.setup-and-interact]
--- */

import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitCol = 'Primary Series Name';

// The server lane of the mixed setup-and-interact scenario: Scenario 3's colour scheme is
// persisted through dapi.layouts and read back after a close + reload, which local mode does
// not serve. Scenarios 1 and 2 are pure client interaction and stay in barchart-setup-interact-spec.ts.

test('Bar Chart — Setup and core interaction (colour coding layout round-trip)', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 4000});
  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  await v.setViewerProps(page, 'Bar chart', [{set: {
    splitColumnName: splitCol, valueColumnName: 'CAST Idea ID', valueAggrType: 'count',
  }, wait: 900}]);

  await softStep('Scenario 3 Step 1: grid color coding on the Split column drives bar colors and survives a layout round-trip', async () => {

    const CEIL = 250;
    const FLOOR = 8000;

    await page.evaluate(({split}) => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      df.selection.setAll(false);
      const col = df.col(split);
      delete col.tags['.color-coding-categorical'];
      delete col.tags['.color-coding-type'];
      col.meta.colors.setCategorical({});
    }, {split: splitCol});

    await v.waitForViewerRendered(page, 'Bar chart', 500);

    await v.waitForCanvasQuiet(page, 'Bar chart', {timeoutMs: 900, optional: true});

    await v.snapshotCanvasColors(page, 'Bar chart');
    const precheck = (await v.diffCanvasColors(page, 'Bar chart')).deltaPx;

    const scheme = await page.evaluate(({split}) => {
      const col = grok.shell.tv.dataFrame.col(split);
      const palette = ['#ff0000', '#00ff00', '#0000ff', '#ff00ff', '#00ffff', '#ffff00', '#ff8000', '#8000ff'];
      const s: Record<string, string> = {};
      col.categories.forEach((c: string, i: number) => { s[c] = palette[i % palette.length]; });
      col.meta.colors.setCategorical(s);
      for (const vw of grok.shell.tv.viewers) if (vw.type !== 'Grid') try { vw.invalidate?.(); } catch (_) {}
      return s;
    }, {split: splitCol});
    // The assertion is about pixels, so the wait is too: waitForViewerRendered resolves on any
    // render in the last 400ms, which can land before the recolour is painted.
    const colorDelta = await v.waitForCanvasChange(page, 'Bar chart', {minDelta: 1, timeoutMs: 3000})
      .catch(() => 0);

    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id;
    });

    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.close();
    });
    await v.pollValue(
      () => page.evaluate(() => !Array.from(grok.shell.tv.viewers).some((x: any) => x.type === 'Bar chart')),
      (closed) => closed, 500, 100);
    await page.evaluate(({split}) => {
      const col = grok.shell.tv.dataFrame.col(split);
      delete col.tags['.color-coding-categorical'];
      delete col.tags['.color-coding-type'];
      col.meta.colors.setCategorical({});
    }, {split: splitCol});
    const clearedKeys = await v.pollValue(() => page.evaluate(({split}) => {
      const col = grok.shell.tv.dataFrame.col(split);
      return Object.keys(JSON.parse(col.tags['.color-coding-categorical'] ?? '{}')).length;
    }, {split: splitCol}), (n) => n === 0, 400, 100);
    await page.evaluate(async ({id}) => {
      const saved = await grok.dapi.layouts.find(id);
      (window as any).__savedLayout = saved;
      grok.shell.tv.loadLayout(saved);
    }, {id: layoutId});
    const rt = await v.pollValue(() => page.evaluate(({split}) => {
      const col2 = grok.shell.tv.dataFrame.col(split);
      const bc2 = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      return {restored: JSON.parse(col2.tags['.color-coding-categorical'] ?? '{}'), reopened: !!bc2};
    }, {split: splitCol}),
    (x) => x.reopened && Object.keys(x.restored).length > 0, 3000, 150);
    await page.evaluate(async () => {
      const w = window as any;
      await grok.dapi.layouts.delete(w.__savedLayout);
      delete w.__savedLayout;
    });

    expect(precheck).toBeGreaterThanOrEqual(0); 
    expect(precheck).toBeLessThan(CEIL);
    expect(colorDelta).toBeGreaterThan(FLOOR);
    expect(clearedKeys).toBe(0);
    expect(rt.restored).toEqual(scheme);
    expect(rt.reopened).toBe(true);
  });

  await softStep('No page errors', async () => {
    expect(pageErrors).toEqual([]);
    expect(consoleErrors).toEqual([]);
  });

  v.finishSpec();
});
