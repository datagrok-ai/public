/* ---
realizes: [viewers.line-chart]
--- */
import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

// The layout and project round-trips live in line-chart-server-spec.ts.
test.use(specTestOptions);

test('Line chart legend', async ({page}) => {
  test.setTimeout(900_000);

  await openDatagrok(page);
  await v.openTable(page);
  await v.installEventWaits(page);

  await softStep('Setup + Sc1 steps 1-2: add Line chart, Split=Series → categorical legend', async () => {
    const items = await page.evaluate(async () => {
      const w = window as any;
      const tv = (window as any).grok.shell.tv;
      tv.addViewer('Line chart');
      await w.__poll(() => (window as any).grok.shell.tv.viewers.filter((x: any) => x.type === 'Line chart').length,
        (c: number) => c > 0, 1000);
      const lc = tv.viewers.find((x: any) => x.type === 'Line chart');
      lc.props.splitColumnName = 'Series';
      try { lc.props.legendVisibility = 'Always'; } catch (_) {}
      let prev = -1;
      await w.__poll(() => lc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
        (c: number) => { const settled = c > 0 && c === prev; prev = c; return settled; }, 1500);
      return lc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
    });
    expect(items).toBeGreaterThan(0);
  });

  await softStep('Sc1 verification: legend hover surface available (real DOM)', async () => {
    const legend = page.locator('[name="viewer-Line-chart"] [name="legend"]').first();
    if (await legend.count() === 0) {
      const alt = page.locator('[name="viewer-Line chart"] [name="legend"]').first();
      await alt.waitFor({timeout: 10000});
      await alt.hover();
    } else {
      await legend.waitFor({timeout: 10000});
      await legend.hover();
    }
  });

  await softStep('Sc2 steps 1-2: enable Multi Axis (multiAxis=true)', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const lc = w.grok.shell.tv.viewers.find((x: any) => x.type === 'Line chart');
      const quiet = w.__quiet('viewer:Line chart.onViewerRendered', 150, 1500);
      try { lc.props.multiAxis = true; } catch (e) { return {multiAxis: false, err: String(e)}; }
      await quiet;
      return {multiAxis: lc.props.multiAxis};
    });
    expect(res.multiAxis).toBe(true);
  });

  await softStep('Sc4 steps 1-2: yColumnNames = [Average Mass, TPSA]', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const lc = w.grok.shell.tv.viewers.find((x: any) => x.type === 'Line chart');
      const items = () => lc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      const quiet = w.__quiet('viewer:Line chart.onViewerRendered', 150, 2000);
      lc.props.yColumnNames = ['Average Mass', 'TPSA'];
      await quiet;
      return {yCols: lc.props.yColumnNames, totalItems: await w.__settledFor(items, 150, 2000, 25)};
    });
    expect(res.yCols).toEqual(['Average Mass', 'TPSA']);
    expect(res.totalItems).toBeGreaterThan(0);
  });

  await softStep('Sc4 step 3: replace Y column → NIBR logP', async () => {
    const res = await page.evaluate(async () => {
      const w = window as any;
      const lc = w.grok.shell.tv.viewers.find((x: any) => x.type === 'Line chart');
      const items = () => lc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      const quiet = w.__quiet('viewer:Line chart.onViewerRendered', 150, 1800);
      lc.props.yColumnNames = ['Average Mass', 'NIBR logP'];
      await quiet;
      return {yCols: lc.props.yColumnNames, items: await w.__settledFor(items, 150, 1800, 25)};
    });
    expect(res.yCols).toEqual(['Average Mass', 'NIBR logP']);
  });

  await softStep('Cleanup', async () => { await v.cleanupShell(page); });

  v.finishSpec();
});
