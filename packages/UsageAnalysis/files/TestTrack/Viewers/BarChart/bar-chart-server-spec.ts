/* ---
realizes: [barchart.cp.setup-and-interact]
--- */

import {test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitCol = 'Primary Series Name';

// The server lane of BOTH mixed bar-chart scenarios, in one file so the section loads its
// datasets once: the Data panel step of bar-chart.md and Scenario 3 of barchart-setup-interact.md
// are the only steps whose subject is a layout saved through dapi.layouts, closed and reloaded.
// Everything else about the bar chart is client behaviour and lives on the local lane.

test('Bar chart tests — data panel and colour-coding layout round-trip', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) ||
    /404 \(\)/.test(text) ||
    /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (e) => pageErrors.push(String(e)));

  await openDatagrok(page);
  await v.installEventWaits(page);
  await v.closeAllAndWait(page);

  await page.evaluate(async ({demog, spgi}) => {
    const w = window as any;
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.addTableView(await w.__readCsv(demog));
    await w.__tableReady();
    const df2 = await w.__readCsv(spgi);
    df2.name = 'SPGI';
    grok.shell.addTableView(df2);
    await w.__tableReady();
  }, {demog: demogPath, spgi: spgiPath});

  const showTable = (name: string | null) => page.evaluate(async (n: string | null) => {
    const w = window as any;
    const view = Array.from(grok.shell.views).find((x: any) =>
      x.type === 'TableView' && (n === null ? x.dataFrame.name !== 'SPGI' : x.dataFrame.name === n));
    await w.__settled('grok.events.onCurrentViewChanged', () => { grok.shell.v = view; }, 500);
  }, name);

  await softStep('Data panel', async () => {
    await showTable(null);

    const result = await page.evaluate(async () => {
      const w = window as any;
      const icon = document.querySelector('[name="icon-bar-chart"]') as HTMLElement;
      await w.__settled('grok.events.onViewerAdded', () => icon.click(), 1000);

      const bc = Array.from(grok.shell.tv.viewers).find((view: any) => view.type === 'Bar chart') as any;
      const r: any[] = [];

      for (const src of ['Filtered', 'Selected', 'All']) {
        await w.__settled('viewer:Bar chart.onViewerRendered', () => {
          bc.props.rowSource = src;
        }, 200);
        r.push(bc.props.rowSource);
      }
      bc.props.rowSource = 'All';

      const spgi = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      await w.__settled('viewer:Bar chart.onViewerRendered', () => { bc.dataFrame = spgi; }, 500);
      r.push(bc.dataFrame.name);

      await w.__settled('viewer:Bar chart.onViewerRendered', () => {
        bc.props.filter = '${CAST Idea ID} < 634835';
      }, 500);
      r.push(bc.props.filter);

      await w.__settled('viewer:Bar chart.onViewerRendered', () => {
        bc.props.colorColumnName = 'Chemical Space Y';
      }, 300);
      r.push(bc.props.colorColumnName);

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      const saved = await w.__findSaved(() => grok.dapi.layouts.find(layoutId), 1000);
      try {
        await w.__settled('grok.events.onViewerClosed', () => bc.close(), 500);
        await w.__settled('grok.events.onViewLayoutApplied', () => grok.shell.tv.loadLayout(saved), 3000);

        const bc2 = Array.from(grok.shell.tv.viewers).find((view: any) => view.type === 'Bar chart') as any;
        r.push(bc2 ? bc2.props.colorColumnName : 'NOT_RESTORED');
        r.push(bc2 ? bc2.props.filter : 'NOT_RESTORED');
        if (bc2) await w.__settled('grok.events.onViewerClosed', () => bc2.close(), 500);
      }
      finally {
        await grok.dapi.layouts.delete(saved);
      }
      return r;
    });
    expect(result.slice(0, 3)).toEqual(['Filtered', 'Selected', 'All']);
    expect(result[3]).toBe('SPGI');
    expect(result[4]).toBe('${CAST Idea ID} < 634835');
    expect(result[5]).toBe('Chemical Space Y');
    expect(result[6]).toBe('Chemical Space Y');
    expect(result[7]).toBe('${CAST Idea ID} < 634835');
  });

  await softStep('Scenario 3 Step 1: grid color coding on the Split column drives bar colors and survives a layout round-trip', async () => {
    const CEIL = 250;
    const FLOOR = 8000;

    await showTable('SPGI');
    await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');
    await v.waitForViewerRendered(page, 'Bar chart', 500);
    await v.setViewerProps(page, 'Bar chart', [{set: {
      splitColumnName: splitCol, valueColumnName: 'CAST Idea ID', valueAggrType: 'count',
    }, wait: 900}]);

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
    await v.waitForViewerQuiet(page, 'Bar chart', {gapMs: 300, capMs: 900});

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
    // the assertion is about pixels, so the wait is too: onViewerRendered lands before the paint
    const colorDelta = await v.waitForCanvasChange(page, 'Bar chart', {minDelta: FLOOR + 1, timeoutMs: 3000})
      .catch(async () => (await v.diffCanvasColors(page, 'Bar chart')).deltaPx);

    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id;
    });
    try {
      await page.evaluate(async ({split}) => {
        const w = window as any;
        const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
        await w.__settled('grok.events.onViewerClosed', () => bc.close(), 500);
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
        const w = window as any;
        w.__savedLayout = await w.__findSaved(() => grok.dapi.layouts.find(id), 1000);
        grok.shell.tv.loadLayout(w.__savedLayout);
      }, {id: layoutId});
      const rt = await v.pollValue(() => page.evaluate(({split}) => {
        const col2 = grok.shell.tv.dataFrame.col(split);
        const bc2 = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
        return {restored: JSON.parse(col2.tags['.color-coding-categorical'] ?? '{}'), reopened: !!bc2};
      }, {split: splitCol}),
      (x) => x.reopened && Object.keys(x.restored).length > 0, 3000, 150);

      expect(precheck).toBeGreaterThanOrEqual(0);
      expect(precheck).toBeLessThan(CEIL);
      expect(colorDelta).toBeGreaterThan(FLOOR);
      expect(clearedKeys).toBe(0);
      expect(rt.restored).toEqual(scheme);
      expect(rt.reopened).toBe(true);
    }
    finally {
      await page.evaluate(async ({id}) => {
        const w = window as any;
        const saved = w.__savedLayout ?? await grok.dapi.layouts.find(id);
        delete w.__savedLayout;
        if (saved) await grok.dapi.layouts.delete(saved);
      }, {id: layoutId});
    }
  });

  await softStep('No page errors', async () => {
    expect(pageErrors).toEqual([]);
    expect(consoleErrors).toEqual([]);
  });

  await v.closeAllAndWait(page);
  v.finishSpec();
});
