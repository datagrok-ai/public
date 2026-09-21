/* ---
realizes: [linechart.cp.analytical-overlays, linechart.cp.legend-color-and-persistence]
--- */
import {expect, type Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

// The server lane of the LineChart section: every step whose subject is a layout or a project
// surviving a round-trip through the server. The client-side ladders live in line-chart-spec.ts,
// analytical-overlays-spec.ts and legend-color-and-persistence-spec.ts on the local lane; the
// configured state is set up here directly through the API rather than re-driven.
test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitColumn = 'Stereo Category';
const baselineColors: Record<string, number> = {
  R_ONE: 0xFFFF0000,
  S_ABS: 0xFF00FF00,
  S_ACHIR: 0xFF0000FF,
  S_PART: 0xFFFFFF00,
  S_UNKN: 0xFFFF00FF,
};

async function lcProps(page: Page, ...names: string[]): Promise<Record<string, any>> {
  return page.evaluate((ns) => {
    const lc = (Array.from(grok.shell.tv.viewers) as any[]).find((x: any) => x.type === 'Line chart') as any;
    const out: Record<string, any> = {};
    for (const n of ns) out[n] = lc.props[n];
    return out;
  }, names);
}

async function lcSetProps(page: Page, props: Record<string, any>) {
  await v.setViewerProps(page, 'Line chart', [{set: props}], 300);
}

async function addLineChart(page: Page) {
  await v.addViewerByIcon(page, 'line-chart', 'Line-chart', 15_000, 'Line chart');
}

// awaiting layouts.save IS the completion signal, the same call the Save to Gallery menu makes
async function saveLayout(page: Page): Promise<string> {
  return page.evaluate(async () => String((await grok.dapi.layouts.save(grok.shell.tv.saveLayout())).id));
}

async function closeLineChart(page: Page) {
  await page.evaluate(() => {
    (Array.from(grok.shell.tv.viewers) as any[]).find((x: any) => x.type === 'Line chart')?.close();
  });
  await page.locator('[name="viewer-Line-chart"]').first().waitFor({state: 'detached', timeout: 2000});
}

async function loadLayout(page: Page, layoutId: string) {
  await page.evaluate(async (id) => {
    const w = window as any;
    const saved = await grok.dapi.layouts.find(id);
    const applied = w.__armed('grok.events.onViewLayoutApplied', 3000);
    grok.shell.tv.loadLayout(saved);
    await applied;
    await w.__poll(() => !!(Array.from(grok.shell.tv.viewers) as any[]).find((x: any) => x.type === 'Line chart')?.props,
      (ok: boolean) => ok, 3000, 50);
  }, layoutId);
}

// Detached, the way deleteProjectWithCleanup already is: the find+delete pair costs a second
// per layout on dev and this spec saves five of them, none of which it ever reads again. The
// worker fixture drains __pendingDeletes before it closes the page, so the server state still goes.
async function deleteLayout(page: Page, layoutId: string) {
  await page.evaluate((id) => {
    const w = window as any;
    w.__pendingDeletes = w.__pendingDeletes ?? [];
    w.__pendingDeletes.push((async () => {
      try {
        const saved = await grok.dapi.layouts.find(id);
        if (saved) await grok.dapi.layouts.delete(saved);
      } catch (_) {}
    })());
  }, layoutId);
}

async function readROneColor(page: Page): Promise<number | null> {
  return page.evaluate((col) => {
    const tv = grok.shell.tv;
    const cat = tv.dataFrame.col(col);
    for (let i = 0; i < tv.dataFrame.rowCount; i++)
      if (cat.get(i) === 'R_ONE') return cat.meta.colors.getColor(i, cat);
    return null;
  }, splitColumn);
}

async function formulaLinesCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const lc = (Array.from(grok.shell.tv.viewers) as any[]).find((x: any) => x.type === 'Line chart') as any;
    try {
      const parsed = JSON.parse(lc.props.formulaLines || '[]');
      return Array.isArray(parsed) ? parsed.length : -1;
    } catch (e) {
      return -1;
    }
  });
}

test('Line chart — layout and project persistence', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  const errorCount = () => pageErrors.length + consoleErrors.length;

  await openDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await addLineChart(page);

  await softStep('Layout save and restore', async () => {
    await lcSetProps(page, {
      xColumnName: 'STARTED', yColumnNames: ['AGE', 'HEIGHT'],
      splitColumnName: 'SEX', multiAxis: true, lineWidth: 3, interpolation: 'Spline',
    });
    const layoutId = await saveLayout(page);
    try {
      await closeLineChart(page);
      await loadLayout(page, layoutId);
      const restored = await lcProps(page, 'xColumnName', 'yColumnNames', 'splitColumnName', 'multiAxis', 'lineWidth', 'interpolation');
      expect(restored.xColumnName).toBe('STARTED');
      expect(restored.yColumnNames).toEqual(['AGE', 'HEIGHT']);
      expect(restored.splitColumnName).toBe('SEX');
      expect(restored.multiAxis).toBe(true);
      expect(restored.lineWidth).toBe(3);
      expect(restored.interpolation).toBe('Spline');
    } finally {
      await deleteLayout(page, layoutId);
    }
  });

  await softStep('Selection checkboxes — the four row-marker flags survive a layout round-trip', async () => {
    const combo = {showCurrentRowLine: true, showMouseOverCategory: false,
      showSelectedRows: false, showMouseOverRowLine: false};
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['HEIGHT'], splitColumnName: '', multiAxis: false,
      lineWidth: 1, interpolation: 'None', rowSource: 'All', ...combo});
    const layoutId = await saveLayout(page);
    try {
      await closeLineChart(page);
      await loadLayout(page, layoutId);
      expect(await lcProps(page, ...Object.keys(combo))).toEqual(combo);
    } finally {
      await deleteLayout(page, layoutId);
    }
    await lcSetProps(page, {showCurrentRowLine: false, showMouseOverCategory: true,
      showSelectedRows: true, showMouseOverRowLine: true});
  });

  await softStep('Data panel checkboxes — packCategories and multiAxis survive a layout round-trip', async () => {
    const combo = {packCategories: false, multiAxis: true};
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['AGE', 'HEIGHT'], ...combo});
    const layoutId = await saveLayout(page);
    try {
      await closeLineChart(page);
      await loadLayout(page, layoutId);
      expect(await lcProps(page, ...Object.keys(combo))).toEqual(combo);
    } finally {
      await deleteLayout(page, layoutId);
    }
    await lcSetProps(page, {packCategories: true, multiAxis: false});
  });

  await v.openTable(page, {path: spgiPath, semTypeTimeoutMs: 3000});
  await addLineChart(page);
  await page.evaluate((args) => {
    const lc = (Array.from(grok.shell.tv.viewers) as any[]).find((x: any) => x.type === 'Line chart') as any;
    lc.props.xColumnName = 'Chemical Space X';
    lc.props.yColumnNames = ['Chemical Space Y'];
    lc.props.splitColumnNames = [args.col];
    grok.shell.tv.dataFrame.col(args.col).meta.colors.setCategorical(args.colors);
  }, {col: splitColumn, colors: baselineColors});
  await v.waitForViewerRendered(page, 'Line chart', 1200);
  const legend = page.locator('[name="viewer-Line-chart"] [name="legend"]');
  await legend.waitFor({timeout: 10000});
  await expect(legend.locator('.d4-legend-item.d4-legend-text-item')).toHaveCount(5);
  expect(await readROneColor(page)).toBe(0xFFFF0000);

  await softStep('S2 steps 6-8: save layout, clear, reapply — formula lines restored', async () => {
    const before = errorCount();
    const formulaLinesSpec = JSON.stringify([
      {type: 'line', formula: '${Chemical Space X} = 500', title: 'const-line', color: '#FF0000'},
      {type: 'band', formula: '${Chemical Space X} in(400, 600)', title: 'const-band', color: '#00FF00'},
    ]);
    await lcSetProps(page, {formulaLines: formulaLinesSpec});
    expect(await formulaLinesCount(page)).toBe(2);

    const layoutId = await saveLayout(page);
    try {
      await lcSetProps(page, {formulaLines: ''});
      expect(await formulaLinesCount(page)).toBe(0);

      await loadLayout(page, layoutId);
      expect(await v.pollValue(() => formulaLinesCount(page), (n) => n === 2, 3000, 100)).toBe(2);
      expect(errorCount()).toBe(before);
    } finally {
      await lcSetProps(page, {formulaLines: ''});
      await deleteLayout(page, layoutId);
    }
  });

  await softStep('S2: category color persists through layout round-trip (GROK-17278)', async () => {
    const before = errorCount();
    await page.evaluate((col) => {
      grok.shell.tv.dataFrame.col(col).meta.colors.setCategorical({R_ONE: 0xFF00AAFF});
    }, splitColumn);
    const expected = await v.pollValue(() => readROneColor(page), (c) => c === 0xFF00AAFF, 400, 50);

    const layoutId = await saveLayout(page);
    try {
      await page.evaluate((col) => {
        grok.shell.tv.dataFrame.col(col).meta.colors.setCategorical({});
      }, splitColumn);
      const cleared = await v.pollValue(() => readROneColor(page), (c) => c !== expected, 400, 50);
      await loadLayout(page, layoutId);
      const restored = await v.pollValue(() => readROneColor(page), (c) => c === expected, 3000, 150);
      expect(cleared).not.toBe(expected);
      expect(restored).toBe(expected);
      expect(errorCount()).toBe(before);
    } finally {
      await deleteLayout(page, layoutId);
    }
  });

  await softStep('S2 Steps 9-13: color AND markers legend survive a project save/close/reopen via the SAVE button (GROK-17278, GROK-19825)', async () => {
    await page.evaluate((col) => {
      grok.shell.tv.dataFrame.col(col).meta.colors.setCategorical({R_ONE: 0xFF00AAFF});
    }, splitColumn);
    const expected = await v.pollValue(() => readROneColor(page), (c) => c === 0xFF00AAFF, 500, 50);

    // every assertion below is @Prop / column-tag state, so the API save path applies (see
    // helpers-registry.yaml for the saveProjectViaApi / saveProjectViaUI boundary)
    let projId: string | null = null;
    try {
      projId = (await saveProjectViaApi(page, 'zz-linechart-color-persist-' + Date.now())).projectId;
      expect(!!projId).toBe(true);

      await v.closeAllAndWait(page);
      await page.evaluate(async (id) => {
        const full = await grok.dapi.projects.find(id);
        await full.open();
      }, projId);
      const lcRestored = await v.pollValue(() => page.evaluate(() => {
        const tv = grok.shell.tv;
        return !!tv && Array.from(tv.viewers).some((x: any) => x.type === 'Line chart');
      }), (restored) => restored, 4500, 150);
      // the reopened view builds its legend a frame after the viewer itself lands, so a
      // one-shot read here races the render and sees no legend at all
      const legendDom = await v.pollValue(() => page.evaluate(() => {
        const el = document.querySelector('[name="viewer-Line-chart"] [name="legend"]') as HTMLElement | null;
        return !!el && getComputedStyle(el).display !== 'none';
      }), (visible) => visible, 4500, 150);
      const color = await page.evaluate((col) => {
        const tv = grok.shell.tv;
        if (!tv) return null;
        const cat = tv.dataFrame.col(col);
        for (let i = 0; i < tv.dataFrame.rowCount; i++)
          if (cat.get(i) === 'R_ONE') return cat.meta.colors.getColor(i, cat);
        return null;
      }, splitColumn);

      expect(lcRestored).toBe(true);
      expect(color).toBe(expected);
      expect(legendDom).toBe(true);
    } finally {
      if (projId) await deleteProjectWithCleanup(page, {projectId: projId});
    }
  });

  // cleanupShell's own 500ms sleep, replaced by the condition it stood for: closeAllAndWait
  // waits for the table views to actually go away
  await page.evaluate(() => {
    const col = (window as any).grok.shell.tv?.dataFrame.col('Stereo Category');
    if (col) {
      delete col.tags['.color-coding-categorical'];
      delete col.tags['.color-coding-type'];
    }
  });
  await v.closeAllAndWait(page);
  v.finishSpec();
});
