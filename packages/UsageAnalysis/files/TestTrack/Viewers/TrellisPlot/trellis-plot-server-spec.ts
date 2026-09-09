/* ---
realizes: [trellisplot.cp.split-and-pick-inner, trellisplot.int.undo-redo-viewer-lifecycle]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

// The server lane of the trellis scenario: the layout / project round-trip, the Multi Curve inner
// viewer (curves.csv has no local copy and the viewer is a package one) and To Script (the menu
// group is package-contributed and absent from the local client). The client-side ladder is
// trellis-plot-spec.ts on the local lane.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const curvesPath = 'System:DemoFiles/curves.csv';

function axisViewportCount(n: number, oneColumnOnly: boolean): number {
  return Math.min(oneColumnOnly ? 5 : (n < 5 * 1.5 ? n : 5), n);
}

async function installTrellisWaits(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__tpApply) return;
    w.__tp = () => Array.from(w.grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot');
    w.__tpCells = () => document.querySelectorAll('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').length;
    // A trellis prop change repaints once, synchronously inside the setter (measured 2026-09-03: every
    // settle armed after the act saw 0 renders and burned its cap), so the subscription is armed first.
    w.__tpRendered = (act: () => any, cap = 1500, gap = 250, viewer?: any) => new Promise<number>((resolve) => {
      const vw = viewer ?? w.__tp();
      let seen = 0;
      let timer: any = null;
      let sub: any = null;
      const done = () => {
        clearTimeout(timer);
        clearTimeout(capT);
        try { sub?.unsubscribe(); } catch (_) {}
        resolve(seen);
      };
      const capT = setTimeout(done, cap);
      try { sub = vw.onViewerRendered.subscribe(() => { seen++; clearTimeout(timer); timer = setTimeout(done, gap); }); }
      catch (_) {}
      act();
    });
    w.__tpApply = async (act: () => any, read?: () => any, cap = 1500, viewer?: any) => {
      const t0 = Date.now();
      await w.__tpRendered(act, cap, 250, viewer);
      if (!read) return undefined;
      return w.__settledFor(read, 150, Math.max(50, cap - (Date.now() - t0)), 25);
    };
  });
}

test('Trellis plot — layout round-trip, project save dialog, Multi Curve inner viewer, To Script', async ({page}) => {
  test.setTimeout(240_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  const onPageError = (e: Error) => { pageErrors.push(String(e)); };
  const onConsole = (m: any) => { if (m.type() === 'error') consoleErrors.push(m.text()); };
  page.on('pageerror', onPageError);
  page.on('console', onConsole);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath});
  await installTrellisWaits(page);

  const setup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    return {rowCount: df.rowCount, sex: df.col('SEX').categories.length, race: df.col('RACE').categories.length};
  });
  expect(setup).toEqual({rowCount: 5850, sex: 2, race: 4});
  const canonicalCellCount = axisViewportCount(setup.sex, false) * axisViewportCount(setup.race, false);

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);
  await page.evaluate(() => {
    const w = window as any;
    return w.__tpApply(() => {
      const tp = w.__tp();
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
    }, w.__tpCells, 1500);
  });
  await expect(page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell')).toHaveCount(canonicalCellCount);

  await softStep('Layout and Project save/restore', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const r: any = {};
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      r.viewersAtSave = Array.from(grok.shell.tv.viewers).map((x: any) => x.type).sort();

      const beforeAddCount = grok.shell.tv.viewers.length;
      grok.shell.tv.addViewer('Histogram');
      grok.shell.tv.addViewer('Bar chart');
      await w.__poll(() => grok.shell.tv.viewers.length, (n: number) => n >= beforeAddCount + 2, 1200, 60);
      r.viewersBefore = Array.from(grok.shell.tv.viewers).map((x: any) => x.type).sort();

      const saved = await w.__findSaved(() => grok.dapi.layouts.find(layoutId), 2000);
      await w.__settled('grok.events.onViewLayoutApplied', () => grok.shell.tv.loadLayout(saved), 3000);
      r.viewersAfter = Array.from(grok.shell.tv.viewers).map((x: any) => x.type).sort();

      await grok.dapi.layouts.delete(saved);
      return r;
    });
    expect(result.viewersBefore.length).toBeGreaterThan(result.viewersAfter.length);
    expect(result.viewersAfter).toEqual(result.viewersAtSave);
    expect(result.viewersAfter).toContain('Trellis plot');

    const errBeforeSave = consoleErrors.length;
    const pageErrBeforeSave = pageErrors.length;
    const saveBtn = page.locator('[name="button-Save"]').first();
    await saveBtn.waitFor({state: 'visible', timeout: 15000});
    await saveBtn.click();

    await page.locator('.d4-dialog').first().waitFor({state: 'visible', timeout: 2000}).catch(() => {});
    const dialogOpen = await page.evaluate(() => !!document.querySelector('.d4-dialog'));
    const errAfterSave = consoleErrors.length;
    const pageErrAfterSave = pageErrors.length;

    await page.locator('[name="button-CANCEL"]').first().click({timeout: 5000}).catch(() => {});
    await page.locator('.d4-dialog').first().waitFor({state: 'detached', timeout: 500}).catch(() => {});
    expect(dialogOpen).toBe(true);
    expect(errAfterSave).toBe(errBeforeSave);
    expect(pageErrAfterSave).toBe(pageErrBeforeSave);
  });

  await softStep('Multi Curve inner viewer', async () => {
    const result = await page.evaluate(async (cPath) => {
      const w = window as any;
      const demogView = grok.shell.tv;
      const demogName = demogView.dataFrame.name;
      const tp = w.__tp();
      if (!tp) return {error: 'Trellis plot not found on demog view'};
      await w.__tpApply(() => {
        tp.props.viewerType = 'Scatter plot';
        tp.props.xColumnNames = ['SEX'];
        tp.props.yColumnNames = ['RACE'];
      }, w.__tpCells, 1200);

      const dfCurves = await grok.dapi.files.readCsv(cPath);
      grok.shell.addTableView(dfCurves);
      await w.__poll(() => grok.shell.tv?.dataFrame?.name, (n: string) => n === dfCurves.name, 1800, 50);
      grok.shell.v = demogView;
      await w.__poll(() => grok.shell.tv?.dataFrame?.name, (n: string) => n === demogName, 800, 50);

      let switchError: string | null = null;
      await w.__tpRendered(() => {
        try { tp.props.table = dfCurves.name; } catch (e) { switchError = String(e); }
      }, 1500, 300);

      const boundToCurves = (() => {
        try { const d = tp.dataFrame; return {name: d?.name ?? null, rows: d?.rowCount ?? -1}; }
        catch { return {name: null, rows: -1}; }
      })();

      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const vs = root.querySelector('[name="viewer selector"]') as HTMLElement;
      vs.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
      await w.__poll(() => document.querySelector('.d4-combo-drop-down'), (e: Element | null) => !!e, 600, 40);
      const mc = document.querySelector('.d4-combo-drop-down [name="icon-multicurveviewer"]');
      const mcClicked = !!mc;
      await w.__tpRendered(() => (mc?.closest('.d4-list-item') as HTMLElement | null)?.click(), 1800, 300);
      const vt = tp.props.viewerType;

      let gearOpenedPropGrid = false;
      const panel = root.closest('.panel-base') as HTMLElement | null;
      const gear = panel?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
      if (gear) {
        gear.click();
        await w.__poll(() => document.querySelector('.property-grid'), (e: Element | null) => !!e, 1200, 60);
        gearOpenedPropGrid = !!document.querySelector('.property-grid');
      }

      let restoreError: string | null = null;
      await w.__tpRendered(() => {
        try { tp.props.table = demogName; } catch (e) { restoreError = String(e); }
      }, 1800, 300);
      const restoredCells = await w.__tpApply(() => {
        tp.props.viewerType = 'Scatter plot';
        tp.props.xColumnNames = ['SEX'];
        tp.props.yColumnNames = ['RACE'];
      }, w.__tpCells, 1800);
      const boundBack = (() => { try { return tp.dataFrame?.name ?? null; } catch { return null; } })();
      return {viewerType: vt, mcClicked, switchError, boundToCurves, restoreError, boundBack,
        curvesName: dfCurves.name, curvesRows: dfCurves.rowCount, demogName,
        restoredCells, gearOpenedPropGrid};
    }, curvesPath);
    expect(result.mcClicked).toBe(true);
    expect(['MultiCurveViewer', 'Multi curve viewer', 'Curves'].includes(result.viewerType)).toBe(true);
    expect(result.gearOpenedPropGrid).toBe(true);

    console.log(`[Multi Curve] table switch: error=${result.switchError} ` +
      `bound=${JSON.stringify(result.boundToCurves)} expected={"name":"${result.curvesName}","rows":${result.curvesRows}}`);
    console.log(`[Multi Curve] restore: error=${result.restoreError} boundBack=${result.boundBack} expected=${result.demogName}`);
    expect(result.switchError).toBeNull();
    expect(result.boundToCurves.name).toBe(result.curvesName);
    expect(result.boundToCurves.rows).toBe(result.curvesRows);
    expect(result.restoreError).toBeNull();
    expect(result.boundBack).toBe(result.demogName);
    expect(result.restoredCells).toBe(canonicalCellCount);
  });

  await softStep('To Script', async () => {
    await page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').first().click({button: 'right', position: {x: 6, y: 6}});
    await page.locator('.d4-menu-popup').last().waitFor({timeout: 10000});
    await page.evaluate(() => (window as any).__menuLeaf('To Script', 'To JavaScript'));
    const result = await page.evaluate(async () => {
      const balloon = await (window as any).__poll(() => document.querySelector('.d4-balloon'),
        (e: Element | null) => !!e, 1500, 60) as Element | null;
      const generated = !!balloon;
      if (balloon) {
        const close = balloon.querySelector('.close') || balloon.querySelector('[name="icon-times"]');
        if (close) (close as HTMLElement).click();
      }
      return {scriptGenerated: generated};
    });
    expect(result.scriptGenerated).toBe(true);
  });

  page.off('pageerror', onPageError);
  page.off('console', onConsole);
  await v.closeAllAndWait(page);
  v.finishSpec();
});
