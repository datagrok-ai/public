/* ---
realizes: [pcplot.cp.transformation-and-filter-integrity]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {knownOpenBug} from '../../helpers/known-open-bug';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('PC Plot — Transformation and Filter/Selection Integrity', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error')
      consoleErrors.push(m.text());
  });

  const isBugRelevantError = (s: string): boolean =>
    /transform|aggregat|GroupAggregation|FormatException|pc[_-]?plot|PcPlot|is not valid JSON/i.test(s);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon)
      (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 15000});

  await v.installEventWaits(page);

  await page.evaluate(() => {
    grok.shell.tv.getFiltersGroup();
  });
  await page.locator('.d4-filter-group-header').waitFor({timeout: 15000});

  await page.evaluate(() => {
    (window as any).__settle = async (stream: any, capMs: number, act?: () => void) => {
      const done = new Promise<void>((resolve) => {
        let sub: any = null;
        try { sub = stream.subscribe(() => { sub.unsubscribe(); resolve(); }); }
        catch (_) {  }
        setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, capMs);
      });
      if (act) act();
      await done;
    };
  });

  const fullCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);

  await softStep('Scenario 1 (GROK-17306) — Reset filters restores filter, keeps selection', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const pc = tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const settle = (window as any).__settle;

      await settle(pc.onViewerRendered, 900, () => {
        pc.props.transformation = JSON.stringify([
          {'#type': 'GroupAggregation', aggType: 'key', colName: 'SEX'},
          {'#type': 'GroupAggregation', aggType: 'pivot', colName: 'DIS_POP'},
          {'#type': 'GroupAggregation', aggType: 'avg', colName: 'WEIGHT'},
        ]);
      });

      await settle(df.onRowsFiltered, 600, () =>
        tv.getFiltersGroup().updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 50}));
      const filteredCount = df.filter.trueCount;
      await settle(df.onSelectionChanged, 200, () => df.selection.init((i: number) => i < 10));
      const selCount = df.selection.trueCount;

      const btn = document.querySelector(
        '.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement;
      await settle(df.onRowsFiltered, 700, () => btn.click());
      const out = {
        filteredCount, selCount,
        restoredFilter: df.filter.trueCount,
        selAfterReset: df.selection.trueCount,
      };
      await settle(pc.onViewerRendered, 400, () => { pc.props.transformation = ''; });
      return out;
    });

    expect(result.filteredCount).toBeLessThan(fullCount);
    expect(result.selCount).toBeGreaterThan(0);

    expect(result.restoredFilter).toBe(fullCount);

    await knownOpenBug('GROK-17306', () => {
      expect(result.selAfterReset).toBe(result.selCount);
    });
  });

  await softStep('Scenario 2 (GROK-18489) — second filter after DateTime color split works', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const pc = tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const settle = (window as any).__settle;
      await settle(df.onRowsFiltered, 600, () => (document.querySelector(
        '.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement).click());

      await settle(pc.onViewerRendered, 500, () => { pc.props.colorColumnName = 'STARTED'; });
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const rect = viewer.getBoundingClientRect();
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
      const revealAndSlider = async () => {
        viewer.dispatchEvent(new MouseEvent('mousemove', {
          bubbles: true, clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
        }));

        return w.__poll(() => document.querySelector('[name="axis-slider-AGE"]'),
          (el: Element | null) => el !== null, 400);
      };
      const dragHandle = async (handleName: string, dir: number, dist: number) => {
        const svg = await revealAndSlider();
        const h = svg.querySelector(`[name="${handleName}"]`)!;
        const hr = h.getBoundingClientRect();
        const cx = hr.x + hr.width / 2;
        const cy = hr.y + hr.height / 2;

        h.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
        await w.__drag(svg as HTMLElement, {x: cx, y: cy + dir * 20}, {x: cx, y: cy + dir * dist},
          {steps: Math.max(1, Math.ceil((dist - 20) / 30)), stepMs: 20, holdMs: 50});
        await settle(df.onRowsFiltered, 600, () =>
          document.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy + dir * dist))));
      };
      await dragHandle('max-handle', 1, 300);
      const firstFiltered = df.filter.trueCount;
      await settle(df.onRowsFiltered, 700, () => (document.querySelector(
        '.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement).click());
      const afterReset = df.filter.trueCount;

      await dragHandle('min-handle', -1, 250);
      const secondFiltered = df.filter.trueCount;

      await dragHandle('max-handle', 1, 200);
      const thirdFiltered = df.filter.trueCount;
      return {colorCol: pc.props.colorColumnName, firstFiltered, afterReset, secondFiltered, thirdFiltered};
    });
    expect(result.colorCol).toBe('STARTED');
    expect(result.firstFiltered).toBeLessThan(fullCount);
    expect(result.afterReset).toBe(fullCount);

    expect(result.secondFiltered).toBeLessThan(fullCount);
    expect(result.thirdFiltered).not.toBe(result.secondFiltered);
  });

  await softStep('Scenario 3 (github-972) — histogram column change does not reset PC filter', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const settle = (window as any).__settle;
      await settle(df.onRowsFiltered, 600, () => (document.querySelector(
        '.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement).click());
      const hist = tv.addViewer('Histogram');
      await settle(hist.onViewerRendered, 1000);

      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const rect = viewer.getBoundingClientRect();
      viewer.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));

      const svg: Element = await w.__poll(() => document.querySelector('[name="axis-slider-AGE"]'),
        (el: Element | null) => el !== null, 400);
      const mh = svg.querySelector('[name="max-handle"]')!;
      const hr = mh.getBoundingClientRect();
      const cx = hr.x + hr.width / 2;
      const cy = hr.y + hr.height / 2;
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});

      mh.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
      await w.__drag(svg as HTMLElement, {x: cx, y: cy + 20}, {x: cx, y: cy + 300},
        {steps: 10, stepMs: 20, holdMs: 50});
      await settle(df.onRowsFiltered, 600, () =>
        document.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy + 300))));
      const filteredCount = df.filter.trueCount;
      await w.__settled('viewer:Histogram.onViewerRendered',
        () => { hist.props.valueColumnName = 'HEIGHT'; }, 500);
      return {
        filteredCount,
        histColAfter: hist.props.valueColumnName,
        filterAfter: df.filter.trueCount,
      };
    });
    expect(result.filteredCount).toBeLessThan(fullCount);
    expect(result.histColAfter).toBe('HEIGHT');

    expect(result.filterAfter).toBe(result.filteredCount);
  });

  await softStep('Scenario 4 (GROK-18091) — aggregation + Filter Panel close, no broken state', async () => {

    const bugErrsBefore =
      pageErrors.filter(isBugRelevantError).length + consoleErrors.filter(isBugRelevantError).length;
    const state = await page.evaluate(async () => {
      const w = window as any;
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const pc = tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      const settle = (window as any).__settle;

      await settle(pc.onViewerRendered, 800, () => {
        pc.props.transformation = JSON.stringify([{'#type': 'GroupAggregation', aggType: 'avg', colName: 'AGE'}]);
      });

      const fp = document.querySelector('[name="viewer-Filters"]')!;
      const closeIcon = fp.querySelector('[name="icon-times"]') as HTMLElement;
      if (closeIcon)
        await w.__settled('grok.events.onViewerClosed', () => closeIcon.click(), 500);

      const pcStillPresent =
        !!tv.viewers.find((vw: any) => vw.type === 'PC Plot') &&
        !!document.querySelector('[name="viewer-PC-Plot"]');
      return {rowCount: df.rowCount, pcStillPresent, transformationApplied: pc.props.transformation};
    });

    const bugErrsAfter =
      pageErrors.filter(isBugRelevantError).length + consoleErrors.filter(isBugRelevantError).length;
    expect(bugErrsAfter).toBe(bugErrsBefore);
    expect(state.rowCount).toBeGreaterThan(0);
    expect(state.pcStillPresent).toBe(true);
  });

  v.finishSpec();
});
