/* ---
realizes: [pcplot.cp.transformation-and-filter-integrity]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';


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
  // GROK-18091 surfaces as a FormatException in the PC-Plot transformation path
  // (GroupAggregation.fromJsonList → PcPlotCore._transform). The shared dev server
  // emits unrelated console noise, so match only errors naming that surface.
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

  // Open the Filter Panel — its header carries the Reset filters button
  // [name="icon-arrow-rotate-left"] used across the scenarios below.
  await page.evaluate(() => {
    grok.shell.tv.getFiltersGroup();
  });
  await page.locator('.d4-filter-group-header').waitFor({timeout: 15000});

  const fullCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);

  await softStep('Scenario 1 (GROK-17306) — Reset filters restores filter, keeps selection', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const pc = tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      // GROK-17306 reproduces with a transformation present. Set it via the
      // property (Context Panel > Data > Transformation), no dialog. The pivot
      // replaces the PC-plot axes, so the filter is applied on the FILTER PANEL
      // (its widgets persist), not on an in-chart axis slider.
      pc.props.transformation = JSON.stringify([
        {'#type': 'GroupAggregation', aggType: 'key', colName: 'SEX'},
        {'#type': 'GroupAggregation', aggType: 'pivot', colName: 'DIS_POP'},
        {'#type': 'GroupAggregation', aggType: 'avg', colName: 'WEIGHT'},
      ]);
      await new Promise((r) => setTimeout(r, 900));
      // Apply a range filter through the Filter Panel.
      tv.getFiltersGroup().updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 50});
      await new Promise((r) => setTimeout(r, 600));
      const filteredCount = df.filter.trueCount;
      df.selection.init((i: number) => i < 10);
      await new Promise((r) => setTimeout(r, 200));
      const selCount = df.selection.trueCount;
      // Click the Filter Panel "Reset filters" button.
      const btn = document.querySelector(
        '.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement;
      btn.click();
      await new Promise((r) => setTimeout(r, 700));
      const out = {
        filteredCount, selCount,
        restoredFilter: df.filter.trueCount,
        selAfterReset: df.selection.trueCount,
      };
      pc.props.transformation = '';
      await new Promise((r) => setTimeout(r, 400));
      return out;
    });
    // The range filter is active (below full count) with the transformation present.
    expect(result.filteredCount).toBeLessThan(fullCount);
    expect(result.selCount).toBeGreaterThan(0);
    // Reset restores the filter to the full row count …
    expect(result.restoredFilter).toBe(fullCount);
    // … and GROK-17306: resetting the filters must not clear the selection.
    expect(result.selAfterReset).toBe(result.selCount);
  });

  await softStep('Scenario 2 (GROK-18489) — second filter after DateTime color split works', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const pc = tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      (document.querySelector(
        '.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 600));
      // STARTED is the DateTime column.
      pc.props.colorColumnName = 'STARTED';
      await new Promise((r) => setTimeout(r, 500));
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const rect = viewer.getBoundingClientRect();
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
      const revealAndSlider = async () => {
        viewer.dispatchEvent(new MouseEvent('mousemove', {
          bubbles: true, clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
        }));
        await new Promise((r) => setTimeout(r, 400));
        return document.querySelector('[name="axis-slider-AGE"]')!;
      };
      const dragHandle = async (handleName: string, dir: number, dist: number) => {
        const svg = await revealAndSlider();
        const h = svg.querySelector(`[name="${handleName}"]`)!;
        const hr = h.getBoundingClientRect();
        const cx = hr.x + hr.width / 2;
        const cy = hr.y + hr.height / 2;
        h.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
        await new Promise((r) => setTimeout(r, 50));
        for (let d = 20; d <= dist; d += 30) {
          document.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dir * d)));
          svg.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dir * d)));
          await new Promise((r) => setTimeout(r, 20));
        }
        document.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy + dir * dist)));
        await new Promise((r) => setTimeout(r, 600));
      };
      await dragHandle('max-handle', 1, 300);
      const firstFiltered = df.filter.trueCount;
      (document.querySelector(
        '.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 700));
      const afterReset = df.filter.trueCount;
      // A different value window: drag the min-handle up.
      await dragHandle('min-handle', -1, 250);
      const secondFiltered = df.filter.trueCount;
      // One more adjustment — the count must move again, not freeze.
      await dragHandle('max-handle', 1, 200);
      const thirdFiltered = df.filter.trueCount;
      return {colorCol: pc.props.colorColumnName, firstFiltered, afterReset, secondFiltered, thirdFiltered};
    });
    expect(result.colorCol).toBe('STARTED');
    expect(result.firstFiltered).toBeLessThan(fullCount);
    expect(result.afterReset).toBe(fullCount);
    // GROK-18489: a second filter after a reset takes effect, and the count keeps
    // responding to further adjustments.
    expect(result.secondFiltered).toBeLessThan(fullCount);
    expect(result.thirdFiltered).not.toBe(result.secondFiltered);
  });

  await softStep('Scenario 3 (github-972) — histogram column change does not reset PC filter', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      (document.querySelector(
        '.d4-filter-group-header [name="icon-arrow-rotate-left"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 600));
      const hist = tv.addViewer('Histogram');
      await new Promise((r) => setTimeout(r, 1000));
      // Apply a PC-plot AGE range filter.
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const rect = viewer.getBoundingClientRect();
      viewer.dispatchEvent(new MouseEvent('mousemove', {
        bubbles: true, clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      await new Promise((r) => setTimeout(r, 400));
      const svg = document.querySelector('[name="axis-slider-AGE"]')!;
      const mh = svg.querySelector('[name="max-handle"]')!;
      const hr = mh.getBoundingClientRect();
      const cx = hr.x + hr.width / 2;
      const cy = hr.y + hr.height / 2;
      const mk = (x: number, y: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
      mh.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
      await new Promise((r) => setTimeout(r, 50));
      for (let dy = 20; dy <= 300; dy += 30) {
        document.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
        svg.dispatchEvent(new MouseEvent('mousemove', mk(cx, cy + dy)));
        await new Promise((r) => setTimeout(r, 20));
      }
      document.dispatchEvent(new MouseEvent('mouseup', mk(cx, cy + 300)));
      await new Promise((r) => setTimeout(r, 600));
      const filteredCount = df.filter.trueCount;
      hist.props.valueColumnName = 'HEIGHT';
      await new Promise((r) => setTimeout(r, 500));
      return {
        filteredCount,
        histColAfter: hist.props.valueColumnName,
        filterAfter: df.filter.trueCount,
      };
    });
    expect(result.filteredCount).toBeLessThan(fullCount);
    expect(result.histColAfter).toBe('HEIGHT');
    // github-972: changing the histogram column must not reset the PC-plot filter.
    expect(result.filterAfter).toBe(result.filteredCount);
  });

  await softStep('Scenario 4 (GROK-18091) — aggregation + Filter Panel close, no broken state', async () => {
    // Baseline count of transformation / PC-Plot / aggregation errors only.
    const bugErrsBefore =
      pageErrors.filter(isBugRelevantError).length + consoleErrors.filter(isBugRelevantError).length;
    const state = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const pc = tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      // The transformation property is a JSON list of GroupAggregations (parsed by
      // GroupAggregation.fromJsonList); a non-JSON string throws a FormatException.
      pc.props.transformation = JSON.stringify([{aggType: 'avg', colName: 'AGE'}]);
      await new Promise((r) => setTimeout(r, 800));
      // Close the Filter Panel via its close (X) icon.
      const fp = document.querySelector('[name="viewer-Filters"]')!;
      const closeIcon = fp.querySelector('[name="icon-times"]') as HTMLElement;
      if (closeIcon)
        closeIcon.click();
      await new Promise((r) => setTimeout(r, 500));
      // GROK-18091 shows as a broken visual state: the DataFrame goes blank or the
      // PC-Plot viewer disappears.
      const pcStillPresent =
        !!tv.viewers.find((vw: any) => vw.type === 'PC Plot') &&
        !!document.querySelector('[name="viewer-PC-Plot"]');
      return {rowCount: df.rowCount, pcStillPresent, transformationApplied: pc.props.transformation};
    });
    // grok.shell.warnings is undefined on this build, so pageerror/console is the
    // error channel.
    const bugErrsAfter =
      pageErrors.filter(isBugRelevantError).length + consoleErrors.filter(isBugRelevantError).length;
    expect(bugErrsAfter).toBe(bugErrsBefore);
    expect(state.rowCount).toBeGreaterThan(0);
    expect(state.pcStillPresent).toBe(true);
  });

  v.finishSpec();
});
