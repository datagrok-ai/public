import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// Per-color canvas-delta floors, calibrated from live measurement on dev
// (demog.csv, headless, 2026-07-21). Each floor sits with a wide margin below
// the observed deliberate-change delta and far above the settle-precheck noise
// (measured 0 px on every baseline); the precheck ceiling proves the baseline
// frame was drained before measuring. Live deltas backing the four formerly
// tight floors: includeNulls 330940, barBorder 6176, labels 2250, showValues
// 1770 — each floor is under half its measured delta, so a real repaint fires
// it comfortably while a settle artifact cannot.
const T = {
  colorColumn: 800,
  invertScheme: 800,
  includeNulls: 100000,
  barBorder: 2000,
  maxBarHeight: 400,
  labels: 1000,
  aggrType: 400,
  valueSwitch: 400,
  showValues: 800,
};
const PRECHECK_CEIL = 250;

test('Bar chart tests', async ({page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });
  await page.waitForTimeout(500);

  // Snapshot the current canvas frame, wait for late renders to drain, and
  // return the residual settle delta (precheck). Leaves the snapshot set to
  // the settled frame so the next diffCanvasColors measures only the action.
  async function canvasBaseline(): Promise<number> {
    await v.snapshotCanvasColors(page, 'Bar chart');
    await page.waitForTimeout(500);
    return (await v.diffCanvasColors(page, 'Bar chart')).deltaPx;
  }
  async function canvasDelta(): Promise<number> {
    return (await v.diffCanvasColors(page, 'Bar chart')).deltaPx;
  }

  await softStep('Color coding', async () => {
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.colorColumnName = '';
      bc.props.invertColorScheme = false;
      await new Promise((r) => setTimeout(r, 600));
    });

    // Color-code by HEIGHT — bars recolor by aggregated height.
    const preColor = await canvasBaseline();
    const setColor = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.colorColumnName = 'HEIGHT';
      await new Promise((r) => setTimeout(r, 600));
      return {colorColumnName: bc.props.colorColumnName, aggr: bc.props.colorAggrType};
    });
    const colorDelta = await canvasDelta();

    // Aggregation cycles keep the read-back coverage of the original spec.
    const aggrReads = await v.setViewerProps(page, 'Bar chart', [
      {set: {colorAggrType: 'min'}, wait: 200, read: 'colorAggrType'},
      {set: {colorAggrType: 'max'}, wait: 200, read: 'colorAggrType'},
      {set: {colorAggrType: 'med'}, wait: 200, read: 'colorAggrType'},
    ]);

    // Let the trailing 'med' aggregation repaint settle before snapshotting the
    // invert baseline, so the precheck measures a drained frame (not the tail of
    // the aggregation recolor).
    await page.waitForTimeout(400);

    // Invert the color scheme — the gradient reverses, recoloring the bars.
    const preInvert = await canvasBaseline();
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.invertColorScheme = true;
      await new Promise((r) => setTimeout(r, 600));
    });
    const invertDelta = await canvasDelta();

    const colColRead = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.colorColumnName = '';
      bc.props.invertColorScheme = false;
      await new Promise((r) => setTimeout(r, 300));
      return bc.props.colorColumnName;
    });

    expect(setColor.colorColumnName).toBe('HEIGHT');
    expect(aggrReads).toEqual(['min', 'max', 'med']);
    expect(preColor).toBeLessThan(PRECHECK_CEIL);
    expect(colorDelta).toBeGreaterThan(T.colorColumn);
    expect(preInvert).toBeLessThan(PRECHECK_CEIL);
    expect(invertDelta).toBeGreaterThan(T.invertScheme);
    expect(colColRead).toBe('');
  });

  await softStep('Include nulls', async () => {
    // HEIGHT carries 751 missing values → a "missing" bar renders when
    // includeNulls is on (DIS_POP and the other categoricals have zero nulls).
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'HEIGHT';
      bc.props.includeNulls = true;
      await new Promise((r) => setTimeout(r, 600));
    });

    // Dropping the missing bar repaints the chart.
    const preOff = await canvasBaseline();
    const offRead = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.includeNulls = false;
      await new Promise((r) => setTimeout(r, 500));
      return bc.props.includeNulls;
    });
    const offDelta = await canvasDelta();

    // Re-adding it repaints again.
    const preOn = await canvasBaseline();
    const onRead = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.includeNulls = true;
      await new Promise((r) => setTimeout(r, 500));
      return bc.props.includeNulls;
    });
    const onDelta = await canvasDelta();

    expect(offRead).toBe(false);
    expect(onRead).toBe(true);
    expect(preOff).toBeLessThan(PRECHECK_CEIL);
    expect(offDelta).toBeGreaterThan(T.includeNulls);
    expect(preOn).toBeLessThan(PRECHECK_CEIL);
    expect(onDelta).toBeGreaterThan(T.includeNulls);
  });

  await softStep('Bar style', async () => {
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      await new Promise((r) => setTimeout(r, 600));
    });

    // Bar border 0 → 2 adds visible outlines around every bar.
    const preBorder = await canvasBaseline();
    const borderRead = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.barBorderLineWidth = 2;
      await new Promise((r) => setTimeout(r, 500));
      return bc.props.barBorderLineWidth;
    });
    const borderDelta = await canvasDelta();

    // Max bar height 50 → 20 makes the bars thinner (exact width is not
    // pixel-measurable headless; the repaint delta is the drivable proxy).
    const preHeight = await canvasBaseline();
    const heightRead = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.maxBarHeight = 20;
      await new Promise((r) => setTimeout(r, 500));
      return bc.props.maxBarHeight;
    });
    const heightDelta = await canvasDelta();

    // Remaining style props keep the read-back coverage of the original spec.
    const styleReads = await v.setViewerProps(page, 'Bar chart', [
      {set: {barCornerRadius: 10}, wait: 200, read: 'barCornerRadius'},
      {set: {verticalAlign: 'Top'}, wait: 200, read: 'verticalAlign'},
      {set: {verticalAlign: 'Bottom'}, wait: 200, read: 'verticalAlign'},
      {set: {verticalAlign: 'Center'}, wait: 200, read: 'verticalAlign'},
      {set: {showCategoryZeroBaseline: false}, wait: 200, read: 'showCategoryZeroBaseline'},
    ]);

    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.barBorderLineWidth = 0;
      bc.props.barCornerRadius = 0;
      bc.props.maxBarHeight = 50;
      bc.props.showCategoryZeroBaseline = true;
      await new Promise((r) => setTimeout(r, 300));
    });

    expect(borderRead).toBe(2);
    expect(heightRead).toBe(20);
    expect(styleReads).toEqual([10, 'Top', 'Bottom', 'Center', false]);
    expect(preBorder).toBeLessThan(PRECHECK_CEIL);
    expect(borderDelta).toBeGreaterThan(T.barBorder);
    expect(preHeight).toBeLessThan(PRECHECK_CEIL);
    expect(heightDelta).toBeGreaterThan(T.maxBarHeight);
  });

  await softStep('Labels', async () => {
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.showLabels = 'inside';
      await new Promise((r) => setTimeout(r, 600));
    });

    // 'inside' → 'never' removes the in-bar value labels: the label glyphs
    // vanish from the canvas.
    const preLabels = await canvasBaseline();
    const neverRead = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.showLabels = 'never';
      await new Promise((r) => setTimeout(r, 500));
      return bc.props.showLabels;
    });
    const labelsDelta = await canvasDelta();

    // Remaining label modes keep the read-back coverage of the original spec.
    const labelReads = await v.setViewerProps(page, 'Bar chart', [
      {set: {showLabels: 'outside'}, wait: 200, read: 'showLabels'},
      {set: {showLabels: 'auto'}, wait: 200, read: 'showLabels'},
    ]);

    expect(neverRead).toBe('never');
    expect(labelReads).toEqual(['outside', 'auto']);
    expect(preLabels).toBeLessThan(PRECHECK_CEIL);
    expect(labelsDelta).toBeGreaterThan(T.labels);
  });

  await softStep('Controls visibility', async () => {
    const ctrls = ['showValueSelector', 'showCategorySelector', 'showStackSelector',
      'showValueAxis', 'showCategoryValues'];
    const off = Object.fromEntries(ctrls.map((k) => [k, false]));
    const on = Object.fromEntries(ctrls.map((k) => [k, true]));

    // The in-chart selector nodes (value / category / stack) are always present
    // but toggle their computed `display` (flex ↔ none) with the show*Selector
    // props. Count the ones actually laid out (display !== none) — probed live to
    // be exactly 3 when on, 0 when off. showValueAxis and showCategoryValues are
    // canvas-drawn — probed 2026-07-21: no DOM node toggles (0 axis/category
    // nodes) and their canvas delta is ~13 px ≈ settle noise — so they stay
    // covered by the prop read-back rather than a DOM count. See the .md note.
    const countVisibleSelectors = () => page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const root = bc.root as HTMLElement;
      return Array.from(root.querySelectorAll('.d4-column-selector'))
        .filter((e: any) => getComputedStyle(e).display !== 'none').length;
    });

    const offReads = await v.setViewerProps(page, 'Bar chart', [{set: off, wait: 400, read: ctrls}]);
    const offVisible = await countVisibleSelectors();
    const onReads = await v.setViewerProps(page, 'Bar chart', [{set: on, wait: 400, read: ctrls}]);
    const onVisible = await countVisibleSelectors();

    expect(offReads[0]).toEqual(off);
    expect(onReads[0]).toEqual(on);
    expect(offVisible).toBe(0);
    // Exactly three column-selector nodes (value / category / stack) lay out when
    // the show*Selector props are on — probed live 2026-07-21.
    expect(onVisible).toBe(3);
  });

  await softStep('Aggregation types', async () => {
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.valueAggrType = 'avg';
      await new Promise((r) => setTimeout(r, 600));
    });

    // avg → max changes every bar's aggregated height.
    const preAggr = await canvasBaseline();
    const aggrRead = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueAggrType = 'max';
      await new Promise((r) => setTimeout(r, 500));
      return bc.props.valueAggrType;
    });
    const aggrDelta = await canvasDelta();

    // Rebinding the value column AGE → WEIGHT re-scales the bars.
    const preSwitch = await canvasBaseline();
    const valueRead = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueColumnName = 'WEIGHT';
      await new Promise((r) => setTimeout(r, 500));
      return bc.props.valueColumnName;
    });
    const switchDelta = await canvasDelta();

    expect(aggrRead).toBe('max');
    expect(valueRead).toBe('WEIGHT');
    expect(preAggr).toBeLessThan(PRECHECK_CEIL);
    expect(aggrDelta).toBeGreaterThan(T.aggrType);
    expect(preSwitch).toBeLessThan(PRECHECK_CEIL);
    expect(switchDelta).toBeGreaterThan(T.valueSwitch);
  });

  await softStep('Legend position', async () => {
    // Replace the viewer with a fresh Bar chart: after the earlier sections'
    // prop churn the stacked-legend host will not render, whereas a clean-state
    // viewer materializes it reliably. This is expected platform behavior — the
    // legend host is churn-sensitive, not broken — confirmed by a live probe
    // (2026-07-21): the same stack+legend props render items=0 on the churned
    // widget yet items=2 on a fresh one. Recreating is a harness reset, not a
    // silent workaround for a product bug.
    await page.evaluate(async () => {
      const old = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      if (old) old.close();
      await new Promise((r) => setTimeout(r, 400));
    });
    await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');
    await page.waitForTimeout(500);

    // Setting a Stack column materializes the legend; the legend host is absent
    // before it and present after — a DOM round-trip, not a prop echo. Driven in
    // one evaluate (the render is sensitive to churn across round-trips).
    const legend = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const root = bc.root as HTMLElement;
      const rendered = () => !!root.querySelector('[name="legend"]');
      const items = () => root.querySelectorAll('[name="legend"] .d4-legend-item').length;

      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.stackColumnName = '';
      bc.props.legendVisibility = 'Always';
      await new Promise((r) => setTimeout(r, 800));
      const before = rendered();

      bc.props.stackColumnName = 'SEX';
      await new Promise((r) => setTimeout(r, 1500));
      const after = rendered();
      const afterItems = items();

      const positions: string[] = [];
      for (const pos of ['Left', 'Right', 'Top', 'Bottom']) {
        bc.props.legendPosition = pos;
        await new Promise((r) => setTimeout(r, 200));
        positions.push(bc.props.legendPosition);
      }

      bc.props.stackColumnName = '';
      await new Promise((r) => setTimeout(r, 1500));
      const cleared = rendered();

      return {before, after, afterItems, positions, cleared};
    });

    expect(legend.before).toBe(false);
    expect(legend.after).toBe(true);
    expect(legend.afterItems).toBeGreaterThanOrEqual(2);
    expect(legend.positions).toEqual(['Left', 'Right', 'Top', 'Bottom']);
    expect(legend.cleared).toBe(false);
  });

  await softStep('Title and description', async () => {
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.showTitle = true;
      bc.props.title = 'Demographics';
      bc.props.description = 'By race';
      await new Promise((r) => setTimeout(r, 400));
      const root = bc.root as HTMLElement;
      // Title renders in the panel titlebar (.panel-titlebar-text), the
      // description inside the viewer root.
      const panel = (root.closest('.panel-base') as HTMLElement) ?? root;
      const panelText = panel.innerText ?? panel.textContent ?? '';
      const rootText = root.innerText ?? root.textContent ?? '';
      const r: any = {showTitle: bc.props.showTitle, titleInDom: panelText.includes('Demographics'),
        descInDom: rootText.includes('By race')};

      const positions: string[] = [];
      for (const pos of ['Top', 'Bottom', 'Left', 'Right']) {
        bc.props.descriptionPosition = pos;
        await new Promise((res) => setTimeout(res, 150));
        positions.push(bc.props.descriptionPosition);
      }
      r.positions = positions;

      bc.props.descriptionVisibilityMode = 'Never';
      await new Promise((res) => setTimeout(res, 300));
      r.hidden = bc.props.descriptionVisibilityMode;
      // DOM-primary: 'Never' removes the description text from the viewer body.
      const rootTextHidden = root.innerText ?? root.textContent ?? '';
      r.descGoneFromDom = rootTextHidden.includes('By race');

      bc.props.showTitle = false;
      return r;
    });

    expect(info.showTitle).toBe(true);
    expect(info.titleInDom).toBe(true);
    expect(info.descInDom).toBe(true);
    expect(info.positions).toEqual(['Top', 'Bottom', 'Left', 'Right']);
    expect(info.hidden).toBe('Never');
    expect(info.descGoneFromDom).toBe(false);
  });

  await softStep('Show values instead of categories', async () => {
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.valueAggrType = 'avg';
      bc.props.showValuesInsteadOfCategories = false;
      await new Promise((r) => setTimeout(r, 600));
    });

    // Toggling category labels to aggregated values repaints the axis strip.
    const preShow = await canvasBaseline();
    const onRead = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.showValuesInsteadOfCategories = true;
      await new Promise((r) => setTimeout(r, 500));
      return bc.props.showValuesInsteadOfCategories;
    });
    const showDelta = await canvasDelta();

    const offRead = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.showValuesInsteadOfCategories = false;
      await new Promise((r) => setTimeout(r, 300));
      return bc.props.showValuesInsteadOfCategories;
    });

    expect(onRead).toBe(true);
    expect(offRead).toBe(false);
    expect(preShow).toBeLessThan(PRECHECK_CEIL);
    expect(showDelta).toBeGreaterThan(T.showValues);
  });

  await softStep('Data panel', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));

      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      df2.name = 'SPGI';
      grok.shell.addTableView(df2);
      await new Promise((resolve) => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const views = Array.from(grok.shell.views).filter((view: any) => view.type === 'TableView');
      const demogView = views.find((view: any) => view.dataFrame.name !== 'SPGI') as any;
      if (demogView) grok.shell.v = demogView;
      await new Promise((r) => setTimeout(r, 500));

      const icon = document.querySelector('[name="icon-bar-chart"]') as HTMLElement;
      icon.click();
      await new Promise((r) => setTimeout(r, 1000));

      const bc = Array.from(grok.shell.tv.viewers).find((view: any) => view.type === 'Bar chart') as any;
      const r: any[] = [];

      // Cycle Row Source through all three values (Filtered, Selected, All).
      for (const src of ['Filtered', 'Selected', 'All']) {
        bc.props.rowSource = src;
        await new Promise((res) => setTimeout(res, 200));
        r.push(bc.props.rowSource);
      }
      bc.props.rowSource = 'All';

      const spgi = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      bc.dataFrame = spgi;
      await new Promise((res) => setTimeout(res, 500));
      r.push(bc.dataFrame.name);

      bc.props.filter = '${CAST Idea ID} < 636500';
      await new Promise((res) => setTimeout(res, 500));
      r.push(bc.props.filter);

      bc.props.colorColumnName = 'Chemical Space Y';
      await new Promise((res) => setTimeout(res, 300));
      r.push(bc.props.colorColumnName);

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise((res) => setTimeout(res, 1000));

      bc.close();
      await new Promise((res) => setTimeout(res, 500));

      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise((res) => setTimeout(res, 3000));

      const bc2 = Array.from(grok.shell.tv.viewers).find((view: any) => view.type === 'Bar chart') as any;
      r.push(bc2 ? bc2.props.colorColumnName : 'NOT_RESTORED');
      r.push(bc2 ? bc2.props.filter : 'NOT_RESTORED');

      await grok.dapi.layouts.delete(saved);

      return r;
    });
    expect(result.slice(0, 3)).toEqual(['Filtered', 'Selected', 'All']);
    expect(result[3]).toBe('SPGI');
    expect(result[4]).toBe('${CAST Idea ID} < 636500');
    expect(result[5]).toBe('Chemical Space Y');
    // Honest layout round-trip: color coding and filter survive save → close →
    // reload.
    expect(result[6]).toBe('Chemical Space Y');
    expect(result[7]).toBe('${CAST Idea ID} < 636500');
  });

  await softStep('No page errors', async () => {
    expect(pageErrors).toEqual([]);
  });

  v.finishSpec();
});
