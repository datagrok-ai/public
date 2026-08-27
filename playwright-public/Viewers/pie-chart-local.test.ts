import {test, expect} from '@datagrok-libraries/test/src/playwright/local-fixtures';
import {specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

// Local-mode twin of pie-chart.test.ts: same scenarios, no server and no fixed sleeps.
// `?mode=local` removes the token exchange, the boot round-trips and the DemoFiles read
// (core/docs/features/ui2/LOCAL_MODE.md); the client boots once per worker (local-fixtures)
// instead of once per file, and every property step waits for the pie's own
// d4-viewer-rendered instead of 200 ms. The two server-backed scenarios are kept, not
// dropped: layout persistence round-trips through ViewLayout JSON in memory, and table
// switching uses a second generated table.
test.use(specTestOptions);

test('Pie chart tests (local mode)', async ({localPage: page}) => {
  test.setTimeout(120_000);

  await page.evaluate(async () => {
    const df = grok.data.demo.demog(5000);
    df.name = 'demog';
    // The generator emits no missing values; the Include-nulls scenario needs some.
    const disease = df.col('disease')!;
    for (let i = 0; i < df.rowCount; i += 7) disease.set(i, null, false);
    grok.shell.addTableView(df);
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await page.evaluate(async () => {
    const pie = grok.shell.tv.addViewer('Pie chart');
    await new Promise((res) => {
      const sub = pie.onViewerRendered.subscribe(() => { sub.unsubscribe(); res(null); });
      setTimeout(res, 5000);
    });
  });
  await page.locator('[name="viewer-Pie-chart"]').waitFor({timeout: 30_000});
  await v.installRenderWait(page);

  // #### Sorting
  await softStep('Sorting', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {categoryColumnName: 'race', pieSortType: 'by value'}, read: 'pieSortType'},
      {set: {pieSortOrder: 'desc'}, read: 'pieSortOrder'},
      {set: {pieSortOrder: 'asc'}, read: 'pieSortOrder'},
      {set: {pieSortType: 'by category'}, read: 'pieSortType'},
      {set: {pieSortOrder: 'asc'}, read: 'pieSortOrder'},
      {set: {pieSortOrder: 'desc'}, read: 'pieSortOrder'},
    ], 'render');
    expect(result).toEqual(['by value', 'desc', 'asc', 'by category', 'asc', 'desc']);
  });

  // #### Segment angle and length
  await softStep('Segment angle and length', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {categoryColumnName: 'race', segmentAngleColumnName: 'age'}, read: 'segmentAngleColumnName'},
      {set: {segmentAngleAggrType: 'sum'}, read: 'segmentAngleAggrType'},
      {set: {segmentAngleAggrType: 'count'}, read: 'segmentAngleAggrType'},
      {set: {segmentLengthColumnName: 'weight'}, read: 'segmentLengthColumnName'},
      {set: {segmentLengthAggrType: 'max'}, read: 'segmentLengthAggrType'},
      {set: {segmentAngleColumnName: '', segmentLengthColumnName: ''}},
    ], 'render');
    expect(result[0]).toBe('age');
    expect(result[1]).toBe('sum');
    expect(result[2]).toBe('count');
    expect(result[3]).toBe('weight');
    expect(result[4]).toBe('max');
  });

  // #### Appearance
  await softStep('Appearance', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {startAngle: 90}, read: 'startAngle'},
      {set: {startAngle: 180}, read: 'startAngle'},
      {set: {startAngle: 0}, read: 'startAngle'},
      {set: {maxRadius: 100}, read: 'maxRadius'},
      {set: {maxRadius: 150}, read: 'maxRadius'},
      {set: {shift: 10}, read: 'shift'},
      {set: {shift: 0}, read: 'shift'},
    ], 'render');
    expect(result).toEqual([90, 180, 0, 100, 150, 10, 0]);
  });

  // #### Labels
  await softStep('Labels', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {labelPosition: 'Inside'}, read: 'labelPosition'},
      {set: {labelPosition: 'Outside'}, read: 'labelPosition'},
      {set: {labelPosition: 'Auto'}, read: 'labelPosition'},
      {set: {showLabel: false}, read: 'showLabel'},
      {set: {showPercentage: false}, read: 'showPercentage'},
      {set: {showValue: true}, read: 'showValue'},
      {set: {showLabel: true, showPercentage: true}},
    ], 'render');
    expect(result).toEqual(['Inside', 'Outside', 'Auto', false, false, true]);
  });

  // #### Outline
  await softStep('Outline', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {outlineLineWidth: 5}, read: 'outlineLineWidth'},
      {set: {outlineLineWidth: 0}, read: 'outlineLineWidth'},
      {set: {outlineLineWidth: 1}, read: 'outlineLineWidth'},
    ], 'render');
    expect(result).toEqual([5, 0, 1]);
  });

  // #### Include nulls
  await softStep('Include nulls', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {categoryColumnName: 'disease', includeNulls: true}, read: 'includeNulls'},
      {set: {includeNulls: false}, read: 'includeNulls'},
      {set: {includeNulls: true, categoryColumnName: 'race'}},
    ], 'render');
    expect(result).toEqual([true, false]);
  });

  // #### Column selector
  await softStep('Column selector', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {showColumnSelector: false}, read: 'showColumnSelector'},
      {set: {showColumnSelector: true}, read: 'showColumnSelector'},
    ], 'render');
    expect(result).toEqual([false, true]);
  });

  // #### Legend
  await softStep('Legend', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {legendVisibility: 'Always'}, read: 'legendVisibility'},
      {set: {legendPosition: 'LeftTop'}, read: 'legendPosition'},
      {set: {legendPosition: 'RightBottom'}, read: 'legendPosition'},
      {set: {legendVisibility: 'Never'}, read: 'legendVisibility'},
      {set: {legendVisibility: 'Auto'}, read: 'legendVisibility'},
    ], 'render');
    expect(result).toEqual(['Always', 'LeftTop', 'RightBottom', 'Never', 'Auto']);
  });

  // #### Category map (dates)
  await softStep('Category map (dates)', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {categoryColumnName: 'started'}, read: 'categoryMap'},
      {set: {categoryMap: 'month'}, read: 'categoryMap'},
      {set: {categoryMap: 'quarter'}, read: 'categoryMap'},
      {set: {categoryColumnName: 'race'}},
    ], 'render');
    expect(result).toEqual(['year', 'month', 'quarter']);
  });

  // #### Row source
  await softStep('Row source', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.dataFrame.selection.init((i: number) => i < 50);
      await (window as any).__viewerRendered('Pie chart');
    });
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {rowSource: 'Selected'}, read: 'rowSource'},
      {set: {rowSource: 'Filtered'}, read: 'rowSource'},
      {set: {rowSource: 'All'}, read: 'rowSource'},
    ], 'render');
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    expect(result).toEqual(['Selected', 'Filtered', 'All']);
  });

  // #### Aggregation functions
  await softStep('Aggregation functions', async () => {
    const steps = ['avg', 'min', 'max', 'sum', 'med', 'stdev', 'count']
      .map((aggr) => ({set: {segmentAngleAggrType: aggr}, read: 'segmentAngleAggrType'}));
    const result = await v.setViewerProps(page, 'Pie chart',
      [{set: {segmentAngleColumnName: 'age'}}, ...steps, {set: {segmentAngleColumnName: ''}}], 'render');
    expect(result).toEqual(['avg', 'min', 'max', 'sum', 'med', 'stdev', 'count']);
  });

  // #### Title and description
  await softStep('Title and description', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      const r: any[] = [];
      pie.props.showTitle = true;
      r.push(pie.props.showTitle);
      pie.props.title = 'Demographics';
      r.push(pie.props.title);
      pie.props.description = 'By race';
      r.push(pie.props.description);
      pie.props.descriptionPosition = 'Top';
      await (window as any).__viewerRendered('Pie chart');
      r.push(pie.props.descriptionPosition);
      pie.props.descriptionVisibilityMode = 'Never';
      r.push(pie.props.descriptionVisibilityMode);
      pie.props.title = '';
      pie.props.showTitle = false;
      return r;
    });
    expect(result[0]).toBe(true);
    expect(result[1]).toBe('Demographics');
    expect(result[2]).toBe('By race');
    expect(result[3]).toBe('Top');
    expect(result[4]).toBe('Never');
  });

  // #### Layout persistence
  await softStep('Layout persistence', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const pie = Array.from(tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'race';
      pie.props.segmentAngleColumnName = 'age';
      pie.props.startAngle = 45;
      pie.props.shift = 5;
      await (window as any).__viewerRendered('Pie chart');

      const before = {
        cat: pie.props.categoryColumnName,
        angle: pie.props.segmentAngleColumnName,
        startAngle: pie.props.startAngle,
        shift: pie.props.shift,
      };

      // The saved layout goes through the same ViewLayout serialization dapi.layouts would
      // persist; only the storage round-trip is replaced, so the viewer state under test is not.
      const json = tv.saveLayout().toJson();
      pie.close();
      await new Promise((res) => setTimeout(res, 100));
      tv.loadLayout(DG.ViewLayout.fromJson(json));
      await new Promise((res) => setTimeout(res, 500));

      const pie2 = Array.from(tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      const after = pie2 ? {
        cat: pie2.props.categoryColumnName,
        angle: pie2.props.segmentAngleColumnName,
        startAngle: pie2.props.startAngle,
        shift: pie2.props.shift,
      } : null;
      return {before, after};
    });
    expect(result.after).toEqual(result.before);
  });

  // #### Selection and interaction
  await softStep('Selection and interaction', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'race';
      pie.props.segmentAngleColumnName = '';
      await (window as any).__viewerRendered('Pie chart');
      const df = grok.shell.tv.dataFrame;

      const canvas = document.querySelector('[name="viewer-Pie-chart"] canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const x = rect.left + rect.width * 0.65;
      const y = rect.top + rect.height * 0.4;
      const settled = (window as any).__viewerRendered('Pie chart');
      for (const t of ['mousedown', 'mouseup', 'click'])
        canvas.dispatchEvent(new MouseEvent(t, {bubbles: true, clientX: x, clientY: y}));
      await settled;
      const sel1 = df.selection.trueCount;

      pie.props.showSelectedRows = false;
      const sOff = pie.props.showSelectedRows;
      pie.props.showSelectedRows = true;
      const sOn = pie.props.showSelectedRows;
      pie.props.showMouseOverRowGroup = false;
      const mOff = pie.props.showMouseOverRowGroup;
      pie.props.showMouseOverRowGroup = true;
      const mOn = pie.props.showMouseOverRowGroup;

      df.selection.setAll(false);
      return {sel1, sOff, sOn, mOff, mOn};
    });
    expect(result.sel1).toBeGreaterThan(0);
    expect(result.sOff).toBe(false);
    expect(result.sOn).toBe(true);
    expect(result.mOff).toBe(false);
    expect(result.mOn).toBe(true);
  });

  // #### On Click modes
  await softStep('On Click modes', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'race';
      const df = grok.shell.tv.dataFrame;
      const r: any[] = [];

      const canvas = document.querySelector('[name="viewer-Pie-chart"] canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const click = async (px: number, py: number) => {
        const settled = (window as any).__viewerRendered('Pie chart');
        for (const t of ['mousedown', 'mouseup', 'click'])
          canvas.dispatchEvent(new MouseEvent(t, {bubbles: true, clientX: px, clientY: py}));
        await settled;
      };
      const x = rect.left + rect.width * 0.65;
      const y = rect.top + rect.height * 0.4;

      pie.props.onClick = 'Select';
      await click(x, y);
      r.push({mode: 'Select', selected: df.selection.trueCount});

      pie.props.onClick = 'Filter';
      await click(x, y);
      r.push({mode: 'Filter', filtered: df.filter.trueCount, total: df.rowCount});

      await click(rect.left + 5, rect.top + 5);
      r.push({mode: 'clear', filtered: df.filter.trueCount});

      pie.props.onClick = 'Select';
      df.selection.setAll(false);
      return r;
    });
    expect(result[0].selected).toBeGreaterThan(0);
    expect(result[1].filtered).toBeLessThan(result[1].total);
    expect(result[2].filtered).toBe(result[1].total);
  });

  // #### Selection between grid and pie chart
  await softStep('Selection between grid and pie chart', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.init((i: number) => i < 50);
      await (window as any).__viewerRendered('Pie chart');
      const gridSel = df.selection.trueCount;

      const canvas = document.querySelector('[name="viewer-Pie-chart"] canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const x = rect.left + rect.width * 0.65;
      const y = rect.top + rect.height * 0.4;
      const settled = (window as any).__viewerRendered('Pie chart');
      for (const t of ['mousedown', 'mouseup', 'click'])
        canvas.dispatchEvent(new MouseEvent(t, {bubbles: true, clientX: x, clientY: y}));
      await settled;
      const pieSel = df.selection.trueCount;

      df.selection.setAll(false);
      return {gridSel, pieSel};
    });
    expect(result.gridSel).toBe(50);
    expect(result.pieSel).toBeGreaterThan(0);
  });

  // #### Auto layout
  await softStep('Auto layout', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {autoLayout: false}, read: 'autoLayout'},
      {set: {marginLeft: 50}, read: 'marginLeft'},
      {set: {marginTop: 50}, read: 'marginTop'},
      {set: {autoLayout: true}, read: 'autoLayout'},
    ], 'render');
    expect(result).toEqual([false, 50, 50, true]);
  });

  // #### Table switching and row source
  await softStep('Table switching and row source', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tv = grok.shell.tv;
      const df2 = grok.data.demo.molecules(500);
      df2.name = 'SPGI';
      grok.shell.addTableView(df2);
      await new Promise((res) => setTimeout(res, 200));
      grok.shell.v = Array.from(grok.shell.views)
        .find((x: any) => x.type === 'TableView' && x.dataFrame.name !== 'SPGI') as any;

      const pie = Array.from(tv.viewers).find((x: any) => x.type === 'Pie chart') as any;
      const r: any[] = [];

      pie.dataFrame = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      await (window as any).__viewerRendered('Pie chart');
      r.push(pie.dataFrame.name);

      pie.dataFrame = df;
      await (window as any).__viewerRendered('Pie chart');
      r.push(pie.dataFrame.name);

      pie.props.categoryColumnName = 'race';
      pie.props.rowSource = 'Selected';
      df.selection.init((i: number) => i < 100);
      await (window as any).__viewerRendered('Pie chart');
      r.push({rowSource: pie.props.rowSource, selCount: df.selection.trueCount});

      pie.props.rowSource = 'Filtered';
      const fg = tv.getFiltersGroup();
      await new Promise((res) => setTimeout(res, 300));
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'race', selected: ['Asian']});
      await (window as any).__viewerRendered('Pie chart');
      r.push({rowSource: pie.props.rowSource, filtered: df.filter.trueCount, total: df.rowCount});

      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'race', selected: df.col('race').categories});
      df.selection.setAll(false);
      pie.props.rowSource = 'All';
      return r;
    });
    expect(result[0]).toBe('SPGI');
    expect(result[1]).toBe('demog');
    expect(result[2].rowSource).toBe('Selected');
    expect(result[2].selCount).toBe(100);
    expect(result[3].rowSource).toBe('Filtered');
    expect(result[3].filtered).toBeGreaterThan(0);
    expect(result[3].filtered).toBeLessThan(result[3].total);
  });

  // A render wait that times out has silently become a sleep; the whole point of the
  // event ladder is that it does not.
  expect(await v.renderWaitCapHits(page)).toBe(0);

  v.finishSpec();
});
