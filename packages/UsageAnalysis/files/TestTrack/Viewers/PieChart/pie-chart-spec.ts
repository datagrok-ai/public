/* ---
realizes: [piechart.int.in-viewer-column-selector-reconfigures-pie, piechart.int.mouseover-row-group-cross-highlight]
--- */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:DemoFiles/chem/SPGI.csv';

test('Pie chart tests', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'pie-chart', 'Pie-chart');

  await softStep('Sorting', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {categoryColumnName: 'RACE', pieSortType: 'by value'}, read: 'pieSortType'},
      {set: {pieSortOrder: 'desc'}, read: 'pieSortOrder'},
      {set: {pieSortOrder: 'asc'}, read: 'pieSortOrder'},
      {set: {pieSortType: 'by category'}, read: 'pieSortType'},
      {set: {pieSortOrder: 'asc'}, read: 'pieSortOrder'},
      {set: {pieSortOrder: 'desc'}, read: 'pieSortOrder'},
    ]);
    expect(result).toEqual(['by value', 'desc', 'asc', 'by category', 'asc', 'desc']);
  });

  await softStep('Appearance', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {startAngle: 90}, wait: 200, read: 'startAngle'},
      {set: {startAngle: 180}, wait: 200, read: 'startAngle'},
      {set: {startAngle: 0}, wait: 200, read: 'startAngle'},
      {set: {maxRadius: 100}, wait: 200, read: 'maxRadius'},
      {set: {maxRadius: 150}, wait: 200, read: 'maxRadius'},
      {set: {shift: 10}, wait: 200, read: 'shift'},
      {set: {shift: 0}, wait: 200, read: 'shift'},
    ]);
    expect(result).toEqual([90, 180, 0, 100, 150, 10, 0]);
  });

  await softStep('Labels', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: any[] = [];

      for (const pos of ['Inside', 'Outside', 'Auto']) {
        pie.props.labelPosition = pos;
        await new Promise(res => setTimeout(res, 200));
        r.push(pie.props.labelPosition);
      }

      pie.props.showLabel = false;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.showLabel);

      pie.props.showPercentage = false;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.showPercentage);

      pie.props.showValue = true;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.showValue);

      pie.props.showLabel = true;
      pie.props.showPercentage = true;

      return r;
    });
    expect(result).toEqual(['Inside', 'Outside', 'Auto', false, false, true]);
  });

  await softStep('Outline', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {outlineLineWidth: 5}, wait: 200, read: 'outlineLineWidth'},
      {set: {outlineLineWidth: 0}, wait: 200, read: 'outlineLineWidth'},
      {set: {outlineLineWidth: 1}, wait: 200, read: 'outlineLineWidth'},
    ]);
    expect(result).toEqual([5, 0, 1]);
  });

  await softStep('Column selector', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {showColumnSelector: false}, wait: 200, read: 'showColumnSelector'},
      {set: {showColumnSelector: true}, wait: 200, read: 'showColumnSelector'},
    ]);
    expect(result).toEqual([false, true]);
  });

  await softStep('In-viewer column selector re-pick', async () => {
    await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'RACE';
      pie.props.showColumnSelector = true;
      pie.props.legendVisibility = 'Always';
      await new Promise(res => setTimeout(res, 600));
    });
    const readState = () => page.evaluate(() => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const labels = Array.from(pie.root.querySelectorAll('[name="legend"] .d4-legend-item .d4-legend-value'))
        .map((e: any) => (e.textContent ?? '').trim());
      return {
        cat: pie.props.categoryColumnName,
        labels,
        sexCats: df.col('SEX').categories.slice(),
        raceCats: df.col('RACE').categories.slice(),
      };
    });
    const pick = (columnName: string) => v.pickColumnViaSelector(page, {
      comboboxSuffix: 'category', columnName,
      viewerType: 'Pie chart', propName: 'categoryColumnName',
      scopeSelector: '[name="viewer-Pie-chart"]', popupWaitStrategy: 'either',
    });
    await pick('SEX');
    const afterSex = await readState();
    await pick('RACE');
    const afterRace = await readState();
    expect(afterSex.cat).toBe('SEX');
    expect([...afterSex.labels].sort()).toEqual([...afterSex.sexCats].sort());
    expect(afterRace.cat).toBe('RACE');
    expect([...afterRace.labels].sort()).toEqual([...afterRace.raceCats].sort());
    await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      pie.props.legendVisibility = 'Auto';
      await new Promise(res => setTimeout(res, 300));
    });
  });

  await softStep('Legend', async () => {
    const result = await v.setViewerProps(page, 'Pie chart', [
      {set: {legendVisibility: 'Always'}, wait: 200, read: 'legendVisibility'},
      {set: {legendPosition: 'LeftTop'}, wait: 200, read: 'legendPosition'},
      {set: {legendPosition: 'RightBottom'}, wait: 200, read: 'legendPosition'},
      {set: {legendVisibility: 'Never'}, wait: 200, read: 'legendVisibility'},
      {set: {legendVisibility: 'Auto'}, wait: 200, read: 'legendVisibility'},
    ]);
    expect(result).toEqual(['Always', 'LeftTop', 'RightBottom', 'Never', 'Auto']);
  });

  await softStep('Row source', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      df.selection.init((i: number) => i < 50);
      await new Promise(res => setTimeout(res, 300));
      const r: any[] = [];

      pie.props.rowSource = 'Selected';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.rowSource);

      pie.props.rowSource = 'Filtered';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.rowSource);

      pie.props.rowSource = 'All';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.rowSource);

      df.selection.setAll(false);
      return r;
    });
    expect(result).toEqual(['Selected', 'Filtered', 'All']);
  });

  await softStep('Title and description', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: any[] = [];

      pie.props.showTitle = true;
      r.push(pie.props.showTitle);

      pie.props.title = 'Demographics';
      r.push(pie.props.title);

      pie.props.description = 'By race';
      r.push(pie.props.description);

      pie.props.descriptionPosition = 'Top';
      await new Promise(res => setTimeout(res, 200));
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

  await softStep('Layout persistence', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;

      pie.props.categoryColumnName = 'RACE';
      pie.props.segmentAngleColumnName = 'AGE';
      pie.props.startAngle = 45;
      pie.props.shift = 5;
      await new Promise(res => setTimeout(res, 500));

      const before = {
        cat: pie.props.categoryColumnName,
        angle: pie.props.segmentAngleColumnName,
        startAngle: pie.props.startAngle,
        shift: pie.props.shift,
      };

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(res => setTimeout(res, 1000));

      pie.close();
      await new Promise(res => setTimeout(res, 500));

      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise(res => setTimeout(res, 3000));

      const pie2 = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const after = pie2 ? {
        cat: pie2.props.categoryColumnName,
        angle: pie2.props.segmentAngleColumnName,
        startAngle: pie2.props.startAngle,
        shift: pie2.props.shift,
      } : null;

      await grok.dapi.layouts.delete(saved);
      return {before, after};
    });
    expect(result.after).toEqual(result.before);
  });

  await softStep('Selection and interaction', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'RACE';

      pie.props.showSelectedRows = false;
      const sOff = pie.props.showSelectedRows;
      pie.props.showSelectedRows = true;
      const sOn = pie.props.showSelectedRows;

      pie.props.showMouseOverRowGroup = false;
      const mOff = pie.props.showMouseOverRowGroup;
      pie.props.showMouseOverRowGroup = true;
      const mOn = pie.props.showMouseOverRowGroup;

      return {sOff, sOn, mOff, mOn};
    });
    expect(result.sOff).toBe(false);
    expect(result.sOn).toBe(true);
    expect(result.mOff).toBe(false);
    expect(result.mOn).toBe(true);
  });

  await softStep('Mouse-over row group cross-highlight', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      const race = df.col('RACE');
      pie.props.categoryColumnName = 'RACE';
      const snap = () => {
        const cv = pie.root.querySelector('canvas') as HTMLCanvasElement;
        const d = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
        const m: Record<number, number> = {};
        for (let i = 0; i < d.length; i += 4) {
          const k = (d[i] << 16) | (d[i + 1] << 8) | d[i + 2];
          m[k] = (m[k] ?? 0) + 1;
        }
        return m;
      };
      const diff = (a: Record<number, number>, b: Record<number, number>) => {
        let s = 0;
        for (const k of Object.keys(b)) s += Math.abs(b[+k] - (a[+k] ?? 0));
        for (const k of Object.keys(a)) if (!(k in b)) s += a[+k];
        return s;
      };
      const settled = async () => {
        let prev = snap();
        for (let i = 0; i < 8; i++) {
          await new Promise(res => setTimeout(res, 400));
          const cur = snap();
          if (diff(prev, cur) === 0) return cur;
          prev = cur;
        }
        return prev;
      };
      pie.props.showMouseOverRowGroup = false;
      await new Promise(res => setTimeout(res, 600));
      const baseOff = await settled();
      df.rows.highlight((i: number) => race.get(i) === 'Asian');
      await new Promise(res => setTimeout(res, 600));
      const offDelta = diff(baseOff, await settled());
      df.rows.highlight(null);
      await new Promise(res => setTimeout(res, 600));
      await settled();
      pie.props.showMouseOverRowGroup = true;
      await new Promise(res => setTimeout(res, 600));
      const baseOn = await settled();
      df.rows.highlight((i: number) => race.get(i) === 'Asian');
      await new Promise(res => setTimeout(res, 600));
      const onDelta = diff(baseOn, await settled());
      df.rows.highlight(null);
      await new Promise(res => setTimeout(res, 600));
      const clearDelta = diff(baseOn, await settled());
      return {offDelta, onDelta, clearDelta};
    });
    console.log(`Mouse-over row group px: offDelta=${result.offDelta} onDelta=${result.onDelta} clearDelta=${result.clearDelta}`);
    expect(result.offDelta).toBeLessThan(2000);
    expect(result.onDelta).toBeGreaterThan(20000);
    expect(result.clearDelta).toBeLessThan(2000);
  });

  await softStep('Auto layout', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: any[] = [];

      pie.props.autoLayout = false;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.autoLayout);

      pie.props.marginLeft = 50;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.marginLeft);

      pie.props.marginTop = 50;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.marginTop);

      pie.props.autoLayout = true;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.autoLayout);

      return r;
    });
    expect(result).toEqual([false, 50, 50, true]);
  });

  await softStep('Table switching and row source (SPGI)', async () => {
    const demogName = await page.evaluate(async (spgiPath) => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));

      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      df.name = 'demog';
      grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const df2 = await grok.dapi.files.readCsv(spgiPath);
      df2.name = 'SPGI';
      grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const views = Array.from(grok.shell.views).filter((v: any) => v.type === 'TableView');
      const demogView = views.find((v: any) => v.dataFrame.name !== 'SPGI') as any;
      if (demogView) grok.shell.v = demogView;
      await new Promise(r => setTimeout(r, 500));
      return df.name;
    }, spgiPath);

    await v.addViewerByIcon(page, 'pie-chart', 'Pie-chart');

    const switched = await page.evaluate(async (demogName) => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: string[] = [];

      pie.props.table = 'SPGI';
      await new Promise(res => setTimeout(res, 500));
      r.push(pie.dataFrame.name);

      pie.props.table = demogName;
      await new Promise(res => setTimeout(res, 500));
      r.push(pie.dataFrame.name);
      return r;
    }, demogName);
    expect(switched[0]).toBe('SPGI');
    expect(switched[1]).toBe(demogName);

    const selection = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      pie.props.rowSource = 'Selected';
      df.selection.init((i: number) => i < 100);
      await new Promise(res => setTimeout(res, 300));
      return {rowSource: pie.props.rowSource, selCount: df.selection.trueCount};
    });
    expect(selection.rowSource).toBe('Selected');
    expect(selection.selCount).toBe(100);

    await page.evaluate(() => { grok.shell.tv.getFiltersGroup(); });
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 15000});
    const filtered = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      pie.props.rowSource = 'Filtered';
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian']});
      await new Promise(res => setTimeout(res, 500));
      const r = {rowSource: pie.props.rowSource, filtered: df.filter.trueCount};

      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: df.col('RACE').categories});
      df.selection.setAll(false);
      pie.props.rowSource = 'All';
      return r;
    });
    expect(filtered.rowSource).toBe('Filtered');
    expect(filtered.filtered).toBeGreaterThan(0);
  });

  v.finishSpec();
});
