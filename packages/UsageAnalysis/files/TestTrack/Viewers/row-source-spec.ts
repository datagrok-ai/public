/* ---
realizes: [viewers.scatter-plot, viewers.line-chart, viewers.histogram, viewers.bar-chart, viewers.pie-chart, viewers.box-plot, viewers.pc-plot, viewers.filters]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../spec-login';
import * as v from '../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const demogName = 'demog';
const spgiName = 'spgi-100';

const allRowSources = ['Filtered', 'All', 'Selected', 'SelectedOrCurrent', 'FilteredSelected', 'MouseOverGroup', 'CurrentRow', 'MouseOverRow'] as const;

async function openDatasets(page: Page) {
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await page.evaluate(async (args: {spgiPath: string; demogName: string; spgiName: string}) => {
    const w = window as any;
    const demog = grok.shell.t;
    demog.name = args.demogName;
    const demogView = grok.shell.tv;
    const spgi = await w.__readCsv(args.spgiPath);
    spgi.name = args.spgiName;
    grok.shell.addTableView(spgi);
    await w.__tableReady(5000);
    grok.shell.v = demogView;
    await w.__poll(() => grok.shell.tv?.dataFrame?.name, (n: string) => n === args.demogName, 2000, 25);
  }, {spgiPath, demogName, spgiName});
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
}

async function testDemogRowSources(page: Page, viewerType: string, filterPanelSex: string) {
  const result = await page.evaluate(async (args: {viewerType: string; filterPanelSex: string}) => {
    const w = window as any;
    const tv = grok.shell.tv;
    const df = tv.dataFrame;
    const v = tv.viewers.find((vw: any) => vw.type === args.viewerType)!;
    const rendered = (act: () => void, capMs: number) =>
      w.__settled('viewer:' + args.viewerType + '.onViewerRendered', act, capMs);
    const r: any = {};

    const fg = tv.getFiltersGroup();
    await w.__filtered(() => fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: [args.filterPanelSex]}),
      300, df.filter.trueCount);
    r.filteredCount = df.filter.trueCount;
    r.rowSourceFiltered = v.props.rowSource;

    await w.__filtered(() => fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: df.col('SEX').categories}),
      200, df.filter.trueCount);

    await rendered(() => { v.props.rowSource = 'All'; }, 200);
    r.rsAll = v.props.rowSource;

    df.selection.setAll(false);
    await rendered(() => { v.props.rowSource = 'Selected'; }, 200);
    r.rsSelected = v.props.rowSource;
    r.selEmpty = df.selection.trueCount;

    await rendered(() => df.selection.init((i: number) => i < 100), 200);
    const ageCol = df.col('AGE');
    let selWithAge44 = 0;
    for (let i = 0; i < 100; i++)
      if (ageCol.get(i) > 44) selWithAge44++;
    r.selCount = df.selection.trueCount;
    r.selWithAgeGt44 = selWithAge44;

    await rendered(() => { v.props.rowSource = 'SelectedOrCurrent'; }, 200);
    r.rsSelectedOrCurrent = v.props.rowSource;

    df.selection.setAll(false);
    let curIdx = -1;
    for (let i = 0; i < df.rowCount; i++)
      if (ageCol.get(i) > 44) { curIdx = i; break; }
    await rendered(() => { df.currentRowIdx = curIdx; }, 200);
    r.currentAge = ageCol.get(df.currentRowIdx);

    await rendered(() => {
      v.props.rowSource = 'FilteredSelected';
      df.selection.setAll(false);
    }, 200);
    r.rsFilteredSelected = v.props.rowSource;

    await rendered(() => df.selection.init((i: number) => {
      const a = ageCol.get(i);
      return a >= 42 && a <= 47;
    }), 200);
    r.filtSelCount = df.selection.trueCount;
    let ageGt44InSel = 0;
    for (let i = 0; i < df.rowCount; i++) {
      const a = ageCol.get(i);
      if (a >= 42 && a <= 47 && df.filter.get(i) && a > 44) ageGt44InSel++;
    }
    r.filtSelAgeGt44 = ageGt44InSel;
    df.selection.setAll(false);

    await rendered(() => { v.props.rowSource = 'MouseOverGroup'; }, 200);
    r.rsMouseOverGroup = v.props.rowSource;

    await rendered(() => {
      df.mouseOverGroup = DG.BitSet.create(df.rowCount, (i: number) => df.col('RACE').get(i) === 'Asian');
    }, 300);
    r.mogCount = df.mouseOverGroup.trueCount;
    df.mouseOverGroup = null;

    await rendered(() => {
      v.props.rowSource = 'CurrentRow';
      df.currentRowIdx = curIdx;
    }, 200);
    r.rsCurrentRow = v.props.rowSource;

    await rendered(() => { v.props.rowSource = 'MouseOverRow'; }, 200);
    r.rsMouseOverRow = v.props.rowSource;

    return r;
  }, {viewerType, filterPanelSex});

  expect(result.rowSourceFiltered).toBe('Filtered');
  expect(result.rsAll).toBe('All');
  expect(result.rsSelected).toBe('Selected');
  expect(result.selEmpty).toBe(0);
  expect(result.selCount).toBe(100);
  expect(result.rsSelectedOrCurrent).toBe('SelectedOrCurrent');
  expect(result.rsFilteredSelected).toBe('FilteredSelected');
  expect(result.rsMouseOverGroup).toBe('MouseOverGroup');
  expect(result.mogCount).toBeGreaterThan(0);
  expect(result.rsCurrentRow).toBe('CurrentRow');
  expect(result.rsMouseOverRow).toBe('MouseOverRow');

  return result;
}

async function testSpgiRowSources(page: Page, viewerType: string, colConfig: Record<string, any>) {
  const result = await page.evaluate(async (args: {viewerType: string; colConfig: Record<string, any>; spgiName: string}) => {
    const w = window as any;
    const tv = grok.shell.tv;
    const v = tv.viewers.find((vw: any) => vw.type === args.viewerType)!;
    const rendered = (act: () => void, capMs: number) =>
      w.__settled('viewer:' + args.viewerType + '.onViewerRendered', act, capMs);
    const spgiDf = (Array.from(grok.shell.tableViews) as any[])
      .find((vw: any) => vw.dataFrame?.name === args.spgiName)?.dataFrame;
    if (!spgiDf) return {error: 'spgi not found'};

    await rendered(() => { v.props.table = spgiDf.name; }, 500);

    await rendered(() => {
      for (const [key, val] of Object.entries(args.colConfig))
        v.props[key] = val;
      v.props.rowSource = 'Filtered';
      v.props.filter = '${Stereo Category} in ["R_ONE", "S_UNKN"]';
    }, 500);

    const r: any = {tableSet: v.dataFrame?.name};
    for (const rs of ['Filtered', 'All', 'Selected', 'SelectedOrCurrent', 'FilteredSelected', 'MouseOverGroup', 'CurrentRow', 'MouseOverRow']) {
      v.props.rowSource = rs;
      r['rs_' + rs] = v.props.rowSource;
    }

    v.props.rowSource = 'MouseOverGroup';
    const scCol = spgiDf.col('Stereo Category');
    await rendered(() => {
      spgiDf.mouseOverGroup = DG.BitSet.create(spgiDf.rowCount, (i: number) => scCol.get(i) === 'R_ONE');
    }, 300);
    r.mogCount = spgiDf.mouseOverGroup.trueCount;
    spgiDf.mouseOverGroup = null;

    v.close();
    r.closed = await w.__poll(() => !tv.viewers.some((vw: any) => vw.type === args.viewerType), (c: boolean) => c, 300, 25);

    return r;
  }, {viewerType, colConfig, spgiName});

  expect(result.tableSet).toBeTruthy();
  for (const rs of allRowSources)
    expect(result['rs_' + rs]).toBe(rs);
  expect(result.mogCount).toBeGreaterThan(0);
  expect(result.closed).toBe(true);

  return result;
}

async function addViewerWithFilter(page: Page, viewerType: string, props: Record<string, any>, filter: string) {
  const result = await page.evaluate(async (args: {viewerType: string; props: Record<string, any>; filter: string}) => {
    const w = window as any;
    const tv = grok.shell.tv;
    const v = tv.addViewer(args.viewerType);
    const rendered = (act: () => void, capMs: number) =>
      w.__settled('viewer:' + args.viewerType + '.onViewerRendered', act, capMs);
    await rendered(() => {
      for (const [key, val] of Object.entries(args.props))
        v.props[key] = val;
    }, 500);
    await rendered(() => { v.props.filter = args.filter; }, 300);
    return {type: v.type, rowSource: v.props.rowSource, filter: v.props.filter};
  }, {viewerType, props, filter});

  expect(result.type).toBe(viewerType);
  expect(result.rowSource).toBe('Filtered');
  expect(result.filter).toBe(filter);
}

test('Row Source tests', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await openDatasets(page);

  await page.evaluate(() => {
    grok.shell.tv.getFiltersGroup();
  });
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  await softStep('Scatter Plot row sources (demog)', async () => {
    await addViewerWithFilter(page, 'Scatter plot',
      {xColumnName: 'AGE', yColumnName: 'HEIGHT', colorColumnName: 'RACE'},
      '${AGE} > 44');
    await testDemogRowSources(page, 'Scatter plot', 'M');
  });

  await softStep('Scatter Plot row sources (spgi-100)', async () => {
    await testSpgiRowSources(page, 'Scatter plot',
      {xColumnName: 'Chemical Space X', yColumnName: 'Chemical Space Y', colorColumnName: 'Stereo Category'});
  });

  await softStep('Line Chart row sources (demog)', async () => {
    await addViewerWithFilter(page, 'Line chart',
      {xColumnName: 'AGE', yColumnNames: ['HEIGHT']},
      '${AGE} > 44');
    await testDemogRowSources(page, 'Line chart', 'F');
  });

  await softStep('Line Chart row sources (spgi-100)', async () => {
    await testSpgiRowSources(page, 'Line chart',
      {xColumnName: 'Chemical Space X', yColumnNames: ['TPSA']});
  });

  await softStep('Histogram row sources (demog)', async () => {
    await addViewerWithFilter(page, 'Histogram',
      {valueColumnName: 'AGE'},
      '${AGE} > 44');
    await testDemogRowSources(page, 'Histogram', 'M');
  });

  await softStep('Histogram row sources (spgi-100)', async () => {
    await testSpgiRowSources(page, 'Histogram',
      {valueColumnName: 'TPSA'});
  });

  await softStep('Bar Chart row sources (demog)', async () => {
    await addViewerWithFilter(page, 'Bar chart',
      {valueColumnName: 'AGE', splitColumnName: 'RACE'},
      '${AGE} > 44');
    await testDemogRowSources(page, 'Bar chart', 'M');
  });

  await softStep('Bar Chart row sources (spgi-100)', async () => {
    await testSpgiRowSources(page, 'Bar chart',
      {valueColumnName: 'TPSA', splitColumnName: 'Stereo Category'});
  });

  await softStep('Pie Chart row sources (demog)', async () => {
    await addViewerWithFilter(page, 'Pie chart',
      {categoryColumnName: 'RACE'},
      '${AGE} > 44');
    await testDemogRowSources(page, 'Pie chart', 'M');
  });

  await softStep('Pie Chart row sources (spgi-100)', async () => {
    await testSpgiRowSources(page, 'Pie chart',
      {categoryColumnName: 'Stereo Category'});
  });

  await softStep('Box Plot row sources (demog)', async () => {
    await addViewerWithFilter(page, 'Box plot',
      {valueColumnName: 'AGE', category1ColumnName: 'RACE'},
      '${AGE} > 44');
    await testDemogRowSources(page, 'Box plot', 'F');
  });

  await softStep('Box Plot row sources (spgi-100)', async () => {
    await testSpgiRowSources(page, 'Box plot',
      {valueColumnName: 'TPSA', category1ColumnName: 'Stereo Category'});
  });

  await softStep('PC Plot row sources (demog)', async () => {
    await addViewerWithFilter(page, 'PC Plot',
      {columnNames: ['AGE', 'HEIGHT', 'WEIGHT']},
      '${AGE} > 44');
    await testDemogRowSources(page, 'PC Plot', 'M');
  });

  await softStep('PC Plot row sources (spgi-100)', async () => {
    await testSpgiRowSources(page, 'PC Plot',
      {columnNames: ['Chemical Space X', 'Chemical Space Y', 'TPSA']});
  });

  await v.cleanupShell(page);

  v.finishSpec();
});
