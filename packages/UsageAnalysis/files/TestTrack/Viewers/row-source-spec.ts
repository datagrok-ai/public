import {test, expect, Page} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.BASE_URL ?? 'https://dev.datagrok.ai';
const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

const allRowSources = ['Filtered', 'All', 'Selected', 'SelectedOrCurrent', 'FilteredSelected', 'MouseOverGroup', 'CurrentRow', 'MouseOverRow'] as const;

async function openDatasetsAndSetup(page: Page) {
  await page.waitForFunction(() => {
    try { return typeof grok !== 'undefined' && grok.shell &&
      typeof grok.shell.settings?.showFiltersIconsConstantly === 'boolean'; }
    catch (e) { return false; }
  }, {timeout: 30000});

  await page.evaluate(async (args: {demogPath: string; spgiPath: string}) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
    grok.shell.closeAll();

    // Open demog
    const df1 = await grok.dapi.files.readCsv(args.demogPath);
    const tv1 = grok.shell.addTableView(df1);
    await new Promise(resolve => {
      const sub = df1.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });

    // Open spgi-100
    const df2 = await grok.dapi.files.readCsv(args.spgiPath);
    const tv2 = grok.shell.addTableView(df2);
    await new Promise(resolve => {
      const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
    const hasBioChem = Array.from({length: df2.columns.length}, (_, i) => df2.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }

    // Switch back to demog
    const views = Array.from(grok.shell.views);
    const demogView = views.find((v: any) => v.name === 'Table');
    if (demogView) grok.shell.v = demogView;
    await new Promise(r => setTimeout(r, 500));
  }, {demogPath, spgiPath});

  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
}

/** Test all row sources for a viewer on demog, including filter panel, selection, and mouseOverGroup */
async function testDemogRowSources(page: Page, viewerType: string, filterPanelSex: string) {
  const result = await page.evaluate(async (args: {viewerType: string; filterPanelSex: string}) => {
    const tv = grok.shell.tv;
    const df = tv.dataFrame;
    const v = tv.viewers.find((vw: any) => vw.type === args.viewerType)!;
    const r: any = {};

    // Step 3.1: Filtered + Filter Panel
    const fg = tv.getFiltersGroup();
    await new Promise(res => setTimeout(res, 300));
    fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: [args.filterPanelSex]});
    await new Promise(res => setTimeout(res, 300));
    r.filteredCount = df.filter.trueCount;
    r.rowSourceFiltered = v.props.rowSource;
    // Reset filter panel
    fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: df.col('SEX').categories});
    await new Promise(res => setTimeout(res, 200));

    // Step 3.2: All
    v.props.rowSource = 'All';
    await new Promise(res => setTimeout(res, 200));
    r.rsAll = v.props.rowSource;

    // Step 3.3: Selected (empty)
    df.selection.setAll(false);
    v.props.rowSource = 'Selected';
    await new Promise(res => setTimeout(res, 200));
    r.rsSelected = v.props.rowSource;
    r.selEmpty = df.selection.trueCount;

    // Step 3.4: Select 100 rows
    df.selection.init((i: number) => i < 100);
    await new Promise(res => setTimeout(res, 200));
    const ageCol = df.col('AGE');
    let selWithAge44 = 0;
    for (let i = 0; i < 100; i++)
      if (ageCol.get(i) > 44) selWithAge44++;
    r.selCount = df.selection.trueCount;
    r.selWithAgeGt44 = selWithAge44;

    // Step 3.5: SelectedOrCurrent
    v.props.rowSource = 'SelectedOrCurrent';
    await new Promise(res => setTimeout(res, 200));
    r.rsSelectedOrCurrent = v.props.rowSource;

    // Step 3.6: Reset selection, set current row
    df.selection.setAll(false);
    let curIdx = -1;
    for (let i = 0; i < df.rowCount; i++)
      if (ageCol.get(i) > 44) { curIdx = i; break; }
    df.currentRowIdx = curIdx;
    await new Promise(res => setTimeout(res, 200));
    r.currentAge = ageCol.get(df.currentRowIdx);

    // Step 3.7: FilteredSelected (empty)
    v.props.rowSource = 'FilteredSelected';
    df.selection.setAll(false);
    await new Promise(res => setTimeout(res, 200));
    r.rsFilteredSelected = v.props.rowSource;

    // Step 3.8: Select AGE 42-47
    df.selection.init((i: number) => {
      const a = ageCol.get(i);
      return a >= 42 && a <= 47;
    });
    await new Promise(res => setTimeout(res, 200));
    r.filtSelCount = df.selection.trueCount;
    let ageGt44InSel = 0;
    for (let i = 0; i < df.rowCount; i++) {
      const a = ageCol.get(i);
      if (a >= 42 && a <= 47 && df.filter.get(i) && a > 44) ageGt44InSel++;
    }
    r.filtSelAgeGt44 = ageGt44InSel;
    df.selection.setAll(false);

    // Step 3.9: MouseOverGroup
    v.props.rowSource = 'MouseOverGroup';
    await new Promise(res => setTimeout(res, 200));
    r.rsMouseOverGroup = v.props.rowSource;

    // Step 3.10: MouseOverGroup programmatic
    df.mouseOverGroup = DG.BitSet.create(df.rowCount, (i: number) => df.col('RACE').get(i) === 'Asian');
    await new Promise(res => setTimeout(res, 300));
    r.mogCount = df.mouseOverGroup.trueCount;
    df.mouseOverGroup = null;

    // Step 3.11: CurrentRow
    v.props.rowSource = 'CurrentRow';
    df.currentRowIdx = curIdx;
    await new Promise(res => setTimeout(res, 200));
    r.rsCurrentRow = v.props.rowSource;

    // Step 3.12: MouseOverRow
    v.props.rowSource = 'MouseOverRow';
    await new Promise(res => setTimeout(res, 200));
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

/** Switch viewer to spgi-100, configure, and test all row sources */
async function testSpgiRowSources(page: Page, viewerType: string, colConfig: Record<string, any>) {
  const result = await page.evaluate(async (args: {viewerType: string; colConfig: Record<string, any>}) => {
    const tv = grok.shell.tv;
    const v = tv.viewers.find((vw: any) => vw.type === args.viewerType)!;
    const views = Array.from(grok.shell.views);
    const spgiView = views.find((vw: any) => vw.name === 'Table (2)');
    const spgiDf = (spgiView as any)?.table;
    if (!spgiDf) return {error: 'spgi not found'};

    v.props.table = spgiDf.name;
    await new Promise(r => setTimeout(r, 500));

    for (const [key, val] of Object.entries(args.colConfig))
      v.props[key] = val;
    v.props.rowSource = 'Filtered';
    v.props.filter = '${Stereo Category} in ["R_ONE", "S_UNKN"]';
    await new Promise(r => setTimeout(r, 500));

    const r: any = {tableSet: v.dataFrame?.name};
    for (const rs of ['Filtered', 'All', 'Selected', 'SelectedOrCurrent', 'FilteredSelected', 'MouseOverGroup', 'CurrentRow', 'MouseOverRow']) {
      v.props.rowSource = rs;
      r['rs_' + rs] = v.props.rowSource;
    }

    // MouseOverGroup programmatic
    v.props.rowSource = 'MouseOverGroup';
    const scCol = spgiDf.col('Stereo Category');
    spgiDf.mouseOverGroup = DG.BitSet.create(spgiDf.rowCount, (i: number) => scCol.get(i) === 'R_ONE');
    await new Promise(res => setTimeout(res, 300));
    r.mogCount = spgiDf.mouseOverGroup.trueCount;
    spgiDf.mouseOverGroup = null;

    // Close viewer
    v.close();
    await new Promise(res => setTimeout(res, 300));
    r.closed = !tv.viewers.some((vw: any) => vw.type === args.viewerType);

    // Switch back to demog
    const demogView = Array.from(grok.shell.views).find((vw: any) => vw.name === 'Table');
    if (demogView) grok.shell.v = demogView;
    await new Promise(r => setTimeout(r, 200));

    return r;
  }, {viewerType, colConfig});

  expect(result.tableSet).toBeTruthy();
  for (const rs of allRowSources)
    expect(result['rs_' + rs]).toBe(rs);
  expect(result.mogCount).toBeGreaterThan(0);
  expect(result.closed).toBe(true);

  return result;
}

/** Add a viewer with given props, set filter, and return it */
async function addViewerWithFilter(page: Page, viewerType: string, props: Record<string, any>, filter: string) {
  const result = await page.evaluate(async (args: {viewerType: string; props: Record<string, any>; filter: string}) => {
    const tv = grok.shell.tv;
    const v = tv.addViewer(args.viewerType);
    for (const [key, val] of Object.entries(args.props))
      v.props[key] = val;
    await new Promise(r => setTimeout(r, 500));
    v.props.filter = args.filter;
    await new Promise(r => setTimeout(r, 300));
    return {type: v.type, rowSource: v.props.rowSource, filter: v.props.filter};
  }, {viewerType, props, filter});

  expect(result.type).toBe(viewerType);
  expect(result.rowSource).toBe('Filtered');
  expect(result.filter).toBe(filter);
}

test('Row Source tests', async ({page}) => {
  test.setTimeout(900_000);

  // Phase 1: Navigate
  await page.goto(baseUrl);
  await openDatasetsAndSetup(page);

  // Phase 2: Open filter panel on demog
  await page.evaluate(() => {
    grok.shell.tv.getFiltersGroup();
  });
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  // #### Scatter Plot
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

  // #### Line Chart
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

  // #### Histogram
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

  // #### Bar Chart
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

  // #### Pie Chart
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

  // #### Box Plot
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

  // #### PC Plot
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

  // #### Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
