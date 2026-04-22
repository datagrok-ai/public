import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';
const datasetPath = 'System:DemoFiles/demog.csv';
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

test('Box plot tests', async ({page}) => {
  test.setTimeout(600_000);

  // Phase 1: Navigate + login
  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Add Box Plot via toolbox icon
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-box-plot"]');
    if (icon) (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});

  // #### Plot style: box vs violin
  await softStep('Plot style: box vs violin', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.valueColumnName = 'AGE';
      bp.props.category1ColumnName = 'RACE';
      const defaultStyle = bp.props.plotStyle;
      bp.props.plotStyle = 'violin';
      const violin = bp.props.plotStyle;
      bp.props.bins = 50;
      const bins50 = bp.props.bins;
      bp.props.bins = 500;
      const bins500 = bp.props.bins;
      bp.props.interquartileLineWidth = 10;
      const iqw = bp.props.interquartileLineWidth;
      bp.props.plotStyle = 'box';
      bp.props.bins = 100; bp.props.interquartileLineWidth = 6;
      return {defaultStyle, violin, bins50, bins500, iqw, final: bp.props.plotStyle};
    });
    expect(result.defaultStyle).toBe('box');
    expect(result.violin).toBe('violin');
    expect(result.bins50).toBe(50);
    expect(result.bins500).toBe(500);
    expect(result.iqw).toBe(10);
    expect(result.final).toBe('box');
  });

  // #### Two-level categories
  await softStep('Two-level categories', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.category1ColumnName = 'SEX';
      bp.props.category2ColumnName = 'RACE';
      const cat2 = bp.props.category2ColumnName;
      bp.props.showMinorCategories = false;
      const minorOff = bp.props.showMinorCategories;
      bp.props.showMinorCategories = true;
      const minorOn = bp.props.showMinorCategories;
      bp.props.showAllCategories = true;
      const showAllOn = bp.props.showAllCategories;
      bp.props.showAllCategories = false;
      bp.props.category2ColumnName = '';
      return {cat2, minorOff, minorOn, showAllOn};
    });
    expect(result.cat2).toBe('RACE');
    expect(result.minorOff).toBe(false);
    expect(result.minorOn).toBe(true);
    expect(result.showAllOn).toBe(true);
  });

  // #### Statistics display
  await softStep('Statistics display', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      const showStatsDefault = bp.props.showStatistics;
      bp.props.showTotalCount = true;
      bp.props.showInliersCount = true;
      bp.props.showOutliersCount = true;
      bp.props.showStdev = true;
      bp.props.showQ1 = true;
      bp.props.showQ3 = true;
      const r = {
        showStatsDefault,
        totalCount: bp.props.showTotalCount,
        inliers: bp.props.showInliersCount,
        outliers: bp.props.showOutliersCount,
        stdev: bp.props.showStdev,
        q1: bp.props.showQ1,
        q3: bp.props.showQ3
      };
      bp.props.showStatistics = false;
      r.statsOff = bp.props.showStatistics;
      bp.props.showStatistics = true;
      bp.props.showTotalCount = false; bp.props.showInliersCount = false;
      bp.props.showOutliersCount = false; bp.props.showStdev = false;
      bp.props.showQ1 = false; bp.props.showQ3 = false;
      return r;
    });
    expect(result.showStatsDefault).toBe(true);
    expect(result.totalCount).toBe(true);
    expect(result.statsOff).toBe(false);
  });

  // #### Markers
  await softStep('Markers', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.valueColumnName = 'AGE';
      bp.props.category1ColumnName = 'SEX';
      bp.props.markerType = 'square';
      const markerType = bp.props.markerType;
      bp.props.markerSize = 10;
      const markerSize = bp.props.markerSize;
      bp.props.markerOpacity = 80;
      const markerOpacity = bp.props.markerOpacity;
      bp.props.markersColumnName = 'RACE';
      const markersCol = bp.props.markersColumnName;
      bp.props.markerSizeColumnName = 'WEIGHT';
      const markerSizeCol = bp.props.markerSizeColumnName;
      bp.props.markersColumnName = '';
      bp.props.markerSizeColumnName = '';
      bp.props.markerType = 'circle';
      bp.props.markerSize = null;
      bp.props.markerOpacity = 35;
      return {markerType, markerSize, markerOpacity, markersCol, markerSizeCol};
    });
    expect(result.markerType).toBe('square');
    expect(result.markerSize).toBe(10);
    expect(result.markerOpacity).toBe(80);
    expect(result.markersCol).toBe('RACE');
    expect(result.markerSizeCol).toBe('WEIGHT');
  });

  // #### Marker and Bin color coding
  await softStep('Marker and Bin color coding', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.markerColorColumnName = 'AGE';
      const markerColorAGE = bp.props.markerColorColumnName;
      bp.props.colorAxisType = 'logarithmic';
      const colorAxisLog = bp.props.colorAxisType;
      bp.props.invertColorScheme = true;
      const invertOn = bp.props.invertColorScheme;
      bp.props.markerColorColumnName = 'RACE';
      const markerColorRACE = bp.props.markerColorColumnName;
      bp.props.invertColorScheme = false;
      bp.props.markerColorColumnName = '';
      bp.props.binColorColumnName = 'WEIGHT';
      const binColorWeight = bp.props.binColorColumnName;
      bp.props.binColorAggrType = 'min';
      const aggrMin = bp.props.binColorAggrType;
      bp.props.binColorAggrType = 'max';
      const aggrMax = bp.props.binColorAggrType;
      bp.props.binColorAggrType = 'med';
      const aggrMed = bp.props.binColorAggrType;
      bp.props.binColorColumnName = '';
      bp.props.colorAxisType = 'linear';
      return {markerColorAGE, colorAxisLog, invertOn, markerColorRACE,
        binColorWeight, aggrMin, aggrMax, aggrMed};
    });
    expect(result.markerColorAGE).toBe('AGE');
    expect(result.colorAxisLog).toBe('logarithmic');
    expect(result.invertOn).toBe(true);
    expect(result.binColorWeight).toBe('WEIGHT');
    expect(result.aggrMin).toBe('min');
    expect(result.aggrMax).toBe('max');
    expect(result.aggrMed).toBe('med');
  });

  // #### Value axis configuration
  await softStep('Value axis configuration', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.valueColumnName = 'AGE';
      bp.props.axisType = 'logarithmic';
      const axisLog = bp.props.axisType;
      bp.props.invertYAxis = true;
      const invertY = bp.props.invertYAxis;
      bp.props.valueMin = 20;
      const min20 = bp.props.valueMin;
      bp.props.valueMax = 60;
      const max60 = bp.props.valueMax;
      bp.props.valueMin = null; bp.props.valueMax = null;
      bp.props.axisType = 'linear'; bp.props.invertYAxis = false;
      return {axisLog, invertY, min20, max60};
    });
    expect(result.axisLog).toBe('logarithmic');
    expect(result.invertY).toBe(true);
    expect(result.min20).toBe(20);
    expect(result.max60).toBe(60);
  });

  // #### Zoom by filter and Show empty categories
  await softStep('Zoom by filter and empty categories', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.valueColumnName = 'AGE';
      const zoomDefault = bp.props.zoomValuesByFilter;
      bp.props.zoomValuesByFilter = false;
      const zoomOff = bp.props.zoomValuesByFilter;
      bp.props.zoomValuesByFilter = true;
      bp.props.category1ColumnName = 'RACE';
      bp.props.showEmptyCategories = false;
      const emptyOff = bp.props.showEmptyCategories;
      bp.props.showEmptyCategories = true;
      const emptyOn = bp.props.showEmptyCategories;
      return {zoomDefault, zoomOff, emptyOff, emptyOn};
    });
    expect(result.zoomDefault).toBe(true);
    expect(result.zoomOff).toBe(false);
    expect(result.emptyOff).toBe(false);
    expect(result.emptyOn).toBe(true);
  });

  // #### Box plot components
  await softStep('Box plot components toggles', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.showMeanCross = false;
      const meanOff = bp.props.showMeanCross;
      bp.props.showMedianDash = false;
      const medianOff = bp.props.showMedianDash;
      bp.props.showUpperDash = false;
      const upperOff = bp.props.showUpperDash;
      bp.props.showLowerDash = false;
      const lowerOff = bp.props.showLowerDash;
      bp.props.showInsideValues = false;
      const insideOff = bp.props.showInsideValues;
      bp.props.showOutsideValues = false;
      const outsideOff = bp.props.showOutsideValues;
      bp.props.showMeanCross = true; bp.props.showMedianDash = true;
      bp.props.showUpperDash = true; bp.props.showLowerDash = true;
      bp.props.showInsideValues = true; bp.props.showOutsideValues = true;
      return {meanOff, medianOff, upperOff, lowerOff, insideOff, outsideOff};
    });
    expect(result.meanOff).toBe(false);
    expect(result.medianOff).toBe(false);
    expect(result.upperOff).toBe(false);
    expect(result.lowerOff).toBe(false);
    expect(result.insideOff).toBe(false);
    expect(result.outsideOff).toBe(false);
  });

  // #### Controls visibility
  await softStep('Controls visibility', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.showValueSelector = false;
      const valSelOff = bp.props.showValueSelector;
      bp.props.showCategorySelector = false;
      const catSelOff = bp.props.showCategorySelector;
      bp.props.showColorSelector = false;
      const colorSelOff = bp.props.showColorSelector;
      bp.props.showSizeSelector = false;
      const sizeSelOff = bp.props.showSizeSelector;
      bp.props.showValueAxis = false;
      const valAxisOff = bp.props.showValueAxis;
      bp.props.showCategoryAxis = false;
      const catAxisOff = bp.props.showCategoryAxis;
      bp.props.showValueSelector = true; bp.props.showCategorySelector = true;
      bp.props.showColorSelector = true; bp.props.showSizeSelector = true;
      bp.props.showValueAxis = true; bp.props.showCategoryAxis = true;
      return {valSelOff, catSelOff, colorSelOff, sizeSelOff, valAxisOff, catAxisOff};
    });
    expect(result.valSelOff).toBe(false);
    expect(result.catSelOff).toBe(false);
    expect(result.colorSelOff).toBe(false);
    expect(result.sizeSelOff).toBe(false);
    expect(result.valAxisOff).toBe(false);
    expect(result.catAxisOff).toBe(false);
  });

  // #### Title and description
  await softStep('Title and description', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.showTitle = true;
      bp.props.title = 'Age by Race';
      bp.props.description = 'Box plot of patient ages';
      bp.props.descriptionVisibilityMode = 'Always';
      bp.props.descriptionPosition = 'Bottom';
      const r = {
        showTitle: bp.props.showTitle, title: bp.props.title,
        description: bp.props.description, descVisMode: bp.props.descriptionVisibilityMode,
        descPos: bp.props.descriptionPosition
      };
      bp.props.descriptionVisibilityMode = 'Never';
      r.descVisModeNever = bp.props.descriptionVisibilityMode;
      bp.props.showTitle = false; bp.props.title = ''; bp.props.description = '';
      bp.props.descriptionVisibilityMode = 'Auto';
      return r;
    });
    expect(result.showTitle).toBe(true);
    expect(result.title).toBe('Age by Race');
    expect(result.description).toBe('Box plot of patient ages');
    expect(result.descVisMode).toBe('Always');
    expect(result.descPos).toBe('Bottom');
    expect(result.descVisModeNever).toBe('Never');
  });

  // #### Date category mapping
  await softStep('Date category mapping', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.category1ColumnName = 'STARTED';
      const catStarted = bp.props.category1ColumnName;
      bp.props.category1Map = 'Month';
      const mapMonth = bp.props.category1Map;
      bp.props.category1Map = 'Quarter';
      const mapQuarter = bp.props.category1Map;
      bp.props.category1ColumnName = 'RACE';
      return {catStarted, mapMonth, mapQuarter};
    });
    expect(result.catStarted).toBe('STARTED');
    expect(result.mapMonth).toBe('Month');
    expect(result.mapQuarter).toBe('Quarter');
  });

  // #### Style customization
  await softStep('Style customization', async () => {
    const result = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.whiskerLineWidth = 4;
      const whiskerWidth = bp.props.whiskerLineWidth;
      bp.props.whiskerWidthRatio = 1.0;
      const ratio1 = bp.props.whiskerWidthRatio;
      bp.props.whiskerWidthRatio = 0.3;
      const ratio03 = bp.props.whiskerWidthRatio;
      bp.props.autoLayout = false;
      const autoLayoutOff = bp.props.autoLayout;
      bp.props.axisUseColumnFormat = false;
      const axisFormatOff = bp.props.axisUseColumnFormat;
      bp.props.autoLayout = true; bp.props.axisUseColumnFormat = true;
      bp.props.whiskerLineWidth = 2; bp.props.whiskerWidthRatio = 0.5;
      return {whiskerWidth, ratio1, ratio03, autoLayoutOff, axisFormatOff};
    });
    expect(result.whiskerWidth).toBe(4);
    expect(result.ratio1).toBe(1);
    expect(result.ratio03).toBe(0.3);
  });

  // #### Viewer filter formula
  await softStep('Filter formula', async () => {
    const result = await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.filter = '${AGE} > 40';
      await new Promise(r => setTimeout(r, 500));
      const filterSet = bp.props.filter;
      bp.props.filter = '';
      const filterCleared = bp.props.filter;
      return {filterSet, filterCleared};
    });
    expect(result.filterSet).toBe('${AGE} > 40');
    expect(result.filterCleared).toBe('');
  });

  // #### P-value (t-test)
  await softStep('P-value toggle', async () => {
    const result = await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.category1ColumnName = 'SEX';
      bp.props.valueColumnName = 'AGE';
      const defaultPVal = bp.props.showPValue;
      bp.props.showPValue = false;
      const off = bp.props.showPValue;
      bp.props.showPValue = true;
      const on = bp.props.showPValue;
      bp.props.category1ColumnName = 'RACE';
      await new Promise(r => setTimeout(r, 300));
      const pvalWithRace = bp.props.showPValue;
      return {defaultPVal, off, on, pvalWithRace};
    });
    expect(result.defaultPVal).toBe(true);
    expect(result.off).toBe(false);
    expect(result.on).toBe(true);
    expect(result.pvalWithRace).toBe(true);
  });

  // #### Legend
  await softStep('Legend visibility and position', async () => {
    const result = await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.category1ColumnName = 'SEX';
      await new Promise(r => setTimeout(r, 200));
      bp.props.markerColorColumnName = 'RACE';
      await new Promise(r => setTimeout(r, 200));
      const markerColor = bp.props.markerColorColumnName;
      bp.props.legendVisibility = 'Never';
      const legendNever = bp.props.legendVisibility;
      bp.props.legendVisibility = 'Always';
      const legendAlways = bp.props.legendVisibility;
      bp.props.legendPosition = 'RightTop';
      const legendRightTop = bp.props.legendPosition;
      bp.props.legendPosition = 'LeftBottom';
      const legendLeftBottom = bp.props.legendPosition;
      bp.props.markerColorColumnName = '';
      bp.props.legendVisibility = 'Auto';
      return {markerColor, legendNever, legendAlways, legendRightTop, legendLeftBottom};
    });
    expect(result.markerColor).toBe('RACE');
    expect(result.legendNever).toBe('Never');
    expect(result.legendAlways).toBe('Always');
    expect(result.legendRightTop).toBe('RightTop');
    expect(result.legendLeftBottom).toBe('LeftBottom');
  });

  // #### Layout save and restore
  await softStep('Layout save and restore', async () => {
    const layoutId = await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      bp.props.valueColumnName = 'WEIGHT';
      bp.props.category1ColumnName = 'RACE';
      await new Promise(r => setTimeout(r, 200));
      bp.props.markerColorColumnName = 'SEX';
      bp.props.showTotalCount = true;
      bp.props.plotStyle = 'violin';
      await new Promise(r => setTimeout(r, 300));
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });

    await page.evaluate(async (id: string) => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      tv.addViewer('Box plot');
      await new Promise(r => setTimeout(r, 500));
      const saved = await grok.dapi.layouts.find(id);
      tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
    }, layoutId);

    const restored = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      return {
        restored: !!bp,
        value: bp?.props?.valueColumnName,
        cat1: bp?.props?.category1ColumnName,
        markerColor: bp?.props?.markerColorColumnName,
        totalCount: bp?.props?.showTotalCount,
        plotStyle: bp?.props?.plotStyle
      };
    });
    expect(restored.restored).toBe(true);
    expect(restored.value).toBe('WEIGHT');
    expect(restored.cat1).toBe('RACE');
    expect(restored.markerColor).toBe('SEX');
    expect(restored.totalCount).toBe(true);
    expect(restored.plotStyle).toBe('violin');

    await page.evaluate(async (id: string) => {
      const saved = await grok.dapi.layouts.find(id);
      if (saved) await grok.dapi.layouts.delete(saved);
    }, layoutId);
  });

  // #### Table switching (spgi-100)
  await softStep('Table switching to spgi-100', async () => {
    // Open spgi-100, switch table on box plot, set props
    const result = await page.evaluate(async (path: string) => {
      const df2 = await grok.dapi.files.readCsv(path);
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 4000);
      });
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 3000));

      // Go back to demog view with box plot
      const views = Array.from(grok.shell.views);
      const bpView = views.find((v: any) =>
        v.type === 'TableView' && Array.from(v.viewers || []).some((vw: any) => vw.type === 'Box plot')
      );
      if (bpView) grok.shell.v = bpView;
      await new Promise(r => setTimeout(r, 500));
      const bp = grok.shell.tv.viewers.find(v => v.type === 'Box plot');
      if (!bp) return {switched: false, error: 'no box plot'};

      // Switch table using the df2 reference we still have
      bp.props.table = df2;
      await new Promise(r => setTimeout(r, 3000));
      bp.props.valueColumnName = 'Average Mass';
      await new Promise(r => setTimeout(r, 500));
      bp.props.category1ColumnName = 'Series';
      await new Promise(r => setTimeout(r, 500));
      bp.props.filter = '${Average Mass} > 225';
      await new Promise(r => setTimeout(r, 500));
      bp.props.binColorColumnName = 'TPSA';
      await new Promise(r => setTimeout(r, 500));

      return {
        switched: true,
        value: bp.props.valueColumnName,
        cat1: bp.props.category1ColumnName,
        filter: bp.props.filter,
        binColor: bp.props.binColorColumnName
      };
    }, spgiPath);
    expect(result.switched).toBe(true);
    expect(result.value).toBe('Average Mass');
    expect(result.cat1).toBe('Series');
    expect(result.filter).toBe('${Average Mass} > 225');
    expect(result.binColor).toBe('TPSA');

    await page.evaluate(() => grok.shell.closeAll());
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
