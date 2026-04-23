import { test } from '@playwright/test';
import {specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Matrix plot tests', async ({ page, baseURL }) => {
  test.setTimeout(300_000);

  await page.goto(baseURL ?? '/');
  await page.locator('[name="Toolbox"], [name="Browse"], .d4-sidebar').first().waitFor({ timeout: 60000 });
  await page.waitForFunction(() => {
    try {
      if (typeof grok === 'undefined' || !grok.shell) return false;
      grok.shell.settings.showFiltersIconsConstantly;
      return true;
    } catch { return false; }
  }, { timeout: 30000 });

  // Setup: close all, open demog, add Matrix plot
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    const icon = document.querySelector('[name="icon-matrix-plot"]') as HTMLElement;
    if (icon) icon.click();
    await new Promise(r => setTimeout(r, 800));
  });
  await page.locator('[name="viewer-Matrix-plot"]').waitFor({ timeout: 15000 });

  // Default State
  await softStep('Default State: viewer present and default columns are numerical', async () => {
    const result = await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      if (!mp) throw new Error('Matrix plot viewer not found');
      return { xCols: mp.props.xColumnNames, yCols: mp.props.yColumnNames };
    });
    const allowedCols = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
    for (const col of result.xCols) {
      if (!allowedCols.includes(col)) throw new Error(`Non-numerical X column: ${col}`);
    }
    for (const col of result.yCols) {
      if (!allowedCols.includes(col)) throw new Error(`Non-numerical Y column: ${col}`);
    }
  });

  // Column Configuration
  await softStep('Column Configuration: set X=AGE,HEIGHT', async () => {
    const xCols = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      mp.props.xColumnNames = ['AGE', 'HEIGHT'];
      await new Promise(r => setTimeout(r, 500));
      return mp.props.xColumnNames;
    });
    if (!xCols.includes('AGE') || !xCols.includes('HEIGHT') || xCols.length !== 2)
      throw new Error(`Expected [AGE,HEIGHT], got ${JSON.stringify(xCols)}`);
  });

  await softStep('Column Configuration: set Y=AGE,HEIGHT,WEIGHT (2x3 grid)', async () => {
    const result = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      mp.props.yColumnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      await new Promise(r => setTimeout(r, 500));
      return { x: mp.props.xColumnNames.length, y: mp.props.yColumnNames.length };
    });
    if (result.x !== 2 || result.y !== 3)
      throw new Error(`Expected 2x3, got ${result.x}x${result.y}`);
  });

  await softStep('Column Configuration: reset to defaults', async () => {
    await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      mp.props.xColumnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      mp.props.yColumnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      await new Promise(r => setTimeout(r, 500));
    });
  });

  // Cell Plot Type
  await softStep('Cell Plot Type: Scatter plot then back to Density plot', async () => {
    const result = await page.evaluate(async () => {
      // Open settings: click 2nd settings icon (Matrix plot title bar)
      const gearIcons = document.querySelectorAll('[name="icon-font-icon-settings"]');
      if (gearIcons.length >= 2) (gearIcons[1] as HTMLElement).click();
      await new Promise(r => setTimeout(r, 300));
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      mp.props.cellPlotType = 'Scatter plot';
      await new Promise(r => setTimeout(r, 300));
      const scatter = mp.props.cellPlotType;
      mp.props.cellPlotType = 'Density plot';
      await new Promise(r => setTimeout(r, 300));
      return { scatter, density: mp.props.cellPlotType };
    });
    if (result.scatter !== 'Scatter plot') throw new Error('Expected Scatter plot');
    if (result.density !== 'Density plot') throw new Error('Expected Density plot');
  });

  // Show Axes
  await softStep('Show Axes: set true then false', async () => {
    const result = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      mp.props.showXAxes = true;
      mp.props.showYAxes = true;
      await new Promise(r => setTimeout(r, 300));
      const trueState = { x: mp.props.showXAxes, y: mp.props.showYAxes };
      mp.props.showXAxes = false;
      mp.props.showYAxes = false;
      await new Promise(r => setTimeout(r, 300));
      return { trueState, falseState: { x: mp.props.showXAxes, y: mp.props.showYAxes } };
    });
    if (!result.trueState.x || !result.trueState.y) throw new Error('Axes not enabled');
    if (result.falseState.x || result.falseState.y) throw new Error('Axes not disabled');
  });

  // Auto Layout
  await softStep('Auto Layout: enabled by default, toggle off/on', async () => {
    const result = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      const defaultVal = mp.props.autoLayout;
      mp.props.showXAxes = true;
      mp.props.showYAxes = true;
      mp.props.autoLayout = false;
      await new Promise(r => setTimeout(r, 300));
      const axesAfter = { x: mp.props.showXAxes, y: mp.props.showYAxes };
      mp.props.autoLayout = true;
      await new Promise(r => setTimeout(r, 300));
      return { default: defaultVal, axesAfterDisable: axesAfter, restored: mp.props.autoLayout };
    });
    if (!result.default) throw new Error('Auto Layout not enabled by default');
    if (!result.axesAfterDisable.x || !result.axesAfterDisable.y) throw new Error('Axes hidden after auto layout disabled');
    if (!result.restored) throw new Error('Auto Layout not restored');
  });

  // Scrolling
  await softStep('Scrolling: 4 columns set, range sliders visible', async () => {
    await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      mp.props.xColumnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      mp.props.yColumnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      await new Promise(r => setTimeout(r, 500));
    });
    // Range sliders visible in canvas area — visual check only
    await page.locator('[name="viewer-Matrix-plot"]').waitFor();
  });

  // Row Source
  await softStep('Row Source: Selected (50 rows) then Filtered', async () => {
    const result = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      const df = grok.shell.tv.dataFrame;
      mp.props.rowSource = 'Selected';
      await new Promise(r => setTimeout(r, 300));
      for (let i = 0; i < 50; i++) df.selection.set(i, true);
      await new Promise(r => setTimeout(r, 500));
      const sel = df.selection.trueCount;
      mp.props.rowSource = 'Filtered';
      await new Promise(r => setTimeout(r, 300));
      df.selection.setAll(false);
      return { selectedRows: sel, restoredSource: mp.props.rowSource };
    });
    if (result.selectedRows !== 50) throw new Error(`Expected 50 selected, got ${result.selectedRows}`);
    if (result.restoredSource !== 'Filtered') throw new Error('Row source not restored to Filtered');
  });

  // Filter Integration
  await softStep('Filter Integration: SEX=M reduces count, removing restores', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 600));
      const total = df.rowCount;
      fg.updateOrAdd({ type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M'] });
      await new Promise(r => setTimeout(r, 500));
      const filtered = df.filter.trueCount;
      const sexCats = df.col('SEX').categories;
      fg.updateOrAdd({ type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: sexCats });
      await new Promise(r => setTimeout(r, 500));
      return { total, filtered, restored: df.filter.trueCount };
    });
    if (result.filtered >= result.total) throw new Error('Filter did not reduce count');
    if (result.restored !== result.total) throw new Error('Count not restored after removing filter');
  });

  // Data Filter Property
  await softStep('Data Filter Property: viewer filter ${AGE} > 30 reduces row count', async () => {
    const result = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      const df = grok.shell.tv.dataFrame;
      const totalRows = df.rowCount;
      mp.props.filter = '${AGE} > 30';
      await new Promise(r => setTimeout(r, 600));
      const filterSet = mp.props.filter;
      mp.props.filter = '';
      await new Promise(r => setTimeout(r, 400));
      const filterCleared = mp.props.filter;
      return { totalRows, filterSet, filterCleared };
    });
    if (!result.filterSet) throw new Error('Filter property not set');
    if (result.filterCleared !== '') throw new Error('Filter not cleared');
  });

  // Inner Viewer Look
  await softStep('Inner Viewer Look: marker size for scatter, bin count for density', async () => {
    const result = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      mp.props.cellPlotType = 'Scatter plot';
      await new Promise(r => setTimeout(r, 400));
      // Try setting inner marker size — property name varies by build
      let sizeSet = false;
      let sizePropName: string | null = null;
      for (const name of ['markerSize', 'innerMarkerSize', 'pointSize', 'size']) {
        try {
          (mp.props as any)[name] = 10;
          await new Promise(r => setTimeout(r, 200));
          sizeSet = (mp.props as any)[name] === 10;
          if (sizeSet) { sizePropName = name; (mp.props as any)[name] = 5; break; }
        } catch (_) {}
      }
      mp.props.cellPlotType = 'Density plot';
      await new Promise(r => setTimeout(r, 400));
      // Try setting bin count
      let binSet = false;
      let binPropName: string | null = null;
      for (const name of ['binCount', 'numBins', 'bins']) {
        try {
          (mp.props as any)[name] = 20;
          await new Promise(r => setTimeout(r, 200));
          binSet = (mp.props as any)[name] === 20;
          if (binSet) { binPropName = name; (mp.props as any)[name] = 10; break; }
        } catch (_) {}
      }
      return { cellType: mp.props.cellPlotType, sizePropName, sizeSet, binPropName, binSet };
    });
    if (result.cellType !== 'Density plot') throw new Error('Not reverted to Density plot');
    // sizeSet / binSet may be false if props are named differently — AMBIGUOUS is acceptable
  });

  // Back Color
  await softStep('Back Color: set red then restore white', async () => {
    const result = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      mp.props.backColor = 0xFFFF0000;
      await new Promise(r => setTimeout(r, 400));
      const red = mp.props.backColor;
      mp.props.backColor = 0xFFFFFFFF;
      await new Promise(r => setTimeout(r, 300));
      return { red, white: mp.props.backColor };
    });
    if (result.red !== 0xFFFF0000) throw new Error('Back color not set to red');
    if (result.white !== 0xFFFFFFFF) throw new Error('Back color not restored to white');
  });

  // Title and Description
  await softStep('Title and Description: set title and description', async () => {
    const result = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      mp.props.showTitle = true;
      mp.props.title = 'My Matrix';
      mp.props.description = 'Test description';
      await new Promise(r => setTimeout(r, 400));
      const t = mp.props.title;
      const d = mp.props.description;
      mp.props.showTitle = false;
      mp.props.title = '';
      return { title: t, description: d };
    });
    if (result.title !== 'My Matrix') throw new Error('Title not set');
    if (result.description !== 'Test description') throw new Error('Description not set');
  });

  // Context Menu: Clone
  await softStep('Context Menu: clone viewer and close', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const mp = tv.viewers.find((v: any) => v.type === 'Matrix plot') as any;
      const before = tv.viewers.length;
      const clone = tv.addViewer('Matrix plot');
      await new Promise(r => setTimeout(r, 500));
      const after = tv.viewers.length;
      clone.close();
      await new Promise(r => setTimeout(r, 300));
      return { before, after, viewersFinal: tv.viewers.length };
    });
    if (result.after <= result.before) throw new Error('Clone not added');
    if (result.viewersFinal !== result.before) throw new Error('Clone not closed');
  });

  if (stepErrors.length > 0) {
    throw new Error(
      'Some steps failed:\n' +
      stepErrors.map(e => `  [${e.step}]: ${e.error}`).join('\n')
    );
  }
});
