import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';

test('PC Plot tests', async ({page}) => {
  test.setTimeout(600_000);

  // The canvas-only steps below change nothing but how the plot is painted, so the
  // check is that driving them raises nothing and leaves the viewer alive.
  // grok.shell.warnings is undefined on this build, hence the page/console baseline.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  const errorCount = () => pageErrors.length + consoleErrors.length;
  const viewerAlive = () => page.evaluate(() =>
    !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot')
    && !!document.querySelector('[name="viewer-PC-Plot"]'));

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon) (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 10000});

  await softStep('Menu Ribbon and To Script', async () => {
    const viewerExists = await page.evaluate(() => !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot'));
    expect(viewerExists).toBe(true);

    await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
      }));
      await new Promise(r => setTimeout(r, 500));
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const toScript = items.find(el => el.textContent!.trim() === 'To Script');
      if (toScript) {
        const parent = toScript.closest('.d4-menu-item')!;
        parent.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        parent.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 300));
        const sub = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const js = sub.find(el => el.textContent!.trim() === 'To JavaScript');
        if (js) js.closest('.d4-menu-item')!.click();
      }
      await new Promise(r => setTimeout(r, 500));
    });

    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot');
      if (pc) pc.close();
      await new Promise(r => setTimeout(r, 500));
      const icon = document.querySelector('[name="icon-pc-plot"]');
      if (icon) (icon as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1000));
    });
    const reopened = await page.evaluate(() => !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot'));
    expect(reopened).toBe(true);
  });

  await softStep('Axis scale via the context menu', async () => {
    const menuResult = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      const openMenu = async () => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        await new Promise(r => setTimeout(r, 500));
      };
      const clickSub = async (parent: string, child: string) => {
        const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const p = items.find(el => el.textContent!.trim() === parent);
        if (!p) return false;
        const pm = p.closest('.d4-menu-item')!;
        pm.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        pm.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 300));
        const sub = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const c = sub.find(el => el.textContent!.trim() === child);
        if (c) c.closest('.d4-menu-item')!.click();
        await new Promise(r => setTimeout(r, 300));
        return !!c;
      };
      await openMenu();
      await clickSub('Y Axis', 'Global');
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const afterGlobal = pc.props.normalizeEachColumn;
      await openMenu();
      await clickSub('Y Axis', 'Normalized');
      const afterNorm = pc.props.normalizeEachColumn;
      return { afterGlobal, afterNorm };
    });
    expect(menuResult.afterGlobal).toBe(false);
    expect(menuResult.afterNorm).toBe(true);
  });

  // Which lines get painted — current, mouse-over, all — is a canvas outcome with
  // no DOM counterpart, so the check is the no-error floor plus a live viewer.
  await softStep('Selection & line display', async () => {
    const errBefore = errorCount();
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const wait = () => new Promise(r => setTimeout(r, 150));
      pc.props.showCurrentLine = false; await wait();
      pc.props.showCurrentLine = true; await wait();
      pc.props.showMouseOverLine = false; await wait();
      pc.props.showMouseOverLine = true; await wait();
      pc.props.showMouseOverRowGroup = true; await wait();
      pc.props.showAllLines = false; await wait();
      pc.props.showAllLines = true; await wait();
      pc.props.showMouseOverRowGroup = false;
    });
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  // Line widths, label orientation and margins are pure painting. The axis sliders
  // are read afterwards to confirm the layout pass rebuilt the plot.
  await softStep('Style & layout', async () => {
    const errBefore = errorCount();
    const sliders = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const wait = () => new Promise(r => setTimeout(r, 150));
      pc.props.lineWidth = 3; await wait();
      pc.props.currentLineWidth = 5; await wait();
      pc.props.mouseOverLineWidth = 5; await wait();
      pc.props.labelsOrientation = 'Vert'; await wait();
      pc.props.minMaxOrientation = 'Vert'; await wait();
      pc.props.horzMargin = 60; await wait();
      pc.props.autoLayout = false; await wait();
      pc.props.lineWidth = 0.5; pc.props.currentLineWidth = 2; pc.props.mouseOverLineWidth = 2;
      pc.props.labelsOrientation = 'Auto'; pc.props.minMaxOrientation = 'Auto';
      pc.props.horzMargin = 40; pc.props.autoLayout = true;
      await new Promise(r => setTimeout(r, 400));
      return document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]').length;
    });
    expect(sliders).toBeGreaterThan(0);
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Reset and filter visibility from the context menu', async () => {
    const result = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('dblclick', {
        bubbles: true, cancelable: true,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
      }));
      await new Promise(r => setTimeout(r, 500));
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
      }));
      await new Promise(r => setTimeout(r, 500));
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const rv = items.find(el => el.textContent!.trim() === 'Reset View');
      if (rv) rv.closest('.d4-menu-item')!.click();
      await new Promise(r => setTimeout(r, 300));
      return { resetFound: !!rv };
    });
    expect(result.resetFound).toBe(true);
  });

  await softStep('Filter panel interaction', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 500));
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 50});
      await new Promise(r => setTimeout(r, 500));
      const afterFilter = grok.shell.tv.dataFrame.filter.trueCount;
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 18, max: 89});
      await new Promise(r => setTimeout(r, 500));
      const afterReset = grok.shell.tv.dataFrame.filter.trueCount;
      return { afterFilter, afterReset };
    });
    expect(result.afterFilter).toBeLessThan(5850);
    expect(result.afterReset).toBe(5850);
  });

  // Each axis carries its own range slider named after its column, and the sliders
  // sit in painted order, so their sequence is the rendered axis order.
  await softStep('Column reordering from the Context Panel list', async () => {
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const sliderOrder = () =>
        Array.from(document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]'))
          .map((e) => e.getAttribute('name')!.replace('axis-slider-', ''));
      pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      await new Promise(r => setTimeout(r, 800));
      const before = sliderOrder();
      pc.props.columnNames = ['WEIGHT', 'AGE', 'HEIGHT'];
      await new Promise(r => setTimeout(r, 800));
      const after = sliderOrder();
      pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      await new Promise(r => setTimeout(r, 500));
      return { before, after };
    });
    expect(result.before).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
    expect(result.after).toEqual(['WEIGHT', 'AGE', 'HEIGHT']);
  });

  // Every box-plot component toggle draws to canvas, so this is the no-error floor
  // over the whole component surface.
  await softStep('Density component toggles', async () => {
    const errBefore = errorCount();
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const wait = () => new Promise(r => setTimeout(r, 150));
      pc.props.showDensity = true;
      pc.props.densityStyle = 'box plot'; await wait();
      pc.props.showInterquartileRange = false; await wait();
      pc.props.showInterquartileRange = true;
      pc.props.showUpperDash = false; pc.props.showUpperDash = true; await wait();
      pc.props.showLowerDash = false; pc.props.showLowerDash = true; await wait();
      pc.props.showMeanCross = false; pc.props.showMeanCross = true; await wait();
      pc.props.showMedian = false; pc.props.showMedian = true; await wait();
      pc.props.showCircles = true; await wait();
      pc.props.densityStyle = 'violin plot'; await wait();
      pc.props.bins = 200; await wait();
      pc.props.whiskerLineWidth = 5; await wait();
      pc.props.densityStyle = 'circles';
      pc.props.bins = 100; pc.props.whiskerLineWidth = 2;
      pc.props.showDensity = false;
      await new Promise(r => setTimeout(r, 300));
    });
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  // Legend visibility is a real DOM signal; the gradient options (log axis,
  // inversion, min/max clamps) are canvas-only and are driven under the floor.
  await softStep('Color coding, legend & grid coloring', async () => {
    const errBefore = errorCount();
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const wait = (ms = 600) => new Promise(r => setTimeout(r, ms));
      const legend = () => {
        const el = document.querySelector('[name="viewer-PC-Plot"] .d4-legend') as HTMLElement | null;
        return {present: !!el, labels: el ? (el.innerText || '').split('\n').map(s => s.trim()).filter(Boolean) : []};
      };

      // Numeric colouring: gradient options are canvas-only, just drive them.
      pc.props.colorColumnName = 'AGE'; await wait();
      pc.props.colorAxisType = 'logarithmic'; await wait(300);
      pc.props.invertColorScheme = true; await wait(300);
      pc.props.invertColorScheme = false;
      pc.props.colorMin = 30; pc.props.colorMax = 60; await wait(300);
      pc.props.colorMin = null; pc.props.colorMax = null; pc.props.colorAxisType = 'linear';

      // Categorical colouring: the legend lists the column's categories.
      pc.props.colorColumnName = 'RACE'; await wait();
      const categorical = legend();
      pc.props.legendPosition = 'Left'; await wait(300);
      pc.props.legendPosition = 'Right'; pc.props.legendPosition = 'Top';
      pc.props.legendPosition = 'Bottom'; await wait(300);
      pc.props.legendVisibility = 'Never'; await wait();
      const hidden = legend();
      pc.props.legendVisibility = 'Auto'; await wait();
      const restored = legend();

      // Colour coding set on the grid column, read by the plot.
      pc.props.colorColumnName = 'HEIGHT';
      const df = grok.shell.tv.dataFrame;
      df.col('HEIGHT').meta.colors.setLinear([DG.Color.blue, DG.Color.red]);
      await wait(500);
      df.col('HEIGHT').meta.colors.setConditional({'20-150': DG.Color.green, '150-250': DG.Color.orange});
      await wait(500);
      df.col('HEIGHT').meta.colors.setLinear();
      pc.props.colorColumnName = '';
      pc.props.legendPosition = 'Auto';
      await wait(300);
      return {categorical, hidden, restored};
    });
    expect(result.categorical.present).toBe(true);
    expect(result.hidden.present).toBe(false);
    expect(result.restored.labels).toEqual(result.categorical.labels);
    expect(errorCount()).toBe(errBefore);
  });

  // The description is rendered inside the viewer element and can be read back; the
  // title lives in the surrounding header chrome, so it is driven but not asserted.
  await softStep('Title and description', async () => {
    const errBefore = errorCount();
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const wait = () => new Promise(r => setTimeout(r, 600));
      const shownText = () =>
        ((document.querySelector('[name="viewer-PC-Plot"]') as HTMLElement).innerText || '')
          .replace(/\s+/g, ' ').trim();
      pc.props.title = 'My PC Plot'; await wait();
      pc.props.description = 'Test description'; await wait();
      const withDescription = shownText();
      pc.props.descriptionPosition = 'Bottom'; await wait();
      const moved = shownText();
      pc.props.title = ''; pc.props.description = ''; await wait();
      const cleared = shownText();
      return {withDescription, moved, cleared};
    });
    expect(result.withDescription).toContain('Test description');
    expect(result.moved).toContain('Test description');
    expect(result.cleared).not.toContain('Test description');
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Pick Up / Apply', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const pc1 = tv.viewers.find(v => v.type === 'PC Plot')!;
      pc1.props.columnNames = ['AGE', 'WEIGHT', 'STARTED'];
      pc1.props.logColumnsColumnNames = ['AGE'];
      pc1.props.colorColumnName = 'RACE';
      pc1.props.legendPosition = 'Left';
      pc1.props.title = 'Source Plot';
      tv.addViewer('PC Plot');
      await new Promise(r => setTimeout(r, 500));
      const clickSub = async (idx: number, parent: string, child: string) => {
        const viewers = document.querySelectorAll('[name="viewer-PC-Plot"]');
        const canvas = viewers[idx].querySelector('canvas[name="canvas"]')!;
        const rect = canvas.getBoundingClientRect();
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        await new Promise(r => setTimeout(r, 500));
        const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const p = items.find(el => el.textContent!.trim() === parent);
        if (!p) return;
        const pm = p.closest('.d4-menu-item')!;
        pm.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        pm.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 300));
        const sub = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const c = sub.find(el => el.textContent!.trim() === child);
        if (c) c.closest('.d4-menu-item')!.click();
        await new Promise(r => setTimeout(r, 500));
      };
      await clickSub(0, 'Pick Up / Apply', 'Pick Up');
      await clickSub(1, 'Pick Up / Apply', 'Apply');
      const pcs = tv.viewers.filter(v => v.type === 'PC Plot');
      const r = { pc2Cols: pcs[1]?.props.columnNames?.slice(), pc2Color: pcs[1]?.props.colorColumnName, pc2Title: pcs[1]?.props.title };
      if (pcs[1]) pcs[1].close();
      const pcF = tv.viewers.find(v => v.type === 'PC Plot')!;
      pcF.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      pcF.props.logColumnsColumnNames = []; pcF.props.colorColumnName = '';
      pcF.props.legendPosition = 'Auto'; pcF.props.title = '';
      return r;
    });
    expect(result.pc2Cols).toEqual(['AGE', 'WEIGHT', 'STARTED']);
    expect(result.pc2Color).toBe('RACE');
    expect(result.pc2Title).toBe('Source Plot');
  });

  await softStep('Layout and project save/restore', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));
      tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 500));
      const saved = await grok.dapi.layouts.find(layoutId);
      tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
      const hasScatter = tv.viewers.some(v => v.type === 'Scatter plot');
      const hasPc = tv.viewers.some(v => v.type === 'PC Plot');
      await grok.dapi.layouts.delete(saved);
      return { hasScatter, hasPc };
    });
    expect(result.hasScatter).toBe(false);
    expect(result.hasPc).toBe(true);
  });

  await softStep('Table switching and transformation', async () => {
    const result = await page.evaluate(async (path) => {
      const df2 = await grok.dapi.files.readCsv(path);
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
      const views = Array.from(grok.shell.views);
      const demogView = views.find((v: any) => v.name === 'Table');
      if (demogView) grok.shell.v = demogView;
      await new Promise(r => setTimeout(r, 500));
      const pc = grok.shell.tv.addViewer('PC Plot');
      await new Promise(r => setTimeout(r, 500));
      pc.props.table = df2.name;
      await new Promise(r => setTimeout(r, 500));
      const tableSet = pc.dataFrame?.name;
      const sliderAxes = () =>
        Array.from(document.querySelectorAll('[name="viewer-PC-Plot"] [name^="axis-slider-"]'))
          .map((e) => e.getAttribute('name')!.replace('axis-slider-', ''));
      await new Promise(r => setTimeout(r, 1500));
      const axesBefore = sliderAxes();
      pc.props.transformation = '[{"#type":"GroupAggregation","aggType":"key","colName":"Chemist 521"},{"#type":"GroupAggregation","aggType":"pivot","colName":"Series"},{"#type":"GroupAggregation","aggType":"count","colName":"Id"}]';
      await new Promise(r => setTimeout(r, 3000));
      // The pivot replaces the axes with one generated column per Series value, so
      // the slider names show whether the aggregation was applied.
      const axesAfter = sliderAxes();
      pc.props.transformation = '';
      await new Promise(r => setTimeout(r, 1500));
      const axesReverted = sliderAxes();
      pc.close();
      return { spgiRows: df2.rowCount, tableSet, axesBefore, axesAfter, axesReverted };
    }, spgiPath);
    expect(result.spgiRows).toBe(100);
    expect(result.tableSet).toBeTruthy();
    // Pivoted axes are the Series categories, not the raw numeric columns.
    expect(result.axesAfter).not.toEqual(result.axesBefore);
    expect(result.axesAfter).toContain('Triazoles');
    expect(result.axesReverted).toEqual(result.axesBefore);
  });

  v.finishSpec();
});
