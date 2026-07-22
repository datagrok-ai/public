import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Histogram tests', async ({page}) => {
  test.setTimeout(600_000);

  // Several sections below drive canvas-only properties whose only observable
  // outcome is how the plot is painted. grok.shell.warnings is undefined on this
  // build, so those steps assert a no-error floor built from page/console errors.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  const errorCount = () => pageErrors.length + consoleErrors.length;
  const viewerAlive = () => page.evaluate(() =>
    !!grok.shell.tv.viewers.find(v => v.type === 'Histogram')
    && !!document.querySelector('[name="viewer-Histogram"]'));

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'histogram', 'Histogram');

  // Spline and fill-spline are pure canvas painting with no DOM or dataframe
  // counterpart, so the check is the no-error floor plus a live viewer.
  await softStep('Spline mode', async () => {
    const errBefore = errorCount();
    await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      const wait = () => new Promise(r => setTimeout(r, 200));
      h.props.spline = true; await wait();
      h.props.fillSpline = true; await wait();
      h.props.fillSpline = false; await wait();
      h.props.spline = false; await wait();
    });
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  // Two selector toggles carry a real DOM signal: showColumnSelector controls the
  // value-column combobox ([name="div-column-combobox-value"]) and
  // showSplitSelector controls the split-column combobox
  // ([name="div-column-combobox-split"]) — each combobox's computed visibility
  // flips with its prop. Assert that visibility round-trip.
  await softStep('Appearance — selector visibility', async () => {
    const r = await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      const wait = (ms = 500) => new Promise(res => setTimeout(res, ms));
      const vis = (name: string) => {
        const el = (h.root as HTMLElement).querySelector(`[name="${name}"]`) as HTMLElement | null;
        if (!el) return -1;
        const cs = getComputedStyle(el);
        return (el.offsetParent !== null && cs.display !== 'none' && cs.visibility !== 'hidden') ? 1 : 0;
      };
      const read = () => ({value: vis('div-column-combobox-value'), split: vis('div-column-combobox-split')});
      // Warm up: set the split column and let the combobox settle out of its
      // cold-start layout transient, then assert the OFF→ON toggle round-trip
      // (the settled state is deterministic, the cold-start state is not).
      h.props.splitColumnName = 'SEX';
      h.props.showColumnSelector = true; h.props.showSplitSelector = true; await wait(700);
      h.props.showColumnSelector = false; h.props.showSplitSelector = false; await wait();
      const off = read();
      h.props.showColumnSelector = true; h.props.showSplitSelector = true; await wait();
      const on = read();
      h.props.splitColumnName = '';
      return {off, on};
    });
    expect(r.off.value).toBe(0);
    expect(r.off.split).toBe(0);
    expect(r.on.value).toBe(1);
    expect(r.on.split).toBe(1);
  });

  // The axis toggles, X-axis height, allowColumnSelection, showBinSelector and
  // showRangeSlider have no headless DOM signal — axes and the bin/range strip
  // are canvas-drawn and the bin/range controls do not add, remove, or change the
  // visibility of a DOM node. Drive every one and revert to defaults under the
  // no-error floor.
  await softStep('Appearance — canvas floor', async () => {
    const errBefore = errorCount();
    await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      const wait = () => new Promise(r => setTimeout(r, 150));
      h.props.showXAxis = true; h.props.showYAxis = true; await wait();
      h.props.showXAxis = false; h.props.showYAxis = false; await wait();
      h.props.xAxisHeight = 30; await wait();
      h.props.allowColumnSelection = false; await wait();
      h.props.showBinSelector = false; await wait();
      h.props.showRangeSlider = false; await wait();
      h.props.showXAxis = true; h.props.showYAxis = true; h.props.allowColumnSelection = true;
      h.props.showBinSelector = true; h.props.showRangeSlider = true;
      h.props.xAxisHeight = 20;
      await new Promise(r => setTimeout(r, 300));
    });
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  // Splitting by a categorical column renders a real legend (one item per
  // category) — the DOM signal. The description is drawn inside the viewer and is
  // read back from its innerText. Legend position and the title (header chrome)
  // are driven under the floor.
  await softStep('Labels', async () => {
    const errBefore = errorCount();
    const result = await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      const wait = () => new Promise(r => setTimeout(r, 600));
      const legend = () => {
        const items = Array.from(h.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
        return {
          rendered: !!h.root.querySelector('[name="legend"]'),
          labels: items.map(it => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim()),
        };
      };
      h.props.splitColumnName = 'SEX';
      h.props.legendVisibility = 'Always'; await wait();
      const sexLegend = legend();
      h.props.splitColumnName = 'RACE'; await wait();
      const raceLegend = legend();
      h.props.legendVisibility = 'Never'; await wait();
      const hiddenLegend = legend();
      h.props.legendVisibility = 'Always';
      h.props.legendPosition = 'RightTop'; await wait();
      const movedLegend = legend();
      h.props.splitColumnName = ''; h.props.legendVisibility = 'Auto'; await wait();

      const shownText = () =>
        ((h.root as HTMLElement).innerText || '').replace(/\s+/g, ' ').trim();
      h.props.showTitle = true;
      h.props.title = 'Age Distribution';
      h.props.description = 'Shows distribution of patient ages';
      h.props.descriptionVisibilityMode = 'Always';
      h.props.descriptionPosition = 'Bottom'; await wait();
      const withDescription = shownText();
      h.props.description = ''; h.props.title = ''; h.props.showTitle = false; await wait();
      const cleared = shownText();
      return {sexLegend, raceLegend, hiddenLegend, movedLegend, withDescription, cleared};
    });
    expect(result.sexLegend.rendered).toBe(true);
    expect(result.sexLegend.labels.slice().sort()).toEqual(['F', 'M']);
    expect(result.raceLegend.labels.length).toBe(4);
    expect(result.hiddenLegend.rendered).toBe(false);
    expect(result.movedLegend.rendered).toBe(true);
    expect(result.withDescription).toContain('Shows distribution of patient ages');
    expect(result.cleared).not.toContain('Shows distribution of patient ages');
    expect(errorCount()).toBe(errBefore);
  });

  // Right-clicking the canvas opens a menu that DOES include histogram-specific
  // entries. Assert the promised items are present, then toggle Show Filtered Out
  // Rows from the menu and confirm the prop round-trips. The X-axis-specific
  // right-click from the manual scenario is reduced to the same items surfaced in
  // the main-area menu (Actuation note in histogram.md).
  await softStep('Context menu', async () => {
    const result = await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      const canvas = document.querySelector('[name="viewer-Histogram"] canvas')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      await new Promise(r => setTimeout(r, 600));
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label')).map(e => e.textContent!.trim());
      const before = h.props.showFilteredOutRows;
      const sfor = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(e => e.textContent!.trim() === 'Show Filtered Out Rows');
      if (sfor) (sfor.closest('.d4-menu-item') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 500));
      const after = h.props.showFilteredOutRows;
      h.props.showFilteredOutRows = before;
      return {items, before, after};
    });
    for (const label of ['Show Filtered Out Rows', 'Selection', 'Show Current Row',
      'Show Mouse Over Row', 'Show Mouse Over Row Group', 'Show X Axis', 'Axis Font', 'Controls Font'])
      expect(result.items).toContain(label);
    expect(result.after).toBe(!result.before);
    await page.keyboard.press('Escape');
  });

  await softStep('Layout persistence', async () => {
    const result = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Histogram') as any;
      const tv = grok.shell.tv;

      h.props.valueColumnName = 'WEIGHT';
      h.props.bins = 15;
      h.props.splitColumnName = 'RACE';
      h.props.splitStack = true;
      await new Promise(res => setTimeout(res, 500));

      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(res => setTimeout(res, 1000));

      h.close();
      await new Promise(res => setTimeout(res, 500));

      const saved = await grok.dapi.layouts.find(layoutId);
      tv.loadLayout(saved);
      await new Promise(res => setTimeout(res, 3000));

      const h2 = Array.from(tv.viewers).find((v: any) => v.type === 'Histogram') as any;
      const r: any[] = [];
      r.push(h2 != null);
      r.push(h2 ? h2.props.valueColumnName : 'NOT_RESTORED');
      r.push(h2 ? h2.props.bins : 'NOT_RESTORED');
      r.push(h2 ? h2.props.splitColumnName : 'NOT_RESTORED');
      r.push(h2 ? h2.props.splitStack : 'NOT_RESTORED');

      await grok.dapi.layouts.delete(saved);

      return r;
    });
    expect(result[0]).toBe(true);
    expect(result[1]).toBe('WEIGHT');
    expect(result[2]).toBe(15);
    expect(result[3]).toBe('RACE');
    expect(result[4]).toBe(true);
  });

  // Row Source = Selected paints only the selected rows, so the histogram's
  // canvas content shrinks against the All baseline and restores exactly when
  // switched back. SPGI is used because it carries a selection large enough to
  // measure (Actuation note in histogram.md).
  await softStep('Data — row source', async () => {
    const errBefore = errorCount();
    await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      df.name = 'SPGI';
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });
      (document.querySelector('[name="icon-histogram"]') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1200));
      for (let i = 0; i < 5; i++) df.selection.set(i, true);
      await new Promise(r => setTimeout(r, 400));
      const h = tv.viewers.find(v => v.type === 'Histogram')!;
      h.props.rowSource = 'All';
      await new Promise(r => setTimeout(r, 600));
    });
    const allPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    const allHue = await v.countSelectionHuePixels(page, 'Histogram');
    await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      h.props.rowSource = 'Selected';
      await new Promise(r => setTimeout(r, 600));
    });
    const selectedPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    const selectedHue = await v.countSelectionHuePixels(page, 'Histogram');
    await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      h.props.rowSource = 'All';
      await new Promise(r => setTimeout(r, 500));
    });
    const restoredPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    const restoredHue = await v.countSelectionHuePixels(page, 'Histogram');
    // The orange selection overlay is a second, independent signal: present under
    // All, exactly 0 under Selected (no overlay when the whole chart IS the
    // selection), and restored under All again.
    expect(allPx).toBeGreaterThan(1000);
    expect(allPx - selectedPx).toBeGreaterThan(1000);
    expect(Math.abs(restoredPx - allPx)).toBeLessThan(500);
    expect(allHue).toBeGreaterThan(100);
    expect(selectedHue).toBe(0);
    expect(restoredHue).toBeGreaterThan(100);
    expect(errorCount()).toBe(errBefore);
  });

  // The viewer `filter` formula narrows which rows the histogram paints (it does
  // not touch df.filter.trueCount), so the signal is the canvas content shrinking
  // under the predicate and restoring when cleared. demog is used because its AGE
  // column carries the filter expression; SPGI has no AGE column (Actuation note
  // in histogram.md).
  await softStep('Data — filter formula', async () => {
    const errBefore = errorCount();
    await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });
      (document.querySelector('[name="icon-histogram"]') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1200));
      const h = tv.viewers.find(v => v.type === 'Histogram')!;
      h.props.valueColumnName = 'AGE';
      h.props.filter = '';
      await new Promise(r => setTimeout(r, 600));
    });
    const clearedPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      h.props.filter = '${AGE} > 40';
      await new Promise(r => setTimeout(r, 700));
    });
    const filteredPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      h.props.filter = '';
      await new Promise(r => setTimeout(r, 600));
    });
    const restoredPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    expect(clearedPx).toBeGreaterThan(1000);
    expect(clearedPx - filteredPx).toBeGreaterThan(400);
    expect(Math.abs(restoredPx - clearedPx)).toBeLessThan(300);
    expect(errorCount()).toBe(errBefore);
  });

  // Rebinding the `table` prop to another open table swaps the available columns,
  // so the auto-picked value column changes and the bound dataFrame follows.
  await softStep('Data — table switching', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      const dfSpgi = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      dfSpgi.name = 'SPGI';
      grok.shell.addTableView(dfSpgi);
      await new Promise(resolve => {
        const sub = dfSpgi.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });
      (document.querySelector('[name="icon-histogram"]') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1200));

      const dfDemog = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      dfDemog.name = 'demog';
      grok.shell.addTableView(dfDemog);
      await new Promise(r => setTimeout(r, 800));

      // Return to the SPGI view that hosts the histogram.
      const spgiView = Array.from(grok.shell.tableViews).find((t: any) => t.dataFrame?.name === 'SPGI');
      if (spgiView) grok.shell.v = spgiView as any;
      await new Promise(r => setTimeout(r, 500));

      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      const before = {dfName: h.dataFrame?.name, value: h.props.valueColumnName};
      h.props.table = 'demog';
      await new Promise(r => setTimeout(r, 1000));
      const after = {dfName: h.dataFrame?.name, value: h.props.valueColumnName};
      const demogHasValueCol = dfDemog.columns.names().includes(after.value);
      return {before, after, demogHasValueCol};
    });
    expect(result.before.dfName).toBe('SPGI');
    expect(result.after.dfName).toBe('demog');
    expect(result.after.value).not.toBe(result.before.value);
    expect(result.demogHasValueCol).toBe(true);
  });

  v.finishSpec();
});
