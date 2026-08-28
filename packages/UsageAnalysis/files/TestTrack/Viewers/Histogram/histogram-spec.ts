/* ---
realizes: []
--- */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Histogram tests', async ({page}) => {
  test.setTimeout(600_000);

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

  await softStep('Spline mode', async () => {
    const errBefore = errorCount();
    await v.setViewerProps(page, 'Histogram', [
      {set: {spline: true}},
      {set: {fillSpline: true}},
      {set: {fillSpline: false}},
      {set: {spline: false}},
    ], 200);
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Appearance — selector visibility', async () => {
    const r = await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      const vis = (name: string) => {
        const el = (h.root as HTMLElement).querySelector(`[name="${name}"]`) as HTMLElement | null;
        if (!el) return -1;
        const cs = getComputedStyle(el);
        return (el.offsetParent !== null && cs.display !== 'none' && cs.visibility !== 'hidden') ? 1 : 0;
      };
      const read = () => ({value: vis('div-column-combobox-value'), split: vis('div-column-combobox-split')});

      const settle = async (want: number, cap: number) => {
        const t0 = Date.now();
        let hits = 0;
        while (Date.now() - t0 < cap) {
          const s = read();
          hits = (s.value === want && s.split === want) ? hits + 1 : 0;
          if (hits >= 2) return;
          await new Promise(res => setTimeout(res, 25));
        }
      };

      h.props.splitColumnName = 'SEX';
      h.props.showColumnSelector = true; h.props.showSplitSelector = true;
      await settle(1, 700);
      h.props.showColumnSelector = false; h.props.showSplitSelector = false;
      await settle(0, 500);
      const off = read();
      h.props.showColumnSelector = true; h.props.showSplitSelector = true;
      await settle(1, 500);
      const on = read();
      h.props.splitColumnName = '';
      return {off, on};
    });
    expect(r.off.value).toBe(0);
    expect(r.off.split).toBe(0);
    expect(r.on.value).toBe(1);
    expect(r.on.split).toBe(1);
  });

  await softStep('Appearance — canvas floor', async () => {
    const errBefore = errorCount();
    await v.setViewerProps(page, 'Histogram', [
      {set: {showXAxis: true, showYAxis: true}},
      {set: {showXAxis: false, showYAxis: false}},
      {set: {xAxisHeight: 30}},
      {set: {allowColumnSelection: false}},
      {set: {showBinSelector: false}},
      {set: {showRangeSlider: false}},
      {set: {showXAxis: true, showYAxis: true, allowColumnSelection: true,
        showBinSelector: true, showRangeSlider: true, xAxisHeight: 20}, wait: 300},
    ], 150);
    expect(await viewerAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Labels', async () => {
    const errBefore = errorCount();
    const result = await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      const legend = () => {
        // a hidden legend keeps its element in the DOM with display:none (LegendHost.apply,
        // legend_host.dart) — presence must check the computed style, not existence
        const host = Array.from(h.root.querySelectorAll('[name="legend"]'))
          .find((e: any) => getComputedStyle(e).display !== 'none') as HTMLElement | undefined;
        const items = host
          ? Array.from(host.querySelectorAll('.d4-legend-item')) as HTMLElement[]
          : [];
        return {
          rendered: !!host,
          labels: items.map(it => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim()),
        };
      };
      const until = async (ok: () => boolean, cap: number) => {
        const t0 = Date.now();
        while (!ok() && Date.now() - t0 < cap)
          await new Promise(r => setTimeout(r, 25));
      };
      h.props.splitColumnName = 'SEX';
      h.props.legendVisibility = 'Always';
      await until(() => legend().rendered && legend().labels.length > 0, 600);
      const sexLegend = legend();
      h.props.splitColumnName = 'RACE';
      await until(() => legend().labels.join('|') !== sexLegend.labels.join('|'), 600);
      const raceLegend = legend();
      h.props.legendVisibility = 'Never';
      await until(() => !legend().rendered, 600);
      const hiddenLegend = legend();
      h.props.legendVisibility = 'Always';
      h.props.legendPosition = 'RightTop';
      await until(() => legend().rendered, 600);
      const movedLegend = legend();
      h.props.splitColumnName = ''; h.props.legendVisibility = 'Auto';
      await until(() => legend().labels.length === 0, 600);

      const shownText = () =>
        ((h.root as HTMLElement).innerText || '').replace(/\s+/g, ' ').trim();
      h.props.showTitle = true;
      h.props.title = 'Age Distribution';
      h.props.description = 'Shows distribution of patient ages';
      h.props.descriptionVisibilityMode = 'Always';
      h.props.descriptionPosition = 'Bottom';
      await until(() => shownText().includes('Shows distribution of patient ages'), 600);
      const withDescription = shownText();
      h.props.description = ''; h.props.title = ''; h.props.showTitle = false;
      await until(() => !shownText().includes('Shows distribution of patient ages'), 600);
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

  await softStep('Context menu', async () => {
    const result = await page.evaluate(async () => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      const canvas = document.querySelector('[name="viewer-Histogram"] canvas')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));

      const t0 = Date.now();
      let prev = -1, stable = 0;
      while (Date.now() - t0 < 600) {
        const n = document.querySelectorAll('.d4-menu-item-label').length;
        stable = (n > 0 && n === prev) ? stable + 1 : 0;
        prev = n;
        if (stable >= 2) break;
        await new Promise(r => setTimeout(r, 25));
      }
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label')).map(e => e.textContent!.trim());
      const before = h.props.showFilteredOutRows;
      const sfor = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(e => e.textContent!.trim() === 'Show Filtered Out Rows');
      if (sfor) (sfor.closest('.d4-menu-item') as HTMLElement).click();
      for (let i = 0; i < 20 && h.props.showFilteredOutRows === before; i++)
        await new Promise(r => setTimeout(r, 25));
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

  await softStep('Data — row source', async () => {
    const errBefore = errorCount();
    await v.closeAllAndWait(page);
    await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      df.name = 'SPGI';
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });
      (document.querySelector('[name="icon-histogram"]') as HTMLElement).click();
      for (let i = 0; i < 48; i++) {
        const hv = tv.viewers.find(v => v.type === 'Histogram');
        if (hv && hv.root.querySelector('canvas')) break;
        await new Promise(r => setTimeout(r, 25));
      }
      const selected = new Promise(resolve => {
        const sub = df.onSelectionChanged.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 400);
      });
      for (let i = 0; i < 5; i++) df.selection.set(i, true);
      await selected;
    });
    await v.setViewerProps(page, 'Histogram', [{set: {rowSource: 'All'}, wait: 600}]);
    const allPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    const allHue = await v.countSelectionHuePixels(page, 'Histogram');
    await v.setViewerProps(page, 'Histogram', [{set: {rowSource: 'Selected'}, wait: 600}]);
    const selectedPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    const selectedHue = await v.countSelectionHuePixels(page, 'Histogram');
    await v.setViewerProps(page, 'Histogram', [{set: {rowSource: 'All'}, wait: 500}]);
    const restoredPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    const restoredHue = await v.countSelectionHuePixels(page, 'Histogram');

    expect(allPx).toBeGreaterThan(1000);
    expect(allPx - selectedPx).toBeGreaterThan(1000);
    expect(Math.abs(restoredPx - allPx)).toBeLessThan(500);
    expect(allHue).toBeGreaterThan(100);
    expect(selectedHue).toBe(0);
    expect(restoredHue).toBeGreaterThan(100);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Data — filter formula', async () => {
    const errBefore = errorCount();
    await v.closeAllAndWait(page);
    await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });
      (document.querySelector('[name="icon-histogram"]') as HTMLElement).click();
      for (let i = 0; i < 48; i++) {
        const hv = tv.viewers.find(v => v.type === 'Histogram');
        if (hv && hv.root.querySelector('canvas')) break;
        await new Promise(r => setTimeout(r, 25));
      }
    });
    await v.setViewerProps(page, 'Histogram', [{set: {valueColumnName: 'AGE', filter: ''}, wait: 600}]);

    let prevPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    let clearedPx = prevPx;
    for (let i = 0; i < 5; i++) {

      await page.waitForTimeout(300);
      clearedPx = (await v.countCanvasPixels(page, 'Histogram')).total;
      if (Math.abs(clearedPx - prevPx) < 200) break;
      prevPx = clearedPx;
    }
    expect(Math.abs(clearedPx - prevPx)).toBeLessThan(200);
    const cleared = await page.evaluate(() => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      return {rowCount: grok.shell.tv.dataFrame.rowCount, trueCount: (h as any).filter.trueCount};
    });
    await v.setViewerProps(page, 'Histogram', [{set: {filter: '${AGE} > 40'}, wait: 700}]);
    const filteredPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    const filteredCount = await page.evaluate(() => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      return (h as any).filter.trueCount;
    });
    await v.setViewerProps(page, 'Histogram', [{set: {filter: ''}, wait: 600}]);
    const restoredPx = (await v.countCanvasPixels(page, 'Histogram')).total;
    const restoredCount = await page.evaluate(() => {
      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      return (h as any).filter.trueCount;
    });
    expect(clearedPx).toBeGreaterThan(1000);
    expect(clearedPx - filteredPx).toBeGreaterThan(400);
    expect(Math.abs(restoredPx - clearedPx)).toBeLessThan(300);

    expect(cleared.trueCount).toBe(cleared.rowCount);
    expect(filteredCount).toBeGreaterThan(0);
    expect(filteredCount).toBeLessThan(cleared.rowCount);
    expect(restoredCount).toBe(cleared.rowCount);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Data — table switching', async () => {
    await v.closeAllAndWait(page);
    const result = await page.evaluate(async () => {
      const dfSpgi = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      dfSpgi.name = 'SPGI';
      grok.shell.addTableView(dfSpgi);
      await new Promise(resolve => {
        const sub = dfSpgi.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });
      (document.querySelector('[name="icon-histogram"]') as HTMLElement).click();
      for (let i = 0; i < 48; i++) {
        const hv = grok.shell.tv.viewers.find(v => v.type === 'Histogram');
        if (hv && hv.root.querySelector('canvas')) break;
        await new Promise(r => setTimeout(r, 25));
      }

      const dfDemog = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      dfDemog.name = 'demog';
      grok.shell.addTableView(dfDemog);
      for (let i = 0; i < 32; i++) {
        if (Array.from(grok.shell.tableViews).some((t: any) => t.dataFrame?.name === 'demog')) break;
        await new Promise(r => setTimeout(r, 25));
      }

      const spgiView = Array.from(grok.shell.tableViews).find((t: any) => t.dataFrame?.name === 'SPGI');
      if (spgiView) grok.shell.v = spgiView as any;
      for (let i = 0; i < 20 && grok.shell.tv?.dataFrame?.name !== 'SPGI'; i++)
        await new Promise(r => setTimeout(r, 25));

      const h = grok.shell.tv.viewers.find(v => v.type === 'Histogram')!;
      const before = {dfName: h.dataFrame?.name, value: h.props.valueColumnName};
      h.props.table = 'demog';
      for (let i = 0; i < 40 && h.dataFrame?.name === before.dfName; i++)
        await new Promise(r => setTimeout(r, 25));
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
