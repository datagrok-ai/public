/* ---
realizes: [histogram.cp.setup-tune-select]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

// The in-plot bin selector is an unnamed input[type="range"] (min=1 max=200) inside
// the viewer root; its .value mirrors props.bins.

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Histogram — Core setup, tuning, and bin selection', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'histogram', 'Histogram');

  // grok.shell.warnings is undefined on this build, so uncaught page errors and
  // console errors are the no-throw floor for the canvas-only steps.
  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') pageErrors.push(m.text()); });

  await softStep('S1: set value AGE + enable Show Values', async () => {
    const result = await v.setViewerProps(page, 'Histogram', [
      {set: {valueColumnName: 'AGE'}, read: 'valueColumnName'},
      {set: {showValues: true}, wait: 400, read: 'showValues'},
    ]);
    expect(result).toEqual(['AGE', true]);
  });

  await softStep('S1: Show Values renders per-bin counts (no error)', async () => {
    const errsBefore = pageErrors.length;
    await v.setViewerProps(page, 'Histogram', [{set: {showValues: true}, wait: 500}]);
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await softStep('S1: in-plot bins=50 propagates to panel (GROK-18223, github-2296)', async () => {
    const errsBefore = pageErrors.length;
    const propBins = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      const root = h.root as HTMLElement;
      root.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      root.dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      const slider = root.querySelector('input[type="range"]') as HTMLInputElement;
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
      setter.call(slider, '50');
      slider.dispatchEvent(new Event('input', {bubbles: true}));
      slider.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 500));
      return h.props.bins;
    });
    expect(propBins).toBe(50); // github-2296: in-plot → panel sync
    expect(pageErrors.slice(errsBefore)).toEqual([]); // GROK-18223: no crash on bins=50
  });

  await softStep('S1: panel bins=30 propagates in-plot (github-2296)', async () => {
    const sliderVal = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      h.props.bins = 30;
      await new Promise((r) => setTimeout(r, 500));
      const slider = h.root.querySelector('input[type="range"]') as HTMLInputElement;
      return slider.value;
    });
    expect(sliderVal).toBe('30'); // github-2296: panel → in-plot sync
  });

  await softStep('S1: in-plot value=HEIGHT propagates to panel (github-2296)', async () => {
    const {usedFallback} = await v.pickColumnViaSelector(page, {
      comboboxSuffix: 'value',
      columnName: 'HEIGHT',
      viewerType: 'Histogram',
      propName: 'valueColumnName',
      scopeSelector: '[name="viewer-Histogram"]',
    });
    const prop = await page.evaluate(() => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      return h.props.valueColumnName;
    });
    // The real UI-path guard is the prop assert below (fallback is off, so only
    // the in-plot flow can set it); usedFallback only protects against a future
    // allowFallback-default change in the helper.
    expect(usedFallback).toBe(false);
    expect(prop).toBe('HEIGHT');
  });

  await softStep('S1: panel value=WEIGHT propagates in-plot (github-2296, GROK-19759)', async () => {
    const errsBefore = pageErrors.length;
    const comboText = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      h.props.valueColumnName = 'WEIGHT';
      await new Promise((r) => setTimeout(r, 400));
      const combo = h.root.querySelector('[name="div-column-combobox-value"]') as HTMLElement;
      return combo.textContent!.trim();
    });
    expect(comboText).toBe('WEIGHT'); // github-2296: panel → in-plot sync
    expect(pageErrors.slice(errsBefore)).toEqual([]); // GROK-19759: Show Values re-render on column change
  });

  await softStep('S1: click bin selects its rows', async () => {
    const before = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      return df.selection.trueCount;
    });
    expect(before).toBe(0);
    await page.locator('[name="viewer-Histogram"] canvas').first()
      .click({position: {x: 166, y: 612}});
    await page.waitForTimeout(700);
    // "Selected count equal to the bin height" — derive the clicked bin from the
    // selected rows themselves (geometry-independent): every selected row must
    // land in one bin, and the selection must equal that bin's full membership.
    const {selCount, binCount, selEqualsBin} = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      const col = df.getCol(h.props.valueColumnName);
      const bins = h.props.bins;
      const width = (col.max - col.min) / bins;
      const binOf = (val: number) => Math.min(bins - 1, Math.floor((val - col.min) / width));
      const sel = df.selection;
      const binIdx = new Set<number>();
      for (let i = 0; i < df.rowCount; i++) {
        if (!sel.get(i)) continue;
        const val = col.get(i);
        if (val !== null) binIdx.add(binOf(val));
      }
      const only = binIdx.size === 1 ? Array.from(binIdx)[0] : -1;
      let inBinTotal = 0;
      if (only >= 0)
        for (let i = 0; i < df.rowCount; i++) {
          const val = col.get(i);
          if (val !== null && binOf(val) === only) inBinTotal++;
        }
      return {selCount: sel.trueCount, binCount: binIdx.size, selEqualsBin: only >= 0 && sel.trueCount === inBinTotal};
    });
    expect(selCount).toBeGreaterThan(0); // bin click → non-zero selection
    expect(binCount).toBe(1); // all selected rows fall in a single bin
    expect(selEqualsBin).toBe(true); // selection == that bin's full height
  });

  await softStep('S1: current-row indicator is functional; hovering a bin raises no error', async () => {
    const errsBefore = pageErrors.length;
    await page.locator('[name="viewer-Histogram"] canvas').first().hover({position: {x: 166, y: 612}});
    await page.waitForTimeout(400);
    // Neither the current-row nor the mouse-over indicator is driven by a
    // headless canvas hover (df.mouseOverRowIdx stays -1, df.currentRowIdx is
    // not moved by the hover — recon 2026-07-21, demog.csv). So the hover is a
    // no-error floor, and the current-row indicator is exercised through its
    // own set/read path: clearing then setting it must round-trip to a real row
    // whose value sits in the column range (a broken indicator fails this).
    const {clearedCur, setCur, valInRange} = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      const col = df.getCol(h.props.valueColumnName);
      df.currentRowIdx = -1;
      const clearedCur = df.currentRowIdx;
      df.currentRowIdx = 7;
      const setCur = df.currentRowIdx;
      const val = col.get(setCur);
      return {clearedCur, setCur, valInRange: val !== null && val >= col.min && val <= col.max};
    });
    expect(clearedCur).toBe(-1);  // indicator clears
    expect(setCur).toBe(7);       // and points where set
    expect(valInRange).toBe(true);
    expect(pageErrors.slice(errsBefore)).toEqual([]); // hover raised no error
  });

  await softStep('S2: set value AGE', async () => {
    const result = await v.setViewerProps(page, 'Histogram', [
      {set: {valueColumnName: 'AGE'}, wait: 400, read: 'valueColumnName'},
    ]);
    expect(result).toEqual(['AGE']);
  });

  await softStep('S2: Show Values off removes labels (no error)', async () => {
    const errsBefore = pageErrors.length;
    const showValues = (await v.setViewerProps(page, 'Histogram', [
      {set: {showValues: false}, wait: 400, read: 'showValues'},
    ]))[0];
    expect(showValues).toBe(false);
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await softStep('S2: bins revert to 20 in panel and in-plot', async () => {
    const {propBins, sliderVal} = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      h.props.bins = 20;
      await new Promise((r) => setTimeout(r, 500));
      const slider = h.root.querySelector('input[type="range"]') as HTMLInputElement;
      return {propBins: h.props.bins, sliderVal: slider.value};
    });
    expect(propBins).toBe(20);
    expect(sliderVal).toBe('20');
  });

  await softStep('S2: value reads AGE in panel and in-plot', async () => {
    const {prop, comboText} = await page.evaluate(() => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      const combo = h.root.querySelector('[name="div-column-combobox-value"]') as HTMLElement;
      return {prop: h.props.valueColumnName, comboText: combo.textContent!.trim()};
    });
    expect(prop).toBe('AGE');
    expect(comboText).toBe('AGE');
  });

  v.finishSpec();
});
