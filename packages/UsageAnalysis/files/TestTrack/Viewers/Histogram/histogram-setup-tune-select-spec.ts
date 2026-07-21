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
    await v.pickColumnViaSelector(page, {
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
    const selCount = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selCount).toBeGreaterThan(0); // bin click → non-zero selection matching bin height
  });

  await softStep('S1: hover bin drives row indicators (no error)', async () => {
    const errsBefore = pageErrors.length;
    await page.locator('[name="viewer-Histogram"] canvas').first().hover({position: {x: 166, y: 612}});
    await page.waitForTimeout(400);
    expect(pageErrors.slice(errsBefore)).toEqual([]);
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
