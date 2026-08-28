/* ---
realizes: [densityplot.cp.setup-binning-color, viewers.density-plot]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Unable to find element in cloned iframe/.test(text) ||
  /Stack trace [A-Za-z]+/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text);

async function setPropSelect(page: Page, rowName: string, value: string) {
  await page.evaluate((args: {rowName: string; value: string}) => {
    const row = document.querySelector(`[name="${args.rowName}"]`) as HTMLElement;
    (row?.querySelector('[name^="prop-view-"]') as HTMLElement)?.click();
    const sel = row?.querySelector('select.property-grid-item-editor-spinner') as HTMLSelectElement;
    if (sel) { sel.value = args.value; sel.dispatchEvent(new Event('change', {bubbles: true})); }
  }, {rowName, value});
  await page.keyboard.press('Enter');
  await v.waitForViewerRendered(page, 'Density plot', 300);
}

async function setBinsViaSlider(page: Page, bins: number) {
  await page.evaluate((n: number) => {
    const slider = document.querySelector(
      '[name="viewer-Density-plot"] input[type="range"]') as HTMLInputElement;
    if (slider) {
      slider.value = String(n);
      slider.dispatchEvent(new Event('input', {bubbles: true}));
      slider.dispatchEvent(new Event('change', {bubbles: true}));
    }
  }, bins);
  await v.pollValue(() => page.evaluate(() => {
    const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
    return d?.props.bins;
  }), (b) => b === bins, 400, 100);
}

test('Density Plot — Setup, Axis Columns, Binning, Color Mapping, Persistence', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text()) && !isLocalBootNoise(m.text()))
      consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;

  const settledPx = async () => {
    let prev = (await v.countCanvasPixels(page, 'Density plot')).total;
    let cur = prev;
    for (let i = 0; i < 8; i++) {
      await page.waitForTimeout(350);
      cur = (await v.countCanvasPixels(page, 'Density plot')).total;
      if (Math.abs(cur - prev) < 500) break;
      prev = cur;
    }
    return cur;
  };

  const selectorLabel = (axis: 'x' | 'y') => page.evaluate((a: string) => {
    const el = document.querySelector(
      `[name="viewer-Density-plot"] [name="div-column-combobox-${a}"] .d4-column-selector-column`);
    return (el?.textContent ?? '').trim();
  }, axis);

  const readProp = (prop: string) => page.evaluate((p: string) => {
    const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
    return d?.props[p];
  }, prop);

  const readXY = () => page.evaluate(() => {
    const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
    return {x: d.props.xColumnName, y: d.props.yColumnName};
  });

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'density-plot', 'Density-plot');

  await softStep('Scenario 1 — pick axis columns through the on-viewer selectors (GROK-16612)', async () => {
    const errBefore = errCount();
    // the on-viewer column combobox ignores a synthetic mousedown, so the non-trusted
    // picker silently leaves the column unchanged — drive it with a real click
    const pick = (axis: 'x' | 'y', col: string) => v.pickColumnViaSelectorTrusted(page, {
      role: axis,
      columnName: col,
      scopeSelector: '[name="viewer-Density-plot"]',
      viewerType: 'Density plot',
      propName: axis === 'x' ? 'xColumnName' : 'yColumnName',
    });

    await pick('x', 'WEIGHT');
    expect((await readXY()).x).toBe('WEIGHT');
    expect(await selectorLabel('x')).toBe('WEIGHT');
    await pick('y', 'WEIGHT');
    expect((await readXY()).y).toBe('WEIGHT');
    await pick('y', 'HEIGHT');
    expect((await readXY()).y).toBe('HEIGHT');
    await pick('x', 'AGE');
    const finalXY = await readXY();
    expect(finalXY.x).toBe('AGE');
    expect(finalXY.y).toBe('HEIGHT');
    expect(await selectorLabel('x')).toBe('AGE');
    expect(await selectorLabel('y')).toBe('HEIGHT');

    expect(errCount()).toBe(errBefore);
  });

  await page.evaluate(() => {
    const viewer = document.querySelector('[name="viewer-Density-plot"]') as HTMLElement;
    const gear = viewer?.parentElement?.parentElement
      ?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear?.click();
  });
  await v.pollValue(() => page.locator('.property-grid').count(), (n) => n > 0, 1000, 100);

  await softStep('Scenario 2 — bin count lower boundary (bins=1), 5 vs 200 (strong repaint), and bin shape', async () => {
    const errBefore = errCount();

    await setBinsViaSlider(page, 50);
    const px50ref = await settledPx();
    await setBinsViaSlider(page, 1);
    const bins1Prop = await readProp('bins');
    const px1 = await settledPx();
    console.log(`DensityPlot bins=1 edge: px50ref=${px50ref} px1=${px1} delta=${px1 - px50ref}`);
    expect(bins1Prop).toBe(1);
    expect(px1).toBeGreaterThan(0);

    expect(px1 - px50ref).toBeGreaterThan(70000);
    await setBinsViaSlider(page, 50);
    const px50back = await settledPx();
    expect(await readProp('bins')).toBe(50);

    expect(px1 - px50back).toBeGreaterThan(70000);

    await setBinsViaSlider(page, 5);
    const bins5Prop = await readProp('bins');
    const px5 = await settledPx();
    await setBinsViaSlider(page, 200);
    const bins200Prop = await readProp('bins');
    const px200 = await settledPx();
    console.log(`DensityPlot bins px: px5=${px5} px200=${px200} delta=${px5 - px200}`);
    expect(bins5Prop).toBe(5);
    expect(bins200Prop).toBe(200);

    expect(px5).toBeGreaterThan(0);
    expect(px200).toBeGreaterThan(0);

    expect(px5 - px200).toBeGreaterThan(100000);

    await setPropSelect(page, 'prop-bin-shape', 'rectangle');
    const shape = await readProp('binShape');
    const pxRect = await settledPx();
    console.log(`DensityPlot shape px: pxRect=${pxRect}`);
    expect(shape).toBe('rectangle');
    expect(pxRect).toBeGreaterThan(0);
    expect(errCount()).toBe(errBefore);

  });

  await softStep('Scenario 3 — Invert Color Scheme (strong repaint) and Color Transform Type', async () => {
    const errBefore = errCount();
    const pxBefore = await settledPx();
    await page.evaluate(() => {
      (document.querySelector(
        '[name="prop-invert-color-scheme"] input[type="checkbox"]') as HTMLInputElement)?.click();
    });
    const invert = await v.pollValue(() => readProp('invertColorScheme'), (b) => b === true, 400, 100);
    const pxInvert = await settledPx();
    console.log(`DensityPlot invert px: pxBefore=${pxBefore} pxInvert=${pxInvert} delta=${pxInvert - pxBefore}`);
    expect(invert).toBe(true);

    expect(Math.abs(pxInvert - pxBefore)).toBeGreaterThan(200000);

    await setPropSelect(page, 'prop-color-transform-type', 'logarithmic');
    let ct = await readProp('colorTransformType');
    expect(ct).toBe('logarithmic');
    await setPropSelect(page, 'prop-color-transform-type', 'linear');
    await setPropSelect(page, 'prop-color-transform-type', 'logarithmic');
    ct = await readProp('colorTransformType');
    expect(ct).toBe('logarithmic');
    expect(errCount()).toBe(errBefore);

  });

  await softStep('Scenario 4 — GROK-17118 guard: X selector offers numeric columns only', async () => {
    const shape = await readProp('binShape');
    expect(shape).toBe('rectangle');
    const errBefore = errCount();
    const xBefore = (await readXY()).x;

    await page.evaluate(() => {
      const sel = document.querySelector(
        '[name="viewer-Density-plot"] [name="div-column-combobox-x"]')!;
      (sel.querySelector('.d4-column-selector-column') || sel)
        .dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
    });

    const backdropCount = () => page.locator('.d4-column-selector-backdrop').count();
    await v.pollValue(backdropCount, (n) => n > 0, 3000, 100);

    await page.keyboard.press('s');
    await page.waitForTimeout(100); 
    await page.keyboard.type('ex');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(300); 

    const xAfter = (await readXY()).x;

    await page.evaluate(() => document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true})));
    await v.pollValue(backdropCount, (n) => n === 0, 200, 100);
    expect(xBefore).toBe('AGE');
    expect(xAfter).toBe('AGE');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 5 — degenerate zero-width bin range holds an error-free floor', async () => {

    const constCol = 'DP_ZW_CONST_' + Date.now();
    await page.evaluate((name: string) => {
      const df = grok.shell.tv.dataFrame;
      if (!df.columns.names().includes(name)) df.columns.addNewCalculated(name, '42');
    }, constCol);
    await v.pollValue(
      () => page.evaluate((name: string) => grok.shell.tv.dataFrame.columns.names().includes(name), constCol),
      (present) => present, 800, 100);
    try {
      const errBefore = errCount();

      await v.pickColumnViaSelector(page, {
        comboboxSuffix: 'x',
        columnName: constCol,
        scopeSelector: '[name="viewer-Density-plot"]',
        popupWaitStrategy: 'backdrop',
        viewerType: 'Density plot',
        propName: 'xColumnName',
        allowFallback: false,
      });
      expect((await readXY()).x).toBe(constCol);

      await setPropSelect(page, 'prop-bin-shape', 'rectangle');
      const pxRect = await settledPx();
      const rectAttached = await page.evaluate(() => {
        const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
        return document.body.contains(d.root);
      });
      console.log(`DensityPlot zero-width rect: px=${pxRect} shape=${await readProp('binShape')}`);
      expect(await readProp('binShape')).toBe('rectangle');
      expect(rectAttached).toBe(true);
      expect(pxRect).toBeGreaterThan(0);

      await setPropSelect(page, 'prop-bin-shape', 'hexagon');
      const pxHex = await settledPx();
      const hexAttached = await page.evaluate(() => {
        const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
        return document.body.contains(d.root);
      });
      console.log(`DensityPlot zero-width hex: px=${pxHex} shape=${await readProp('binShape')}`);
      expect(await readProp('binShape')).toBe('hexagon');
      expect(hexAttached).toBe(true);
      expect(pxHex).toBeGreaterThan(0);

      expect(errCount()).toBe(errBefore);
    } finally {

      await page.evaluate(() => {
        const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
        if (d) { d.props.xColumnName = 'AGE'; d.props.binShape = 'rectangle'; }
      });

      await v.waitForViewerRendered(page, 'Density plot', 400);
      await page.evaluate((name: string) => {
        const df = grok.shell.tv?.dataFrame;
        if (df && df.columns.names().includes(name)) df.columns.remove(name);
      }, constCol);
    }
  });

  await page.evaluate(() => grok.shell.closeAll());

  v.finishSpec();
});
