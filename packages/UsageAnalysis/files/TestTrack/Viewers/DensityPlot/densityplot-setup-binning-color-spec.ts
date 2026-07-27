/* ---
realizes: [densityplot.cp.setup-binning-color, viewers.density-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Unable to find element in cloned iframe/.test(text) ||
  /Stack trace [A-Za-z]+/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text);

/** Set a property-grid <select> row: open its inline editor, set value, commit. */
async function setPropSelect(page: Page, rowName: string, value: string) {
  await page.evaluate((args: {rowName: string; value: string}) => {
    const row = document.querySelector(`[name="${args.rowName}"]`) as HTMLElement;
    (row?.querySelector('[name^="prop-view-"]') as HTMLElement)?.click();
    const sel = row?.querySelector('select.property-grid-item-editor-spinner') as HTMLSelectElement;
    if (sel) { sel.value = args.value; sel.dispatchEvent(new Event('change', {bubbles: true})); }
  }, {rowName, value});
  await page.keyboard.press('Enter');
  await page.waitForTimeout(300);
}

/** Drive the on-viewer bin-count range slider (bottom-left) via input+change events. */
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
  await page.waitForTimeout(400);
}

test('Density Plot — Setup, Axis Columns, Binning, Color Mapping, Persistence', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text()))
      consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;

  // Settle-gated canvas ink: repeat until two reads agree, so a delta between
  // measurements is the setter's effect, not a render tail.
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

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'density-plot', 'Density-plot');

  await softStep('Scenario 1 — pick axis columns through the on-viewer selectors (GROK-16612)', async () => {
    const errBefore = errCount();
    const pick = (axis: 'x' | 'y', col: string) => v.pickColumnViaSelector(page, {
      comboboxSuffix: axis,
      columnName: col,
      scopeSelector: '[name="viewer-Density-plot"]',
      popupWaitStrategy: 'backdrop',
      viewerType: 'Density plot',
      propName: axis === 'x' ? 'xColumnName' : 'yColumnName',
      allowFallback: false,
    });
    // Real UI path only (no JS-API fallback): a broken selector fails here
    // instead of being masked. Product read is the prop plus the DOM label.
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
    // GROK-16612: no console errors across the axis switches.
    expect(errCount()).toBe(errBefore);
  });

  // Open the settings panel — the gear sits two parent hops above the viewer root.
  await page.evaluate(() => {
    const viewer = document.querySelector('[name="viewer-Density-plot"]') as HTMLElement;
    const gear = viewer?.parentElement?.parentElement
      ?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear?.click();
  });
  await page.waitForTimeout(1000);

  await softStep('Scenario 2 — bin count lower boundary (bins=1), 5 vs 200 (strong repaint), and bin shape', async () => {
    const errBefore = errCount();

    // Bin-count lower boundary (atlas edge "Bin-count boundary"): bins=1 collapses
    // the plot to a single filled cell, so the ink jumps sharply above the 50-bin
    // grid; reverting to 50 drops it back (round-trip). The bin slider clamps at
    // min=1, so this is the true lower edge of the range.
    await setBinsViaSlider(page, 50);
    const px50ref = await settledPx();
    await setBinsViaSlider(page, 1);
    const bins1Prop = await readProp('bins');
    const px1 = await settledPx();
    console.log(`DensityPlot bins=1 edge: px50ref=${px50ref} px1=${px1} delta=${px1 - px50ref}`);
    expect(bins1Prop).toBe(1);
    expect(px1).toBeGreaterThan(0);
    // The floor keeps a ~2x margin below the smallest observed bins=1 vs bins=50 delta.
    expect(px1 - px50ref).toBeGreaterThan(70000);
    await setBinsViaSlider(page, 50);
    const px50back = await settledPx();
    expect(await readProp('bins')).toBe(50);
    // Reverting to 50 bins drops the fill back sharply (round-trip).
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
    // Bins 5 -> 200 is the strongest signal on this viewer: a few large dark
    // bins ink far more than a fine sparse grid, so px5 dominates px200.
    expect(px5).toBeGreaterThan(0);
    expect(px200).toBeGreaterThan(0);
    // The floor keeps a ~2x margin below the observed 5-vs-200-bin delta.
    expect(px5 - px200).toBeGreaterThan(100000);

    // Bin shape hexagon -> rectangle: modest delta at high bin counts, so the
    // switch holds an error-free floor with the prop confirmed.
    await setPropSelect(page, 'prop-bin-shape', 'rectangle');
    const shape = await readProp('binShape');
    const pxRect = await settledPx();
    console.log(`DensityPlot shape px: pxRect=${pxRect}`);
    expect(shape).toBe('rectangle');
    expect(pxRect).toBeGreaterThan(0);
    expect(errCount()).toBe(errBefore);
    // Leave bins 200, rectangle (peak configuration for the persistence tail).
  });

  await softStep('Scenario 3 — Invert Color Scheme (strong repaint) and Color Transform Type', async () => {
    const errBefore = errCount();
    const pxBefore = await settledPx();
    await page.evaluate(() => {
      (document.querySelector(
        '[name="prop-invert-color-scheme"] input[type="checkbox"]') as HTMLInputElement)?.click();
    });
    await page.waitForTimeout(400);
    const invert = await readProp('invertColorScheme');
    const pxInvert = await settledPx();
    console.log(`DensityPlot invert px: pxBefore=${pxBefore} pxInvert=${pxInvert} delta=${pxInvert - pxBefore}`);
    expect(invert).toBe(true);
    // Inverting the color scheme flips the binned-area background — a very strong
    // repaint; the floor keeps a ~2x margin below the observed invert delta.
    expect(Math.abs(pxInvert - pxBefore)).toBeGreaterThan(200000);

    // Color transform linear -> log -> linear -> log: modest luminance-only
    // change, held as an error-free floor with the prop confirmed.
    await setPropSelect(page, 'prop-color-transform-type', 'logarithmic');
    let ct = await readProp('colorTransformType');
    expect(ct).toBe('logarithmic');
    await setPropSelect(page, 'prop-color-transform-type', 'linear');
    await setPropSelect(page, 'prop-color-transform-type', 'logarithmic');
    ct = await readProp('colorTransformType');
    expect(ct).toBe('logarithmic');
    expect(errCount()).toBe(errBefore);
    // Leave Invert on, transform logarithmic.
  });

  await softStep('Scenario 4 — GROK-17118 guard: X selector offers numeric columns only', async () => {
    const shape = await readProp('binShape');
    expect(shape).toBe('rectangle');
    const errBefore = errCount();
    const xBefore = (await readXY()).x;
    // Open the X selector popup (mousedown on the label opens it).
    await page.evaluate(() => {
      const sel = document.querySelector(
        '[name="viewer-Density-plot"] [name="div-column-combobox-x"]')!;
      (sel.querySelector('.d4-column-selector-column') || sel)
        .dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
    });
    await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
      null, {timeout: 3000}).catch(() => {});
    // Type a non-numeric column name — SEX matches nothing (numeric-only filter).
    await page.keyboard.press('s');
    await page.waitForTimeout(100);
    await page.keyboard.type('ex');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(300);
    // Enter committed nothing — X stays AGE (the numeric-only selector IS the fix).
    const xAfter = (await readXY()).x;
    // Close the (possibly still-open) no-match popup via an outside mousedown.
    await page.evaluate(() => document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true})));
    await page.waitForTimeout(200);
    expect(xBefore).toBe('AGE');
    expect(xAfter).toBe('AGE');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 5 — degenerate zero-width bin range holds an error-free floor', async () => {
    // Atlas edge "Degenerate zero-width bin range": a constant column has
    // min == max, so its bin range has zero width — the bin-size divisor is
    // guarded against divide-by-zero. Build a numeric constant fixture column
    // via the JS API, assign it to X through the real on-viewer selector (a
    // numeric column IS offered), and confirm both bin shapes render without
    // error or freeze. The fixture column is removed in finally.
    const constCol = 'DP_ZW_CONST_' + Date.now();
    await page.evaluate(async (name: string) => {
      const df = grok.shell.tv.dataFrame;
      if (!df.columns.names().includes(name)) df.columns.addNewCalculated(name, '42');
      await new Promise((r) => setTimeout(r, 800));
    }, constCol);
    try {
      const errBefore = errCount();
      // Assign the constant column on X through the on-viewer selector (real UI path).
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

      // Rectangle over the zero-width range: no error, viewer stays attached,
      // and a settled repaint completes (degenerate strip).
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

      // Hexagon over the same zero-width range: same error-free floor.
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

      // The zero-width range neither errored nor froze either bin shape.
      expect(errCount()).toBe(errBefore);
    } finally {
      // Restore the peak-config X column and rectangle shape, then drop the fixture.
      await page.evaluate(async (name: string) => {
        const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
        if (d) { d.props.xColumnName = 'AGE'; d.props.binShape = 'rectangle'; }
        await new Promise((r) => setTimeout(r, 400));
        const df = grok.shell.tv?.dataFrame;
        if (df && df.columns.names().includes(name)) df.columns.remove(name);
      }, constCol);
    }
  });

  await softStep('Scenario 6a — layout round-trip restores the saved viewer set and config', async () => {
    // Pin the peak configuration deterministically before persisting.
    const layoutId = await page.evaluate(async () => {
      const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      d.props.xColumnName = 'AGE';
      d.props.yColumnName = 'HEIGHT';
      d.props.bins = 200;
      d.props.binShape = 'rectangle';
      d.props.invertColorScheme = true;
      await new Promise((r) => setTimeout(r, 800));
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });
    try {
      const result = await page.evaluate(async (id) => {
        const tv = grok.shell.tv;
        await new Promise((r) => setTimeout(r, 1000));
        tv.addViewer('Scatter plot');
        await new Promise((r) => setTimeout(r, 600));
        const saved = await grok.dapi.layouts.find(id);
        tv.loadLayout(saved);
        await new Promise((r) => setTimeout(r, 3000));
        const hasScatter = tv.viewers.some((vw) => vw.type === 'Scatter plot');
        const hasDensity = tv.viewers.some((vw) => vw.type === 'Density plot');
        const dp = tv.viewers.find((vw) => vw.type === 'Density plot');
        return {
          hasScatter, hasDensity,
          x: dp?.props.xColumnName, y: dp?.props.yColumnName,
          bins: dp?.props.bins, binShape: dp?.props.binShape,
          invert: dp?.props.invertColorScheme,
        };
      }, layoutId);
      // Restored set equals the SAVED set: Density present, later Scatter absent.
      expect(result.hasDensity).toBe(true);
      expect(result.hasScatter).toBe(false);
      // The restored Density Plot carries the persisted configuration.
      expect(result.x).toBe('AGE');
      expect(result.y).toBe('HEIGHT');
      expect(result.bins).toBe(200);
      expect(result.binShape).toBe('rectangle');
      expect(result.invert).toBe(true);
    } finally {
      await page.evaluate(async (id) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved) await grok.dapi.layouts.delete(saved);
        } catch (_) {}
      }, layoutId);
    }
  });

  await softStep('Scenario 6b — project save / Close All / reopen restores the Density Plot', async () => {
    const projName = 'zz-densityplot-persistence-probe-' + Date.now();
    let projectId: string | null = null;
    try {
      const saved = await saveProjectViaUI(page, projName);
      projectId = saved.projectId;
      expect(projectId).toBeTruthy();

      const result = await page.evaluate(async (id) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const full = await grok.dapi.projects.find(id);
        await full.open();
        // Viewers re-materialize asynchronously — poll for the restored Density plot.
        let types: string[] = [];
        for (let t = 0; t < 20; t++) {
          await new Promise((r) => setTimeout(r, 1000));
          types = [];
          for (const view of grok.shell.tableViews)
            for (const vw of view.viewers) types.push(vw.type);
          if (types.includes('Density plot')) break;
        }
        let dp: any = null;
        for (const view of grok.shell.tableViews)
          for (const vw of view.viewers)
            if (vw.type === 'Density plot') dp = vw;
        return {
          types,
          x: dp?.props.xColumnName, y: dp?.props.yColumnName,
          bins: dp?.props.bins, binShape: dp?.props.binShape,
          invert: dp?.props.invertColorScheme,
        };
      }, projectId);

      // Cross-session round-trip: a Density Plot is restored with its config.
      expect(result.types).toContain('Density plot');
      expect(result.x).toBe('AGE');
      expect(result.y).toBe('HEIGHT');
      expect(result.bins).toBe(200);
      expect(result.binShape).toBe('rectangle');
      expect(result.invert).toBe(true);
    } finally {
      if (projectId)
        await deleteProjectWithCleanup(page, {projectId});
    }
  });

  await page.evaluate(() => grok.shell.closeAll());

  v.finishSpec();
});
