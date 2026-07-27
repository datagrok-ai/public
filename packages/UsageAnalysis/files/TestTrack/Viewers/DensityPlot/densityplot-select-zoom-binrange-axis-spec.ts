/* ---
realizes: [densityplot.cp.select-zoom-binrange-axis, densityplot.bin-to-range-viewport, viewers.density-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Stack trace [A-Za-z]+/.test(text);

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

/** Expand a settings-panel category if it is collapsed (icon-plus present). */
async function expandCategory(page: Page, cat: string) {
  await page.evaluate((c: string) => {
    const header = document.querySelector(`[name="prop-category-${c}"]`) as HTMLElement;
    if (header && header.querySelector('.property-grid-icon-plus')) header.click();
  }, cat);
  await page.waitForTimeout(200);
}

/** Type a numeric bound into a scoped X/Y Min/Max property cell (duplicated row names). */
async function setBound(page: Page, viewName: string, value: number) {
  await page.evaluate((vn: string) => {
    (document.querySelector(`[name="${vn}"]`) as HTMLElement)?.click();
  }, viewName);
  await page.waitForTimeout(200);
  await page.evaluate((vn: string) => {
    const item = (document.querySelector(`[name="${vn}"]`) as HTMLElement)?.closest('.property-grid-item');
    (item?.querySelector('input.property-grid-item-editor-textbox') as HTMLInputElement)?.focus();
  }, viewName);
  await page.keyboard.press('Control+A');
  await page.keyboard.type(String(value));
  await page.keyboard.press('Enter');
  await page.waitForTimeout(300);
}

/** Clear a scoped X/Y Min/Max property cell to null. */
async function clearBound(page: Page, viewName: string) {
  await page.evaluate((vn: string) => {
    (document.querySelector(`[name="${vn}"]`) as HTMLElement)?.click();
  }, viewName);
  await page.waitForTimeout(200);
  await page.evaluate((vn: string) => {
    const item = (document.querySelector(`[name="${vn}"]`) as HTMLElement)?.closest('.property-grid-item');
    (item?.querySelector('input.property-grid-item-editor-textbox') as HTMLInputElement)?.focus();
  }, viewName);
  await page.keyboard.press('Control+A');
  await page.keyboard.press('Delete');
  await page.keyboard.press('Enter');
  await page.waitForTimeout(300);
}

test('Density Plot — Bin Selection, Zoom, Bin To Range, Axis Configuration', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text()))
      consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;

  const readProp = (prop: string) => page.evaluate((p: string) => {
    const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
    return d?.props[p];
  }, prop);

  // Settle-gated canvas ink (total non-white), used to wait out render tails.
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

  // Summed per-color pixel delta between the pre-action settled frame and the
  // post-action settled frame — catches bin recomputes that recolor in place.
  const measureColorDelta = async (action: () => Promise<void>) => {
    await settledPx();
    await v.snapshotCanvasColors(page, 'Density plot');
    await action();
    await settledPx();
    return (await v.diffCanvasColors(page, 'Density plot')).deltaPx;
  };

  const zoomInOneStep = () => page.evaluate(() => {
    const canvas = document.querySelector('[name="viewer-Density-plot"] canvas[name="canvas"]') as HTMLElement;
    const r = canvas.getBoundingClientRect();
    canvas.dispatchEvent(new WheelEvent('wheel', {
      bubbles: true, cancelable: true,
      clientX: r.left + r.width / 2, clientY: r.top + r.height / 2, deltaY: -300,
    }));
  });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'density-plot', 'Density-plot');

  // Setup: X=AGE, Y=HEIGHT (UI selector path owned by the p0 scenario), defaults
  // bins 50 / hexagon confirmed.
  const setup = await page.evaluate(async () => {
    const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
    d.props.xColumnName = 'AGE';
    d.props.yColumnName = 'HEIGHT';
    d.props.bins = 50;
    d.props.binShape = 'hexagon';
    await new Promise((r) => setTimeout(r, 800));
    return {x: d.props.xColumnName, y: d.props.yColumnName, bins: d.props.bins, shape: d.props.binShape};
  });
  expect(setup).toEqual({x: 'AGE', y: 'HEIGHT', bins: 50, shape: 'hexagon'});

  await softStep('Scenario 1 — bin click selection, tooltip, and the Esc round-trip', async () => {
    const box = await page.locator('[name="viewer-Density-plot"] canvas[name="canvas"]').boundingBox();
    expect(box).toBeTruthy();
    const cx = box!.x + box!.width / 2;
    const cy = box!.y + box!.height / 2;

    const baseline = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(baseline).toBe(0);

    // Real (trusted) click at the canvas center — a populated bin at AGE/HEIGHT.
    await page.mouse.click(cx, cy);
    await page.waitForTimeout(500);
    const selected = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    console.log(`DensityPlot bin-click selected=${selected}`);
    expect(selected).toBeGreaterThan(0);

    // Hover the same bin → document-level tooltip carrying both bin ranges + row count.
    await page.mouse.move(cx + 2, cy + 2);
    await page.mouse.move(cx, cy);
    let tooltip = '';
    for (let i = 0; i < 20; i++) {
      await page.waitForTimeout(200);
      tooltip = await page.evaluate(() =>
        (document.querySelector('.d4-tooltip') as HTMLElement)?.innerText ?? '');
      if (/rows/.test(tooltip)) break;
    }
    console.log(`DensityPlot tooltip="${tooltip.replace(/\n/g, ' | ')}"`);
    expect(tooltip).toMatch(/AGE/);
    expect(tooltip).toMatch(/HEIGHT/);
    expect(tooltip).toMatch(/(\d+)\s+rows/);

    // Esc with focus inside the viewer clears the selection (round-trip).
    await page.mouse.move(cx, cy);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(400);
    const cleared = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(cleared).toBe(0);
  });

  await softStep('Scenario 2 — mouse-wheel zoom and Reset View from the context menu', async () => {
    const errBefore = errCount();
    const basePx = await settledPx();
    await zoomInOneStep();
    const zoomPx = await settledPx();
    console.log(`DensityPlot zoom px: basePx=${basePx} zoomPx=${zoomPx} delta=${Math.abs(zoomPx - basePx)}`);
    // The wheel zoom repaints the bins over the zoomed viewport — a clear delta;
    // the floor keeps a wide margin below the observed zoom delta.
    expect(basePx).toBeGreaterThan(0);
    expect(Math.abs(zoomPx - basePx)).toBeGreaterThan(50000);

    // Reset View via the canvas context menu (a text menu item; no name= attr).
    await page.evaluate(() => {
      const canvas = document.querySelector('[name="viewer-Density-plot"] canvas[name="canvas"]') as HTMLElement;
      const r = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: r.left + r.width / 2, clientY: r.top + r.height / 2,
      }));
    });
    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('[role="menuitem"]'))
        .some((el) => el.textContent?.trim() === 'Reset View'), {timeout: 3000});
    await page.evaluate(() => {
      const item = Array.from(document.querySelectorAll('[role="menuitem"]'))
        .find((el) => el.textContent?.trim() === 'Reset View') as HTMLElement;
      item?.click();
    });
    const resetPx = await settledPx();
    console.log(`DensityPlot reset px: resetPx=${resetPx} |reset-base|=${Math.abs(resetPx - basePx)} |zoom-base|=${Math.abs(zoomPx - basePx)}`);
    // Reset is not pixel-exact (bins are recomputed) — assert a DIRECTIONAL
    // restore: the reset frame is closer to the pre-zoom baseline than the
    // zoomed frame was. Double-click is deliberately NOT used (it selects a bin).
    expect(Math.abs(resetPx - basePx)).toBeLessThan(Math.abs(zoomPx - basePx));
    expect(errCount()).toBe(errBefore);
  });

  // Open the settings panel for the remaining scenarios (gear is two hops up).
  await page.evaluate(() => {
    const viewer = document.querySelector('[name="viewer-Density-plot"]') as HTMLElement;
    const gear = viewer?.parentElement?.parentElement
      ?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear?.click();
  });
  await page.waitForTimeout(1000);
  await expandCategory(page, 'x');
  await expandCategory(page, 'y');
  await expandCategory(page, 'misc');

  await softStep('Scenario 3 — Bin To Range toggled while zoomed', async () => {
    const errBefore = errCount();
    await zoomInOneStep();
    await settledPx();
    const toggleBinToRange = () => page.evaluate(() => {
      (document.querySelector('[name="prop-bin-to-range"] input[type="checkbox"]') as HTMLInputElement)?.click();
    });
    const onDelta = await measureColorDelta(async () => { await toggleBinToRange(); await page.waitForTimeout(400); });
    const onProp = await readProp('binToRange');
    const offDelta = await measureColorDelta(async () => { await toggleBinToRange(); await page.waitForTimeout(400); });
    const offProp = await readProp('binToRange');
    console.log(`DensityPlot binToRange color delta: onDelta=${onDelta} offDelta=${offDelta}`);
    expect(onProp).toBe(true);
    expect(offProp).toBe(false);
    // With the plot zoomed, the bins recompute against the full data range and
    // back — a real repaint in each direction; the floor keeps a wide margin
    // below the observed color-delta.
    expect(onDelta).toBeGreaterThan(50000);
    expect(offDelta).toBeGreaterThan(50000);
    expect(errCount()).toBe(errBefore);
    // Return to the full range for the next scenario.
    await page.evaluate(() => {
      const canvas = document.querySelector('[name="viewer-Density-plot"] canvas[name="canvas"]') as HTMLElement;
      const r = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: r.left + r.width / 2, clientY: r.top + r.height / 2,
      }));
    });
    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('[role="menuitem"]'))
        .some((el) => el.textContent?.trim() === 'Reset View'), {timeout: 3000}).catch(() => {});
    await page.evaluate(() => {
      const item = Array.from(document.querySelectorAll('[role="menuitem"]'))
        .find((el) => el.textContent?.trim() === 'Reset View') as HTMLElement;
      item?.click();
    });
    await settledPx();
  });

  await softStep('Scenario 4 — explicit axis bounds with Bin To Range', async () => {
    const errBefore = errCount();
    const basePx = await settledPx();
    await setBound(page, 'prop-view-x-min', 30);
    await setBound(page, 'prop-view-x-max', 60);
    await setBound(page, 'prop-view-y-min', 150);
    await setBound(page, 'prop-view-y-max', 190);
    const bounds = await page.evaluate(() => {
      const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      return {xMin: d.props.xMin, xMax: d.props.xMax, yMin: d.props.yMin, yMax: d.props.yMax};
    });
    expect(bounds).toEqual({xMin: 30, xMax: 60, yMin: 150, yMax: 190});
    const boundedPx = await settledPx();

    // Bin To Range on then off under explicit bounds — held as an error-free
    // floor with the prop round-trip confirmed (the delta can be modest here).
    const toggleBinToRange = () => page.evaluate(() => {
      (document.querySelector('[name="prop-bin-to-range"] input[type="checkbox"]') as HTMLInputElement)?.click();
    });
    await toggleBinToRange();
    await page.waitForTimeout(400);
    const onProp = await readProp('binToRange');
    await toggleBinToRange();
    await page.waitForTimeout(400);
    const offProp = await readProp('binToRange');
    expect(onProp).toBe(true);
    expect(offProp).toBe(false);

    // Clear all four bounds — the render restores directionally toward the
    // full-range baseline.
    await clearBound(page, 'prop-view-x-min');
    await clearBound(page, 'prop-view-x-max');
    await clearBound(page, 'prop-view-y-min');
    await clearBound(page, 'prop-view-y-max');
    const cleared = await page.evaluate(() => {
      const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      // Normalize the "no bound" reading — the prop surfaces as null or undefined.
      const nn = (x: any) => (x === null || x === undefined);
      return {xMin: nn(d.props.xMin), xMax: nn(d.props.xMax), yMin: nn(d.props.yMin), yMax: nn(d.props.yMax)};
    });
    const clearedPx = await settledPx();
    console.log(`DensityPlot bounds px: basePx=${basePx} boundedPx=${boundedPx} clearedPx=${clearedPx}`);
    expect(cleared).toEqual({xMin: true, xMax: true, yMin: true, yMax: true});
    // Clearing restores toward baseline: the cleared frame is closer to the
    // full-range baseline than the bounded frame was.
    expect(Math.abs(clearedPx - basePx)).toBeLessThan(Math.abs(boundedPx - basePx));
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 5 — axis invert and logarithmic scale, including the non-positive edge', async () => {
    const errBefore = errCount();
    // Invert Y axis + logarithmic on HEIGHT (all-positive), then revert.
    await page.evaluate(() => {
      (document.querySelector('[name="prop-invert-y-axis"] input[type="checkbox"]') as HTMLInputElement)?.click();
    });
    await page.waitForTimeout(300);
    expect(await readProp('invertYAxis')).toBe(true);
    await setPropSelect(page, 'prop-y-axis-type', 'logarithmic');
    expect(await readProp('yAxisType')).toBe('logarithmic');
    await setPropSelect(page, 'prop-y-axis-type', 'linear');
    await page.evaluate(() => {
      (document.querySelector('[name="prop-invert-y-axis"] input[type="checkbox"]') as HTMLInputElement)?.click();
    });
    await page.waitForTimeout(300);
    const reverted = await page.evaluate(() => {
      const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      return {invertY: d.props.invertYAxis, yType: d.props.yAxisType};
    });
    expect(reverted).toEqual({invertY: false, yType: 'linear'});
    // The invert + log round-trip on HEIGHT raises no page or console errors.
    expect(errCount()).toBe(errBefore);

    // Non-positive edge: a computed column whose minimum is at or below zero.
    const errBeforeEdge = errCount();
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      if (!df.columns.names().includes('LOG_EDGE'))
        df.columns.addNewCalculated('LOG_EDGE', '${HEIGHT} - 200');
      const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      d.props.yColumnName = 'LOG_EDGE';
      await new Promise((r) => setTimeout(r, 600));
    });
    await setPropSelect(page, 'prop-y-axis-type', 'logarithmic');
    await page.waitForTimeout(600);
    const alive = await page.evaluate(() => {
      const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      return {rootAttached: document.body.contains(d.root), yType: d.props.yAxisType};
    });
    console.log(`DensityPlot log-edge alive=${JSON.stringify(alive)} errDelta=${errCount() - errBeforeEdge}`);
    // The viewer stays alive and renders an on-canvas warning (canvas-only, not
    // asserted) instead of throwing — no page or console errors.
    expect(alive.rootAttached).toBe(true);
    expect(errCount()).toBe(errBeforeEdge);

    // Restore Y to HEIGHT / linear and drop the computed column.
    await setPropSelect(page, 'prop-y-axis-type', 'linear');
    await page.evaluate(async () => {
      const d = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Density plot') as any;
      d.props.yColumnName = 'HEIGHT';
      const df = grok.shell.tv.dataFrame;
      const edge = df.columns.names().includes('LOG_EDGE') ? df.col('LOG_EDGE') : null;
      if (edge) df.columns.remove('LOG_EDGE');
      await new Promise((r) => setTimeout(r, 400));
    });
  });

  await page.evaluate(() => grok.shell.closeAll());

  v.finishSpec();
});
