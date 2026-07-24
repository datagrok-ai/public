/* ---
realizes: [pcplot.cp.normalize-and-density, pcplot-normalize-density-overlay]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('PC Plot — Normalization and Density Overlay', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error')
      consoleErrors.push(m.text());
  });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon)
      (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 15000});

  await page.evaluate(async () => {
    const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
    pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
    await new Promise((r) => setTimeout(r, 500));
  });

  const readRootInDom = () => page.evaluate(() => {
    const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
    return document.body.contains(pc.root);
  });

  // Settle-gated canvas ink count: repeated reads until two consecutive counts
  // agree, so a delta between two measurements is the setter's effect, not a
  // render tail.
  const settledPx = async () => {
    let prev = (await v.countCanvasPixels(page, 'PC Plot')).total;
    let cur = prev;
    for (let i = 0; i < 5; i++) {
      await page.waitForTimeout(300);
      cur = (await v.countCanvasPixels(page, 'PC Plot')).total;
      if (Math.abs(cur - prev) < 200) break;
      prev = cur;
    }
    return cur;
  };

  // Sparse-canvas basis for the density-ink measurements: with all polylines
  // hidden and the selection empty, only axes/labels are painted, so density
  // shapes become the dominant ink instead of paint-over on a full canvas.
  const setAllLines = (on: boolean) => page.evaluate(async (val) => {
    const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
    grok.shell.tv.dataFrame.selection.setAll(false);
    pc.props.showAllLines = val;
    await new Promise((r) => setTimeout(r, 400));
  }, on);

  // Hidden-lines ink floor measured in Scenario 2 and reused by Scenario 3's
  // anti-stale check.
  let floorPx = -1;

  await softStep('Scenario 1 — switch vertical scale global then back', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const setNormalize = (on: boolean) => page.evaluate(async (val) => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.normalizeEachColumn = val;
      await new Promise((r) => setTimeout(r, 300));
    }, on);
    const normalizedPx = await settledPx();
    await setNormalize(false);
    const globalPx = await settledPx();
    await setNormalize(true);
    const backPx = await settledPx();
    // Keep the measured ink values visible on green runs so the fixed
    // thresholds below can be audited against live numbers.
    console.log(`Normalize scale px: normalizedPx=${normalizedPx} globalPx=${globalPx} backPx=${backPx}`);
    // Repaint-delta: the shared global scale redistributes the polylines, so
    // the ink count must move by a real margin (the direction depends on the
    // data), and the round-trip back to Normalized must land near the original.
    expect(normalizedPx).toBeGreaterThanOrEqual(0);
    expect(globalPx).toBeGreaterThanOrEqual(0);
    expect(Math.abs(globalPx - normalizedPx)).toBeGreaterThan(500);
    expect(Math.abs(backPx - normalizedPx)).toBeLessThan(2000);
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Scenario 2 — enable density, cycle circles/box/violin styles, drive per-part toggles', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    expect(await page.evaluate(() =>
      grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!.props.densityStyle)).toBe('circles');
    await setAllLines(false);
    floorPx = await settledPx();
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.showDensity = true;
      pc.props.densityStyle = 'box plot';
      await new Promise((r) => setTimeout(r, 400));
    });
    const densityPx = await settledPx();
    console.log(`Density px: floorPx=${floorPx} densityPx=${densityPx}`);
    // Repaint-delta on the sparse basis, measured on the BOX PLOT style: with
    // the polylines hidden, enabling the overlay must ADD ink over the
    // axes-only floor by a real margin. The circles overlay hugs the axis line
    // and is not pixel-distinguishable from it, so circles stays covered by
    // the no-error cycling floor below instead of a delta assert.
    expect(floorPx).toBeGreaterThanOrEqual(0);
    expect(densityPx - floorPx).toBeGreaterThan(1000);
    await setAllLines(true);
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.showDensity = true;
      await new Promise((r) => setTimeout(r, 300));
      for (const s of ['circles', 'box plot', 'violin plot']) {
        pc.props.densityStyle = s;
        await new Promise((r) => setTimeout(r, 300));
      }
      // Per-part box/violin component toggles named in the .md — Show Median,
      // Interquartile Range, Mean Cross, upper/lower whisker dash, Show Circles
      // and the Bins change. Each redraws to canvas only, so they share the same
      // no-error + DOM-presence floor rather than a per-toggle assertion.
      pc.props.densityStyle = 'box plot';
      await new Promise((r) => setTimeout(r, 200));
      pc.props.showMedian = false; pc.props.showMedian = true;
      pc.props.showInterquartileRange = false; pc.props.showInterquartileRange = true;
      pc.props.showMeanCross = false; pc.props.showMeanCross = true;
      pc.props.showUpperDash = false; pc.props.showUpperDash = true;
      pc.props.showLowerDash = false; pc.props.showLowerDash = true;
      pc.props.showCircles = true;
      pc.props.bins = 200;
      await new Promise((r) => setTimeout(r, 300));
      pc.props.bins = 100;
      pc.props.showDensity = false;
      await new Promise((r) => setTimeout(r, 300));
    });
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Scenario 3 — density recalculates on normalization double-toggle and AGE log scale (github-1546)', async () => {
    // github-1546: the density overlay must recalculate without throwing or freezing
    // after a normalization double-toggle and a log-scale switch.
    const errBefore = pageErrors.length + consoleErrors.length;
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.showDensity = true;
      pc.props.densityStyle = 'box plot';
      await new Promise((r) => setTimeout(r, 300));
      pc.props.normalizeEachColumn = false;
      await new Promise((r) => setTimeout(r, 400));
      pc.props.normalizeEachColumn = true;
      await new Promise((r) => setTimeout(r, 400));
      pc.props.normalizeEachColumn = false;
      await new Promise((r) => setTimeout(r, 400));
      pc.props.normalizeEachColumn = true;
      await new Promise((r) => setTimeout(r, 400));
    });
    // Anti-stale on the sparse basis: hide the polylines so the surviving
    // density ink is measured against the Scenario 2 axes-only floor — the
    // overlay is still painted after the double-toggle, not stale-empty.
    await setAllLines(false);
    const staleGuardPx = await settledPx();
    console.log(`Stale-guard px: staleGuardPx=${staleGuardPx} floorPx=${floorPx}`);
    expect(floorPx).toBeGreaterThanOrEqual(0);
    expect(staleGuardPx - floorPx).toBeGreaterThan(1000);
    await setAllLines(true);
    const alive = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.logColumnsColumnNames = ['AGE'];
      await new Promise((r) => setTimeout(r, 400));
      const stillResponsive = await new Promise((r) => setTimeout(() => r(true), 10));
      pc.props.logColumnsColumnNames = [];
      await new Promise((r) => setTimeout(r, 300));
      pc.props.showDensity = false;
      await new Promise((r) => setTimeout(r, 300));
      return stillResponsive;
    });
    expect(alive).toBe(true);
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  v.finishSpec();
});
