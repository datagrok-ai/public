/* ---
realizes: [pcplot.cp.normalize-and-density, pcplot.int.normalize-density-overlay]
--- */
import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
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
    if (m.type() === 'error' && !isLocalBootNoise(m.text()))
      consoleErrors.push(m.text());
  });

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon)
      (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 15000});

  await v.installEventWaits(page);

  await v.setViewerProps(page, 'PC Plot', [{set: {columnNames: ['AGE', 'HEIGHT', 'WEIGHT']}}], 500);

  const readRootInDom = () => page.evaluate(() => {
    const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
    return document.body.contains(pc.root);
  });

  const settledPx = async () => {
    await v.waitForCanvasQuiet(page, 'PC Plot', {timeoutMs: 1500, optional: true});
    return (await v.countCanvasPixels(page, 'PC Plot')).total;
  };

  const setAllLines = async (on: boolean) => {
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    await v.setViewerProps(page, 'PC Plot', [{set: {showAllLines: on}}], 400);
  };

  let floorPx = -1;

  await softStep('Scenario 1 — switch vertical scale global then back', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const setNormalize = (on: boolean) =>
      v.setViewerProps(page, 'PC Plot', [{set: {normalizeEachColumn: on}}], 300);
    const normalizedPx = await settledPx();
    await setNormalize(false);
    const globalPx = await settledPx();
    await setNormalize(true);
    const backPx = await settledPx();

    console.log(`Normalize scale px: normalizedPx=${normalizedPx} globalPx=${globalPx} backPx=${backPx}`);

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
    await v.setViewerProps(page, 'PC Plot', [{set: {showDensity: true, densityStyle: 'box plot'}}], 400);
    const densityPx = await settledPx();
    console.log(`Density px: floorPx=${floorPx} densityPx=${densityPx}`);

    expect(floorPx).toBeGreaterThanOrEqual(0);
    expect(densityPx - floorPx).toBeGreaterThan(1000);
    await setAllLines(true);
    await v.setViewerProps(page, 'PC Plot', [
      {set: {showDensity: true}, wait: 300},
      {set: {densityStyle: 'circles'}, wait: 300},
      {set: {densityStyle: 'box plot'}, wait: 300},
      {set: {densityStyle: 'violin plot'}, wait: 300},
      {set: {densityStyle: 'box plot'}, wait: 200},
    ]);

    await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.showMedian = false; pc.props.showMedian = true;
      pc.props.showInterquartileRange = false; pc.props.showInterquartileRange = true;
      pc.props.showMeanCross = false; pc.props.showMeanCross = true;
      pc.props.showUpperDash = false; pc.props.showUpperDash = true;
      pc.props.showLowerDash = false; pc.props.showLowerDash = true;
      pc.props.showCircles = true;
      pc.props.bins = 200;
    });
    await v.waitForViewerRendered(page, 'PC Plot', 300);
    await v.setViewerProps(page, 'PC Plot', [{set: {bins: 100, showDensity: false}}], 300);
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Scenario 3 — density recalculates on normalization double-toggle and AGE log scale (github-1546)', async () => {

    const errBefore = pageErrors.length + consoleErrors.length;
    await v.setViewerProps(page, 'PC Plot', [
      {set: {showDensity: true, densityStyle: 'box plot'}, wait: 300},
      {set: {normalizeEachColumn: false}, wait: 400},
      {set: {normalizeEachColumn: true}, wait: 400},
      {set: {normalizeEachColumn: false}, wait: 400},
      {set: {normalizeEachColumn: true}, wait: 400},
    ]);

    await setAllLines(false);
    const staleGuardPx = await settledPx();
    console.log(`Stale-guard px: staleGuardPx=${staleGuardPx} floorPx=${floorPx}`);
    expect(floorPx).toBeGreaterThanOrEqual(0);
    expect(staleGuardPx - floorPx).toBeGreaterThan(1000);
    await setAllLines(true);
    await v.setViewerProps(page, 'PC Plot', [{set: {logColumnsColumnNames: ['AGE']}}], 400);

    const alive = await page.evaluate(() => document.body.isConnected);
    await v.setViewerProps(page, 'PC Plot', [
      {set: {logColumnsColumnNames: []}, wait: 300},
      {set: {showDensity: false}, wait: 300},
    ]);
    expect(alive).toBe(true);
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  v.finishSpec();
});
