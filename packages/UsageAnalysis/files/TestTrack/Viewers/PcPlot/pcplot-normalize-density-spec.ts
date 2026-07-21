/* ---
realizes: [pcplot.cp.normalize-and-density]
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

  await softStep('Scenario 1 — switch vertical scale global then back', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find((vw: any) => vw.type === 'PC Plot')!;
      pc.props.normalizeEachColumn = false;
      await new Promise((r) => setTimeout(r, 300));
      pc.props.normalizeEachColumn = true;
      await new Promise((r) => setTimeout(r, 300));
    });
    expect(await readRootInDom()).toBe(true);
    expect(pageErrors.length + consoleErrors.length).toBe(errBefore);
  });

  await softStep('Scenario 2 — enable density, cycle circles/box/violin styles, drive per-part toggles', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
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
    const alive = await page.evaluate(async () => {
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
