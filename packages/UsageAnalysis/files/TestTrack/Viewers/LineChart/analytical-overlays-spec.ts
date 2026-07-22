/* ---
realizes: [linechart.cp.analytical-overlays]
realizes_atlas: [linechart.cp.analytical-overlays]
--- */
import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

const pageErrors: string[] = [];
const consoleErrors: string[] = [];

function realErrors(): string[] {
  return [...pageErrors, ...consoleErrors];
}

async function setProps(page: Page, props: Record<string, any>) {
  await page.evaluate((p) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    for (const [k, val] of Object.entries(p)) (lc.props as any)[k] = val;
  }, props);
  await page.waitForTimeout(500);
}

async function formulaLinesCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    try {
      const parsed = JSON.parse(lc.props.formulaLines || '[]');
      return Array.isArray(parsed) ? parsed.length : -1;
    } catch (e) {
      return -1;
    }
  });
}

test('Line Chart — Analytical Overlays', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  pageErrors.length = 0;
  consoleErrors.length = 0;

  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 5000));
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60000});

  await page.locator('[name="icon-line-chart"]').click();
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 15000});

  await setProps(page, {xColumnName: 'CAST Idea ID', yColumnNames: ['Chemical Space X']});

  // The overlays are canvas-rendered and grok.shell.warnings is undefined on this
  // build, so page and console errors are the floor for the steps below.
  await softStep('S1 steps 1-3: enable regression line, no-error floor', async () => {
    const before = realErrors().length;
    await setProps(page, {showRegressionLine: true});
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1 steps 4-5: enable rolling average overlay, no split', async () => {
    const before = realErrors().length;
    await setProps(page, {showMovingAverageLine: true});
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1 steps 6-7: enable standard-deviation overlay, no split', async () => {
    const before = realErrors().length;
    await setProps(page, {showMovingAverageDeviation: true});
    expect(realErrors().length).toBe(before);
  });

  // GROK-20218: with a split, the rolling-average and std-dev overlays render per
  // category.
  await softStep('S1 steps 8-9: apply split, overlays render per category', async () => {
    const before = realErrors().length;
    await setProps(page, {splitColumnNames: ['Stereo Category']});
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1 steps 10-11: revert overlays and split', async () => {
    const before = realErrors().length;
    await setProps(page, {
      showRegressionLine: false,
      showMovingAverageLine: false,
      showMovingAverageDeviation: false,
      splitColumnNames: [],
    });
    expect(realErrors().length).toBe(before);
  });

  const formulaLinesSpec = JSON.stringify([
    {type: 'line', formula: '${Chemical Space X} = 500', title: 'const-line', color: '#FF0000'},
    {type: 'band', formula: '${Chemical Space X} in(400, 600)', title: 'const-band', color: '#00FF00'},
  ]);
  // The line and the band are canvas-rendered; only the formulaLines config itself
  // can be read back.
  await softStep('S2 steps 1-5: add a formula line and a formula band', async () => {
    const before = realErrors().length;
    await setProps(page, {formulaLines: formulaLinesSpec});
    expect(await formulaLinesCount(page)).toBe(2);
    expect(realErrors().length).toBe(before);
  });

  // GROK-19943, GROK-19949: the formula-lines config must survive a layout save and
  // reapply.
  await softStep('S2 steps 6-8: save layout, clear, reapply — formula lines restored', async () => {
    const before = realErrors().length;
    const countBeforeSave = await formulaLinesCount(page);
    expect(countBeforeSave).toBe(2);

    const layoutId = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1200));
      return layout.id;
    });

    // Clear the formula lines to prove reapply restores them.
    await setProps(page, {formulaLines: ''});
    expect(await formulaLinesCount(page)).toBe(0);

    await page.evaluate(async (id) => {
      const tv = grok.shell.tv;
      const saved = await grok.dapi.layouts.find(id);
      tv.loadLayout(saved);
      await new Promise((r) => setTimeout(r, 3000));
    }, layoutId);

    expect(await formulaLinesCount(page)).toBe(2);
    expect(realErrors().length).toBe(before);

    await setProps(page, {formulaLines: ''});
    await page.evaluate(async (id) => {
      try {
        const saved = await grok.dapi.layouts.find(id);
        if (saved) await grok.dapi.layouts.delete(saved);
      } catch (e) { /* best-effort cleanup */ }
    }, layoutId);
  });

  if (stepErrors.length > 0)
    throw new Error(`Line Chart analytical-overlays failures:\n${stepErrors.join('\n')}`);
});
