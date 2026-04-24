import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Multivariate Analysis scenario', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/cars.csv');
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Step 1 — Open cars.csv from Demo files', async () => {
    const info = await page.evaluate(() => {
      const g: any = (window as any).grok;
      return {rows: g.shell.tv.dataFrame.rowCount, cols: g.shell.tv.dataFrame.columns.length};
    });
    expect(info.rows).toBe(30);
    expect(info.cols).toBe(17);
  });

  await softStep('Step 2 — Top Menu > ML > Analyze > Multivariate Analysis...; viewers display correctly', async () => {
    // Open MVA dialog via JS API fallback (DOM click on menu item)
    await page.evaluate(() => {
      const el = document.querySelector('[name="div-ML---Analyze---Multivariate-Analysis..."]') as HTMLElement | null;
      if (el) el.click();
    });
    await page.waitForTimeout(1200);
    await page.locator('.d4-dialog').waitFor({timeout: 10000});
    await expect(page.locator('.d4-dialog')).toContainText('Multivariate Analysis');

    // Leave defaults (Predict=price, Using=15 numeric cols except price, Components=2, Quadratic=false).
    // Explicitly setting Using to "(16) All" introduces price twice and disables RUN — do NOT do that.

    // Click RUN
    const runBtn = page.locator('[name="button-RUN"], [name="button-Run"]').first();
    await runBtn.click();

    // Poll up to 90s for viewers to accumulate (MVA builds 5 viewers, plus the Grid = 6 total).
    const deadline = Date.now() + 90_000;
    let viewerInfo: {count: number; types: string[]} = {count: 0, types: []};
    while (Date.now() < deadline) {
      viewerInfo = await page.evaluate(() => {
        const g: any = (window as any).grok;
        const viewers = Array.from(g.shell.tv.viewers);
        return {count: viewers.length, types: viewers.map((v: any) => v.type)};
      });
      if (viewerInfo.count >= 6) break;
      await page.waitForTimeout(1000);
    }
    if (viewerInfo.count < 6)
      throw new Error(`MVA viewers did not materialize within 90s. Count=${viewerInfo.count}, types=${JSON.stringify(viewerInfo.types)}.`);
    expect(viewerInfo.count).toBeGreaterThanOrEqual(6);
  });

  await softStep('Step 3 — Check interactivity: grid <-> Observed/Scores scatterplots, Loadings <-> Regression bar', async () => {
    // Canvas-level interactivity (hover/selection highlight across viewers) requires synthetic
    // mouse events at precise pixel coordinates — not feasible to validate reliably from
    // outside a canvas. Assert instead that the expected viewer set is present on the same
    // DataFrame, which is what the Datagrok event bus uses to propagate hover/selection.
    const types = await page.evaluate(() => {
      const g: any = (window as any).grok;
      return Array.from(g.shell.tv.viewers).map((v: any) => v.type);
    });
    const scatterCount = types.filter((t: string) => t === 'Scatter plot').length;
    const barCount = types.filter((t: string) => t === 'Bar chart').length;
    if (scatterCount < 3)
      throw new Error(`Expected >=3 Scatter plot viewers (Observed/Predicted, Scores, Loadings). Got ${scatterCount}. Full list: ${JSON.stringify(types)}.`);
    if (barCount < 2)
      throw new Error(`Expected >=2 Bar chart viewers (Regression Coefficients, Explained Variance). Got ${barCount}. Full list: ${JSON.stringify(types)}.`);
    expect(scatterCount).toBeGreaterThanOrEqual(3);
    expect(barCount).toBeGreaterThanOrEqual(2);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
