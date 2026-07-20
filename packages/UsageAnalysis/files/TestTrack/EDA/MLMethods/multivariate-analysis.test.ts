import { test, expect } from './helpers';
import {
  clickDialogPrimary, clickTopMenuLeaf, currentViewerTypes, openDemoCsv,
  resetShell, waitForDialog,
} from './helpers';

// Test Track scenario: EDA/multivariate-analysis.md
// 1. Open cars.csv from Demo files.
// 2. Top Menu > ML > Analyze > Multivariate Analysis... — viewers display.
// 3. Check interactivity: Grid vs Observed/Predicted/Scores scatterplots,
//    Loadings scatterplot vs Regression coefficients bar chart.
//
// Step 3 has no scriptable success criterion in the scenario text. We assert the
// structural prerequisite for interactivity (the full viewer set is attached to the
// same TableView, so the standard selection bus wires them up automatically) and
// leave per-pixel hover propagation to manual verification.

test.describe.serial('EDA / Multivariate Analysis', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('MVA on cars.csv produces Grid + 3 Scatter plots + 2 Bar charts', async ({ page }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, 'cars.csv');

    await clickTopMenuLeaf(page, 'div-ML---Analyze---Multivariate-Analysis...');
    await waitForDialog(page, 'Multivariate Analysis');

    // Defaults (Predict=price, Using=15 numeric cols excluding price, Components=2)
    // produce a valid run on cars.csv — explicitly selecting "All" in Using puts
    // `price` in both Predict and Using and disables RUN. See run notes in
    // public/packages/UsageAnalysis/files/TestTrack/EDA/multivariate-analysis-run.md.
    await clickDialogPrimary(page, ['Run', 'RUN', 'OK']);

    await expect(page.locator('.d4-dialog .d4-dialog-title', { hasText: /Multivariate Analysis/i }))
      .toHaveCount(0, { timeout: 15_000 });

    // MVA renders six analytic viewers: Grid + Observed/Predicted + Scores +
    // Loadings + Regression coefficients + Explained variance. Wait in-browser until the
    // expected analytic viewers (3 scatter + 2 bar, plus the table grid) appear —
    // `waitForFunction` polls inside the page and rejects only once on timeout, so
    // intermediate "not ready yet" attempts are not logged as report errors.
    await page.waitForFunction(() => {
      const tv = (window as any).grok?.shell?.tv;
      if (!tv) return false;
      const lc = Array.from(tv.viewers).map((v: any) => String(v?.type ?? '').toLowerCase());
      return lc.filter((t) => t.includes('scatter')).length >= 3
        && lc.filter((t) => t.includes('bar')).length >= 2
        && lc.some((t) => t === 'grid');
    }, undefined, { timeout: 120_000 });

    // Final structural assertion — runs once, after the wait already confirmed readiness.
    const lc = (await currentViewerTypes(page)).map((t) => t.toLowerCase());
    expect(lc.filter((t) => t.includes('scatter')).length, 'scatter viewers').toBeGreaterThanOrEqual(3);
    expect(lc.filter((t) => t.includes('bar')).length, 'bar viewers').toBeGreaterThanOrEqual(2);
    expect(lc.some((t) => t === 'grid'), 'grid viewer').toBe(true);

    // Interactivity precondition: all viewers share the same DataFrame, so the
    // standard selection bus is in place. Validate by reading the DF row count
    // from the active TableView — it must be the original 30 cars rows.
    const rowCount = await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.rowCount);
    expect(rowCount).toBe(30);
  });
});
