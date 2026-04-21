import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
});

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e.message ?? String(e)}); }
}

test('ANOVA scenario', async ({page}) => {
  test.setTimeout(300_000);

  const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
  const login = process.env.DATAGROK_LOGIN ?? 'admin';
  const password = process.env.DATAGROK_PASSWORD ?? 'admin';

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 60000});

  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
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

  await softStep('Step 1 — Open demog.csv from Demo Files', async () => {
    const info = await page.evaluate(() => {
      const g: any = (window as any).grok;
      return {rows: g.shell.tv.dataFrame.rowCount, cols: g.shell.tv.dataFrame.columns.length};
    });
    expect(info.rows).toBe(5850);
    expect(info.cols).toBeGreaterThanOrEqual(11);
  });

  await softStep('Step 2 — Top Menu > ML > Analyze > ANOVA...', async () => {
    await page.evaluate(() => {
      const el = document.querySelector('[name="div-ML---Analyze---ANOVA..."]') as HTMLElement | null;
      if (el) el.click();
    });
    await page.waitForTimeout(1200);
    await page.locator('.d4-dialog').waitFor({timeout: 10000});
    await expect(page.locator('.d4-dialog')).toContainText('ANOVA');
    await expect(page.locator('[name="input-host-Category"]')).toBeVisible();
    await expect(page.locator('[name="input-host-Feature"]')).toBeVisible();
    await expect(page.locator('[name="input-host-Alpha"]')).toBeVisible();
  });

  await softStep('Step 3 — Click Run; expect Box plot + Analysis + F-test tabs', async () => {
    // Dialog button is named "Run" (mixed case) — NOT "button-RUN" or "button-OK".
    await page.locator('.d4-dialog [name="button-Run"]').last().click();

    // Poll up to 30s for the Box plot viewer + tab host with Analysis/F-test tabs.
    const deadline = Date.now() + 30_000;
    let hasBoxPlot = false;
    let hasAnalysisTab = false;
    let hasFTestTab = false;
    while (Date.now() < deadline) {
      const snap = await page.evaluate(() => {
        const g: any = (window as any).grok;
        const viewers = Array.from(g.shell.tv?.viewers ?? []);
        const boxPlot = viewers.some((v: any) => v.type === 'Box plot');
        const tabHost = document.querySelector('.d4-tab-host');
        const tabLabels = tabHost
          ? Array.from(tabHost.querySelectorAll('.d4-tab-header'))
              .map((el) => (el as HTMLElement).textContent?.trim() ?? '')
          : [];
        return {
          boxPlot,
          analysis: tabLabels.includes('Analysis'),
          ftest: tabLabels.includes('F-test'),
        };
      });
      hasBoxPlot = snap.boxPlot;
      hasAnalysisTab = snap.analysis;
      hasFTestTab = snap.ftest;
      if (hasBoxPlot && hasAnalysisTab && hasFTestTab) break;
      await page.waitForTimeout(500);
    }

    if (!hasBoxPlot)
      throw new Error('Box plot viewer did not appear within 30s.');
    if (!hasAnalysisTab)
      throw new Error('Analysis tab did not appear within 30s.');
    if (!hasFTestTab)
      throw new Error('F-test tab did not appear within 30s.');

    expect(hasBoxPlot).toBe(true);
    expect(hasAnalysisTab).toBe(true);
    expect(hasFTestTab).toBe(true);

    // Final assertion via Playwright locator too.
    await expect(page.locator('[name="viewer-Box-plot"]')).toBeVisible({timeout: 5000});
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
