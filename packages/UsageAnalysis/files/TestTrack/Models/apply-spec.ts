import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Apply predictive model', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup: open demog.csv, apply selenium class / Tabs mode, wait for grid
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    const g: any = (window as any).grok;
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
    const tv = g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    (window as any).__initialColCount = df.columns.length;
    (window as any).__initialColNames = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await softStep('Open demog.csv', async () => {
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return {rows: df?.rowCount ?? 0, cols: df?.columns?.length ?? 0};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.cols).toBeGreaterThan(0);
  });

  await softStep('Go to ML > Models > Apply Model...', async () => {
    await page.locator('[name="div-ML"]').click();
    // Submenus: dispatch mouseenter via page.evaluate (Dart ignores Playwright hover)
    await page.evaluate(() => {
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
      if (!models) throw new Error('Models menu item not found');
      const r = models.getBoundingClientRect();
      const ev = (type: string) => new MouseEvent(type, {
        bubbles: true, cancelable: true, view: window,
        clientX: r.left + 5, clientY: r.top + 5,
      });
      models.dispatchEvent(ev('mouseover'));
      models.dispatchEvent(ev('mouseenter'));
      models.dispatchEvent(ev('mousemove'));
    });
    await page.locator('[name="div-ML---Models---Apply-Model..."]').click();
    await page.locator('[name="dialog-Apply-predictive-model"]').waitFor({timeout: 10_000});
    const title = await page.locator('[name="dialog-Apply-predictive-model"] .d4-dialog-title').innerText();
    expect(title.trim()).toBe('Apply predictive model');
  });

  await softStep('Select the TestDemog model', async () => {
    // Wait up to 10s for the Model dropdown to populate
    await page.waitForFunction(() => {
      const sel = document.querySelector('[name="input-host-Model"] select') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 10_000});
    // Select "TestDemog" from the Model select
    const option = await page.evaluate(() => {
      const sel = document.querySelector('[name="input-host-Model"] select') as HTMLSelectElement | null;
      if (!sel) return null;
      const opts = Array.from(sel.options).map((o) => o.textContent || '');
      const idx = opts.findIndex((t) => t.trim() === 'TestDemog');
      if (idx < 0) return {found: false, opts};
      sel.selectedIndex = idx;
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      return {found: true, opts};
    });
    expect(option?.found, `TestDemog not found in Model options: ${JSON.stringify(option?.opts)}`).toBe(true);
  });

  await softStep('Verify prediction column added to the table', async () => {
    // Confirm the dialog
    await page.locator('[name="button-OK"]').click();
    // Allow server round-trip
    await page.waitForTimeout(5_000);
    const result = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      const initial: number = (window as any).__initialColCount ?? 0;
      const initialNames: string[] = (window as any).__initialColNames ?? [];
      const cols = df ? Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i).name) : [];
      const newCols = cols.filter((n: string) => !initialNames.includes(n));
      return {initial, now: cols.length, newCols};
    });
    expect(result.now).toBeGreaterThan(result.initial);
    expect(result.newCols.length).toBeGreaterThan(0);
  });

  if (stepErrors.length > 0)
    throw new Error(`Soft step failures:\n${stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n')}`);
});
