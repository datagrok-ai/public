import { test, expect, Page } from '@playwright/test';

const baseUrl = 'https://dev.datagrok.ai';
const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Radar viewer (Charts package)', async ({ page }) => {
  test.setTimeout(120_000);
  // Phase 1: Navigate
  await page.goto(baseUrl);
  await page.waitForFunction(() => {
    return typeof grok !== 'undefined' && grok.shell && grok.shell.views;
  }, {timeout: 30000});

  // Phase 2: Open earthquakes.csv
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/geo/earthquakes.csv');
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Step 1: Add Radar viewer for earthquakes.csv
  await softStep('Step 1: Add Radar viewer for earthquakes.csv', async () => {
    await page.evaluate(async () => {
      await grok.shell.tv.addViewer('Radar');
    });
    await page.locator('[name="viewer-Radar"]').waitFor({timeout: 10000});
    const radarVisible = await page.locator('[name="viewer-Radar"]').isVisible();
    expect(radarVisible).toBe(true);
  });

  // Step 2: Open demog.csv and add Radar viewer
  await softStep('Step 2: Add Radar viewer for demog.csv', async () => {
    const result = await page.evaluate(async () => {
      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      await tv2.addViewer('Radar');
      await new Promise(r => setTimeout(r, 2000));
      return { rows: df2.rowCount };
    });
    expect(result.rows).toBe(5850);
    const radarVisible = await page.locator('[name="viewer-Radar"]').isVisible();
    expect(radarVisible).toBe(true);
  });

  // Step 3a: Open properties panel and switch tables
  await softStep('Step 3a: Switch tables via properties', async () => {
    // Click gear icon on the Radar viewer dock panel title bar
    const switched = await page.evaluate(async () => {
      const radar = document.querySelector('[name="viewer-Radar"]');
      if (!radar) return { error: 'no radar' };
      const dockPanel = radar.closest('.panel-content');
      const titleBar = dockPanel?.parentElement?.querySelector('.panel-titlebar');
      const settingsIcon = titleBar?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      if (settingsIcon) settingsIcon.click();
      await new Promise(r => setTimeout(r, 500));

      // Find the Table combobox and switch to "Table" (earthquakes)
      const combo = document.querySelector('.grok-prop-panel select') as HTMLSelectElement;
      if (!combo) return { error: 'no combo' };
      combo.value = 'Table';
      combo.dispatchEvent(new Event('change', { bubbles: true }));
      await new Promise(r => setTimeout(r, 1000));

      const viewers = Array.from(grok.shell.tv.viewers);
      const radarViewer = viewers.find((v: any) => v.type === 'RadarViewer') as any;
      return { tableName: radarViewer?.dataFrame?.name };
    });
    expect(switched).not.toHaveProperty('error');
  });

  // Step 3b: Test column selection dialog
  await softStep('Step 3b: Column selection (Values)', async () => {
    // Click the "..." button to open Select Columns dialog
    const dialogOpened = await page.evaluate(async () => {
      const propPanel = document.querySelector('.grok-prop-panel');
      const btn = propPanel?.querySelector('button');
      if (btn) btn.click();
      await new Promise(r => setTimeout(r, 500));
      return !!document.querySelector('.d4-dialog');
    });
    expect(dialogOpened).toBe(true);

    // Click "All" then OK
    await page.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog');
      const allLink = Array.from(dialog?.querySelectorAll('*') ?? [])
        .find(el => el.textContent?.trim() === 'All' && el.tagName !== 'OPTION') as HTMLElement;
      if (allLink) allLink.click();
      await new Promise(r => setTimeout(r, 300));
      const okBtn = dialog?.querySelector('[name="button-OK"]') as HTMLElement;
      if (okBtn) okBtn.click();
    });
    await page.waitForTimeout(500);
  });

  // Step 3c: Style/color changes
  await softStep('Step 3c: Style and color changes', async () => {
    const result = await page.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const radar = viewers.find((v: any) => v.type === 'RadarViewer') as any;
      if (!radar) return { error: 'no radar' };

      radar.setOptions({
        currentRowColor: 4294901760,
        showMin: true,
        showMax: true,
        showValues: true,
      });
      await new Promise(r => setTimeout(r, 1000));

      const opts = radar.getOptions();
      return {
        showValues: opts.look?.showValues,
        showMin: opts.look?.showMin,
        showMax: opts.look?.showMax,
      };
    });
    expect(result).not.toHaveProperty('error');
    expect(result.showValues).toBe(true);
    expect(result.showMin).toBe(true);
    expect(result.showMax).toBe(true);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
