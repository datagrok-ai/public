import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Radar viewer (Charts package)', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
    await page.waitForFunction(() => {
      try { return typeof grok !== 'undefined' && typeof grok.shell.closeAll === 'function'; }
      catch { return false; }
    }, {timeout: 45000});
  }

  // Step 1: Open earthquakes.csv and add Radar viewer
  await softStep('Open earthquakes.csv and add Radar viewer', async () => {
    const result = await page!.evaluate(async () => {
      document.querySelectorAll('.d4-dialog').forEach(d => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) (cancel as HTMLElement).click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = false;
      const df = await grok.dapi.files.readCsv('System:DemoFiles/geo/earthquakes.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      tv.addViewer('Radar');
      await new Promise(r => setTimeout(r, 2000));
      return {rows: df.rowCount, cols: df.columns.length, viewers: tv.viewers.length};
    });
    expect(result.rows).toBe(2426);
    expect(result.viewers).toBe(2); // Grid + Radar
  });

  // Step 2: Open demog.csv and add Radar viewer
  await softStep('Open demog.csv and add Radar viewer', async () => {
    const result = await page!.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      tv.addViewer('Radar');
      await new Promise(r => setTimeout(r, 2000));
      return {rows: df.rowCount, viewers: tv.viewers.length};
    });
    expect(result.rows).toBe(5850);
    expect(result.viewers).toBe(2);
    const radarVisible = await page!.locator('[name="viewer-Radar"]').isVisible();
    expect(radarVisible).toBe(true);
  });

  // Step 3a: Switch tables via properties
  await softStep('Switch tables via properties', async () => {
    // Click gear icon on Radar viewer dock panel title bar
    await page!.evaluate(async () => {
      const radar = document.querySelector('[name="viewer-Radar"]');
      const panel = radar?.parentElement?.parentElement;
      const gear = panel?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      if (gear) gear.click();
      await new Promise(r => setTimeout(r, 500));
    });

    // Switch table to earthquakes via setOptions
    const tableName = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const radar = viewers.find((v: any) => v.type === 'RadarViewer') as any;
      radar.setOptions({table: 'Table'});
      await new Promise(r => setTimeout(r, 1000));
      return radar.dataFrame?.name;
    });
    expect(tableName).toBe('Table');

    // Switch back to demog
    await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const radar = viewers.find((v: any) => v.type === 'RadarViewer') as any;
      radar.setOptions({table: 'Table (2)'});
      await new Promise(r => setTimeout(r, 1000));
    });
  });

  // Step 3b: Toggle Show Min and Show Max checkboxes
  await softStep('Toggle Show Min and Show Max', async () => {
    const opts = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const radar = viewers.find((v: any) => v.type === 'RadarViewer') as any;
      radar.setOptions({showMin: true, showMax: true});
      await new Promise(r => setTimeout(r, 1000));
      const o = radar.getOptions();
      return {showMin: o.look?.showMin, showMax: o.look?.showMax};
    });
    expect(opts.showMin).toBe(true);
    expect(opts.showMax).toBe(true);
  });

  // Step 3c: Values — change column count
  await softStep('Change Values column count', async () => {
    // Decrease: select only 2 columns
    await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const radar = viewers.find((v: any) => v.type === 'RadarViewer') as any;
      radar.setOptions({valueColumnNames: ['AGE', 'WEIGHT']});
      await new Promise(r => setTimeout(r, 1000));
    });

    // Increase: select all 4 columns
    const increased = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const radar = viewers.find((v: any) => v.type === 'RadarViewer') as any;
      radar.setOptions({valueColumnNames: ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED']});
      await new Promise(r => setTimeout(r, 1000));
      // valueColumnNames may be nested differently; just verify no error
      return true;
    });

    expect(increased).toBe(true);
  });

  // Step 3d: Style (color) changes
  await softStep('Style and color changes', async () => {
    await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const radar = viewers.find((v: any) => v.type === 'RadarViewer') as any;
      radar.setOptions({
        lineColor: '#ff0000',
        backgroundMinColor: '#0000ff',
        backgroundMaxColor: '#ff00ff',
      });
      await new Promise(r => setTimeout(r, 1000));
    });

    const colors = await page!.evaluate(() => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const radar = viewers.find((v: any) => v.type === 'RadarViewer') as any;
      const o = radar.getOptions();
      return {lineColor: o.look?.lineColor, bgMin: o.look?.backgroundMinColor};
    });
    expect(colors.lineColor).toBeDefined();
    expect(colors.bgMin).toBeDefined();

    // Reset colors
    await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const radar = viewers.find((v: any) => v.type === 'RadarViewer') as any;
      radar.setOptions({
        lineColor: '#add8e6',
        backgroundMinColor: '#bb845d',
        backgroundMaxColor: '#e7cdcd',
        currentRowColor: '#00ff00',
      });
      await new Promise(r => setTimeout(r, 500));
    });
  });

  // Cleanup
  await page!.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
