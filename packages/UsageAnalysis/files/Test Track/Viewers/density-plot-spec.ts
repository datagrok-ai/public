import {test, expect, chromium} from '@playwright/test';

const baseUrl = 'http://localhost:8888';
const datasetPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Density plot', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  const page = context.pages()[0] || await context.newPage();

  // Phase 1: Navigate
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell, {timeout: 15000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(() => resolve(undefined), 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Step 1: Verify dataset loaded
  await softStep('Open demog dataset', async () => {
    const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(rowCount).toBe(5850);
  });

  // Step 2: Add Density plot from Toolbox
  await softStep('Add Density plot from Viewers tab', async () => {
    await page.locator('[name="icon-density-plot"]').click();
    await page.locator('[name="viewer-Density-plot"]').waitFor({timeout: 10000});
    const hasViewer = await page.evaluate(() => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++)
        if (grok.shell.tv.viewers[i].type === 'Density plot') return true;
      return false;
    });
    expect(hasViewer).toBe(true);
  });

  // Step 3a: Change X and Y axes
  await softStep('Change X and Y axes on the viewer', async () => {
    await page.evaluate(async () => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++) {
        if (grok.shell.tv.viewers[i].type === 'Density plot') {
          grok.shell.tv.viewers[i].props.xColumnName = 'WEIGHT';
          grok.shell.tv.viewers[i].props.yColumnName = 'HEIGHT';
          break;
        }
      }
    });
    await page.waitForTimeout(1000);

    const axes = await page.evaluate(() => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++) {
        if (grok.shell.tv.viewers[i].type === 'Density plot')
          return {x: grok.shell.tv.viewers[i].props.xColumnName, y: grok.shell.tv.viewers[i].props.yColumnName};
      }
    });
    expect(axes?.x).toBe('WEIGHT');
    expect(axes?.y).toBe('HEIGHT');
  });

  // Step 3b: Zoom in and out with mouse wheel
  await softStep('Zoom in and out with mouse wheel', async () => {
    const canvas = page.locator('[name="viewer-Density-plot"] canvas');
    const box = await canvas.boundingBox();
    expect(box).toBeTruthy();

    const cx = box!.x + box!.width * 0.5;
    const cy = box!.y + box!.height * 0.5;

    await page.mouse.move(cx, cy);
    await page.mouse.wheel(0, -120);
    await page.waitForTimeout(300);
    await page.mouse.wheel(0, -120);
    await page.waitForTimeout(300);
    await page.mouse.wheel(0, 120);
    await page.waitForTimeout(300);
    await page.mouse.wheel(0, 120);
    await page.waitForTimeout(300);

    await expect(canvas).toBeVisible();
  });

  // Step 4-5: Modify properties via JS API
  await softStep('Modify properties — change axes again', async () => {
    await page.evaluate(async () => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++) {
        if (grok.shell.tv.viewers[i].type === 'Density plot') {
          grok.shell.tv.viewers[i].props.xColumnName = 'HEIGHT';
          grok.shell.tv.viewers[i].props.yColumnName = 'WEIGHT';
          break;
        }
      }
    });
    await page.waitForTimeout(1000);

    const axes = await page.evaluate(() => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++) {
        if (grok.shell.tv.viewers[i].type === 'Density plot')
          return {x: grok.shell.tv.viewers[i].props.xColumnName, y: grok.shell.tv.viewers[i].props.yColumnName};
      }
    });
    expect(axes?.x).toBe('HEIGHT');
    expect(axes?.y).toBe('WEIGHT');
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
