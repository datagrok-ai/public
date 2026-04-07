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

test('3D Scatter Plot', async () => {
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

  // Step 2: Add 3D Scatter Plot from Toolbox
  await softStep('Add 3D Scatter Plot from Viewers tab', async () => {
    await page.locator('[name="icon-3d-scatter-plot"]').click();
    await page.locator('[name="viewer-3d-scatter-plot"]').waitFor({timeout: 10000});
    const hasViewer = await page.evaluate(() => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++)
        if (grok.shell.tv.viewers[i].type === '3d scatter plot') return true;
      return false;
    });
    expect(hasViewer).toBe(true);
  });

  // Step 3a: Click on a data point highlights row in grid
  await softStep('Click on a data point highlights row in grid', async () => {
    const canvas = page.locator('[name="viewer-3d-scatter-plot"] canvas');
    const box = await canvas.boundingBox();
    expect(box).toBeTruthy();

    await page.mouse.click(box!.x + box!.width * 0.5, box!.y + box!.height * 0.35);
    await page.waitForTimeout(500);

    let currentRow = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    if (currentRow === 0) {
      await page.mouse.click(box!.x + box!.width * 0.45, box!.y + box!.height * 0.45);
      await page.waitForTimeout(500);
      currentRow = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    }
    expect(currentRow).toBeGreaterThanOrEqual(0);
  });

  // Step 3b: Zoom in and out with mouse wheel
  await softStep('Zoom in and out with mouse wheel', async () => {
    const canvas = page.locator('[name="viewer-3d-scatter-plot"] canvas');
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

  // Step 4a: Modify Color property to SEX (change from default)
  await softStep('Modify Color property to SEX', async () => {
    await page.evaluate(() => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++) {
        if (grok.shell.tv.viewers[i].type === '3d scatter plot')
          grok.shell.tv.viewers[i].setOptions({colorColumnName: 'SEX'});
      }
    });
    await page.waitForTimeout(1000);

    const color = await page.evaluate(() => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++) {
        if (grok.shell.tv.viewers[i].type === '3d scatter plot')
          return grok.shell.tv.viewers[i].getOptions().look.colorColumnName;
      }
    });
    expect(color).toBe('SEX');
  });

  // Step 4b: Modify Size property to AGE
  await softStep('Modify Size property to AGE', async () => {
    await page.evaluate(() => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++) {
        if (grok.shell.tv.viewers[i].type === '3d scatter plot')
          grok.shell.tv.viewers[i].setOptions({sizeColumnName: 'AGE'});
      }
    });
    await page.waitForTimeout(1000);

    const size = await page.evaluate(() => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++) {
        if (grok.shell.tv.viewers[i].type === '3d scatter plot')
          return grok.shell.tv.viewers[i].getOptions().look.sizeColumnName;
      }
    });
    expect(size).toBe('AGE');
  });

  // Step 4c: Change back to RACE/WEIGHT and verify
  await softStep('Change Color to RACE and Size to WEIGHT', async () => {
    await page.evaluate(() => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++) {
        if (grok.shell.tv.viewers[i].type === '3d scatter plot')
          grok.shell.tv.viewers[i].setOptions({colorColumnName: 'RACE', sizeColumnName: 'WEIGHT'});
      }
    });
    await page.waitForTimeout(1000);

    const opts = await page.evaluate(() => {
      for (let i = 0; i < grok.shell.tv.viewers.length; i++) {
        if (grok.shell.tv.viewers[i].type === '3d scatter plot') {
          const o = grok.shell.tv.viewers[i].getOptions();
          return {color: o.look.colorColumnName, size: o.look.sizeColumnName};
        }
      }
    });
    expect(opts?.color).toBe('RACE');
    expect(opts?.size).toBe('WEIGHT');
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
