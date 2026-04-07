import {test, expect, chromium} from '@playwright/test';

const baseUrl = 'http://localhost:8888';
const datasetPath = 'System:DemoFiles/SPGI.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Color consistency across viewers', async () => {
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

  // Phase 3: Test case

  // Step 1-2: Add viewers
  await softStep('Add viewers (histogram, line chart, bar chart, pie chart, trellis plot, box plot)', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const viewers = ['Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'];
      for (const v of viewers)
        tv.addViewer(v);
      await new Promise(r => setTimeout(r, 2000));
      return {count: tv.viewers.length};
    });
    expect(result.count).toBe(7);
  });

  // Step 3: Set legend to Stereo Category for each viewer
  await softStep('Set legend to Stereo Category for each viewer', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      for (let i = 0; i < tv.viewers.length; i++) {
        const v = tv.viewers[i];
        if (v.type === 'Grid') continue;
        if (['Histogram', 'Line chart', 'Bar chart'].includes(v.type))
          v.setOptions({split: 'Stereo Category'});
        else if (v.type === 'Pie chart')
          v.setOptions({category: 'Stereo Category'});
        else
          v.setOptions({color: 'Stereo Category'});
      }
      await new Promise(r => setTimeout(r, 1000));
    });
  });

  // Step 4: Enable color coding in grid and set custom colors
  await softStep('Enable grid color coding and set custom categorical colors', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const col = tv.dataFrame.columns.byName('Stereo Category');
      const gridCol = tv.grid.columns.byName('Stereo Category');
      if (gridCol) gridCol.isTextColorCoded = true;

      col.meta.colors.setCategorical({
        'R_ONE': DG.Color.fromHtml('#0000FF'),
        'S_ABS': DG.Color.fromHtml('#FF00FF'),
        'S_ACHIR': DG.Color.fromHtml('#00AA00'),
        'S_PART': DG.Color.fromHtml('#FF8800'),
        'S_UNKN': DG.Color.fromHtml('#888888'),
      });
      tv.grid.invalidate();
      await new Promise(r => setTimeout(r, 1000));

      const colorTag = col.getTag('.color-coding-categorical');
      return {colorType: col.meta.colors.getType(), hasTag: !!colorTag};
    });
    expect(result.colorType).toBe('Categorical');
    expect(result.hasTag).toBe(true);
  });

  // Step 5: Change a category color and verify propagation
  await softStep('Change R_ONE to red and verify propagation to all viewers', async () => {
    const result = await page.evaluate(async () => {
      const col = grok.shell.tv.dataFrame.columns.byName('Stereo Category');
      col.meta.colors.setCategorical({
        'R_ONE': DG.Color.fromHtml('#FF0000'),
        'S_ABS': DG.Color.fromHtml('#FF00FF'),
        'S_ACHIR': DG.Color.fromHtml('#00AA00'),
        'S_PART': DG.Color.fromHtml('#FF8800'),
        'S_UNKN': DG.Color.fromHtml('#888888'),
      });
      grok.shell.tv.grid.invalidate();
      await new Promise(r => setTimeout(r, 1000));

      const colorTag = col.getTag('.color-coding-categorical');
      const parsed = JSON.parse(colorTag);
      return {r_one_color: parsed['R_ONE']};
    });
    expect(result.r_one_color).toBe(4294901760);
  });

  // Step 6: Save and apply layout
  await softStep('Save layout, close viewers, restore — verify colors preserved', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;

      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1500));

      for (let i = tv.viewers.length - 1; i >= 0; i--) {
        if (tv.viewers[i].type !== 'Grid') tv.viewers[i].close();
      }
      await new Promise(r => setTimeout(r, 500));

      const saved = await grok.dapi.layouts.find(layoutId);
      tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const col = tv.dataFrame.columns.byName('Stereo Category');
      const viewerCount = tv.viewers.length;
      const colorTag = col.getTag('.color-coding-categorical');
      const parsed = JSON.parse(colorTag);

      await grok.dapi.layouts.delete(saved);
      return {viewerCount, r_one_color: parsed['R_ONE']};
    });
    expect(result.viewerCount).toBe(7);
    expect(result.r_one_color).toBe(4294901760);
  });

  // Step 7: Close all, reopen fresh dataset + layout — verify color consistency
  await softStep('Fresh dataset with saved layout — verify color consistency', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1500));

      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 1000));

      const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(() => resolve(undefined), 3000);
      });

      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 4000));

      const col = df.columns.byName('Stereo Category');
      const viewerCount = tv2.viewers.length;
      const colorTag = col.getTag('.color-coding-categorical');
      const parsed = colorTag ? JSON.parse(colorTag) : null;

      await grok.dapi.layouts.delete(saved);
      return {viewerCount, r_one_color: parsed?.['R_ONE'], hasTag: !!colorTag};
    });
    expect(result.viewerCount).toBe(7);
    expect(result.hasTag).toBe(true);
    expect(result.r_one_color).toBe(4294901760);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
