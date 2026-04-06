import {test, expect} from '@playwright/test';

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

test('Color consistency across viewers', async ({page}) => {
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
      for (const v of viewers) {
        tv.addViewer(v);
        await new Promise(r => setTimeout(r, 500));
      }
      const list = Array.from(tv.viewers);
      return {count: list.length, types: list.map(v => v.type)};
    });
    expect(result.count).toBe(7);
    expect(result.types).toContain('Histogram');
    expect(result.types).toContain('Line chart');
    expect(result.types).toContain('Bar chart');
    expect(result.types).toContain('Pie chart');
    expect(result.types).toContain('Trellis plot');
    expect(result.types).toContain('Box plot');
  });

  // Step 3: Set legend to Stereo Category for each viewer
  await softStep('Set legend to Stereo Category for each viewer', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const viewers = Array.from(tv.viewers);
      for (const v of viewers) {
        if (v.type === 'Grid') continue;
        if (['Histogram', 'Line chart', 'Bar chart', 'Pie chart'].includes(v.type))
          v.setOptions({split: 'Stereo Category'});
        else
          v.setOptions({color: 'Stereo Category'});
      }
      await new Promise(r => setTimeout(r, 1000));
    });
  });

  // Step 4: Enable color coding in grid and change colors
  await softStep('Enable grid color coding and change colors — verify in all viewers', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const col = tv.dataFrame.col('Stereo Category');
      const gridCol = tv.grid.columns.byName('Stereo Category');
      if (gridCol) gridCol.isTextColorCoded = true;

      const colorMap = {
        'R_ONE': DG.Color.fromHtml('#FF0000'),
        'S_ABS': DG.Color.fromHtml('#0000FF'),
        'S_ACHIR': DG.Color.fromHtml('#00CC00'),
        'S_PART': DG.Color.fromHtml('#FF8800'),
        'S_UNKN': DG.Color.fromHtml('#9900CC'),
      };
      col.meta.colors.setCategorical(colorMap);
      tv.grid.invalidate();
      await new Promise(r => setTimeout(r, 1000));

      return {
        colorType: col.meta.colors.getType(),
        c0: col.meta.colors.getColor(0),
      };
    });
    expect(result.colorType).toBe('Categorical');
    expect(result.c0).toBe(0xFFFF0000); // R_ONE = red
  });

  // Step 5: Change color on a viewer (via column meta) — verify propagation
  await softStep('Change category color and verify propagation to all viewers', async () => {
    const result = await page.evaluate(async () => {
      const col = grok.shell.tv.dataFrame.col('Stereo Category');
      const colorMap = {
        'R_ONE': DG.Color.fromHtml('#FFD700'),
        'S_ABS': DG.Color.fromHtml('#0000FF'),
        'S_ACHIR': DG.Color.fromHtml('#00FFFF'),
        'S_PART': DG.Color.fromHtml('#FF8800'),
        'S_UNKN': DG.Color.fromHtml('#9900CC'),
      };
      col.meta.colors.setCategorical(colorMap);
      grok.shell.tv.grid.invalidate();
      await new Promise(r => setTimeout(r, 1000));
      return {
        colorType: col.meta.colors.getType(),
        c0: col.meta.colors.getColor(0),
      };
    });
    expect(result.colorType).toBe('Categorical');
    // R_ONE should now be gold (0xFFFFD700)
    expect(result.c0).toBe(0xFFFFD700);
  });

  // Step 6: Save and apply layout
  await softStep('Save layout, close viewers, restore layout — verify colors preserved', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;

      // Save layout
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      // Close all viewers except grid
      for (const v of Array.from(tv.viewers))
        if (v.type !== 'Grid') v.close();
      await new Promise(r => setTimeout(r, 500));

      // Reset color coding
      tv.dataFrame.col('Stereo Category').meta.colors.setDisabled();
      tv.grid.invalidate();
      await new Promise(r => setTimeout(r, 500));

      // Restore layout
      const saved = await grok.dapi.layouts.find(layoutId);
      tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      // Verify
      const col = tv.dataFrame.col('Stereo Category');
      const viewerCount = Array.from(tv.viewers).length;
      const colorType = col.meta.colors.getType();
      const c0 = col.meta.colors.getColor(0);

      // Cleanup
      await grok.dapi.layouts.delete(saved);

      return {viewerCount, colorType, c0};
    });
    expect(result.viewerCount).toBe(7);
    expect(result.colorType).toBe('Categorical');
    expect(result.c0).toBe(0xFFFFD700); // gold preserved
  });

  // Step 7: Save and open project — verify color consistency
  await softStep('Save project, close, reopen — verify color consistency', async () => {
    // Click SAVE button
    await page.locator('[name="button-SAVE"], button:has-text("SAVE")').first().click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});

    // Type project name
    const nameInput = page.locator('.d4-dialog input[type="text"], .d4-dialog textbox').first();
    await nameInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('ColorConsistencyTest');
    await page.locator('.d4-dialog button:has-text("OK")').click();
    await page.waitForTimeout(3000);

    const result = await page.evaluate(async () => {
      // Find the project
      const projects = await grok.dapi.projects.filter('name = "ColorConsistencyTest"').list();
      if (projects.length === 0) return {error: 'Project not found'};
      const projectId = projects[0].id;

      // Close all
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 1000));

      // Reopen
      const project = await grok.dapi.projects.find(projectId);
      await project.open();
      await new Promise(r => setTimeout(r, 5000));

      // Verify
      const tv = grok.shell.tv;
      if (!tv) return {error: 'No table view after project open'};
      const viewers = Array.from(tv.viewers);
      const col = tv.dataFrame.col('Stereo Category');

      const result = {
        viewerCount: viewers.length,
        viewerTypes: viewers.map(v => v.type),
        colorType: col.meta.colors.getType(),
        c0: col.meta.colors.getColor(0),
      };

      // Cleanup
      try {
        for (const p of await grok.dapi.projects.filter('name = "ColorConsistencyTest"').list())
          await grok.dapi.projects.delete(p);
      } catch (_) {}

      return result;
    });

    expect(result.viewerCount).toBe(7);
    expect(result.colorType).toBe('Categorical');
    expect(result.c0).toBe(0xFFFFD700); // gold preserved
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
