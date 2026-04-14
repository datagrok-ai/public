import {test, expect} from '@playwright/test';

const baseUrl = 'http://localhost:8888/';
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

test('Legend filtering', async ({page}) => {
  // Phase 1: Navigate and login
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined', {timeout: 120000});
  const loginField = page.getByRole('textbox', {name: 'Login or Email'});
  if (await loginField.count() > 0) {
    await loginField.click({force: true});
    await page.keyboard.type('admin');
    await page.locator('input[type="password"]').last().click({force: true});
    await page.keyboard.type('admin');
    await page.keyboard.press('Enter');
  }
  await page.waitForFunction(() => {
    try { grok.shell.settings.showFiltersIconsConstantly; return true; } catch { return false; }
  }, {timeout: 120000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Step 2: Add viewers via toolbox icons
  await softStep('Add viewers', async () => {
    const icons = ['icon-scatter-plot', 'icon-histogram', 'icon-line-chart',
      'icon-bar-chart', 'icon-pie-chart', 'icon-trellis-plot', 'icon-box-plot'];
    for (const name of icons) {
      await page.locator(`[name="${name}"]`).click();
      await page.waitForTimeout(1000);
    }
    const count = await page.evaluate(() => grok.shell.tv.viewers.length);
    expect(count).toBe(8);
  });

  // Step 3: Set legend to Stereo Category
  await softStep('Set legend to Stereo Category', async () => {
    await page.evaluate(async () => {
      for (const v of Array.from(grok.shell.tv.viewers)) {
        if (v.type === 'Scatter plot') v.setOptions({color: 'Stereo Category', legendVisibility: 'Always'});
        else if (v.type === 'Histogram') v.setOptions({split: 'Stereo Category', legendVisibility: 'Always'});
        else if (v.type === 'Line chart') v.setOptions({split: 'Stereo Category', legendVisibility: 'Always'});
        else if (v.type === 'Bar chart') v.setOptions({stack: 'Stereo Category', legendVisibility: 'Always'});
        else if (v.type === 'Box plot') v.setOptions({color: 'Stereo Category', legendVisibility: 'Always'});
        else if (v.type === 'Pie chart') v.setOptions({legendVisibility: 'Always'});
        else if (v.type === 'Trellis plot') v.setOptions({legendVisibility: 'Always'});
      }
      await new Promise(r => setTimeout(r, 2000));
    });
  });

  // Step 4: Open Filter Panel, apply filter, check legend
  await softStep('Filter to R_ONE+S_ACHIR, check legend shows only 2 categories', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      df.filter.init((i) => {
        const val = df.col('Stereo Category').get(i);
        return val === 'R_ONE' || val === 'S_ACHIR';
      });
      df.filter.fireChanged();
      await new Promise(r => setTimeout(r, 1000));
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot');
      return {
        filteredRows: df.filter.trueCount,
        legend: Array.from(sp.root.querySelectorAll('.d4-legend-value')).map(el => el.textContent),
      };
    });
    expect(result.filteredRows).toBe(2247);
    expect(result.legend).toEqual(['R_ONE', 'S_ACHIR']);
  });

  // Step 5: Save and apply layout
  await softStep('Save and apply layout', async () => {
    const count = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      const saved = await grok.dapi.layouts.find(layout.id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
      await grok.dapi.layouts.delete(saved);
      return grok.shell.tv.viewers.length;
    });
    expect(count).toBeGreaterThanOrEqual(8);
  });

  // Step 6: Reset filters
  await softStep('Reset filters', async () => {
    const rows = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      await new Promise(r => setTimeout(r, 500));
      return df.filter.trueCount;
    });
    expect(rows).toBe(3624);
  });

  // Step 7: Set in-viewer Filter
  await softStep('Set in-viewer filter to R_ONE, S_UNKN', async () => {
    const result = await page.evaluate(async () => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot');
      sp.setOptions({filter: '${Stereo Category} in ["R_ONE", "S_UNKN"]'});
      await new Promise(r => setTimeout(r, 1000));
      return {
        dfRows: grok.shell.tv.dataFrame.filter.trueCount,
        legend: Array.from(sp.root.querySelectorAll('.d4-legend-value')).map(el => el.textContent),
      };
    });
    expect(result.dfRows).toBe(3624);
    expect(result.legend).toEqual(['R_ONE', 'S_UNKN']);
  });

  // Step 8: Apply additional dataframe filter
  await softStep('Additional dataframe filter, check intersection in legend', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.init((i) => {
        const val = df.col('Stereo Category').get(i);
        return val === 'R_ONE' || val === 'S_UNKN' || val === 'S_ABS';
      });
      df.filter.fireChanged();
      await new Promise(r => setTimeout(r, 1000));
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot');
      const hist = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Histogram');
      return {
        spLegend: Array.from(sp.root.querySelectorAll('.d4-legend-value')).map(el => el.textContent),
        histLegend: Array.from(hist.root.querySelectorAll('.d4-legend-value')).map(el => el.textContent),
      };
    });
    // SP has in-viewer filter (R_ONE, S_UNKN) AND df filter (R_ONE, S_UNKN, S_ABS) = intersection
    expect(result.spLegend).toEqual(['R_ONE', 'S_UNKN']);
    // Histogram has only df filter
    expect(result.histLegend).toEqual(['R_ONE', 'S_ABS', 'S_UNKN']);
  });

  // Step 9: Filter via legend click
  await softStep('Legend click filters in-viewer', async () => {
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot');
      sp.setOptions({filter: ''});
      await new Promise(r => setTimeout(r, 500));
    });
    const result = await page.evaluate(async () => {
      const sp = document.querySelector('[name="viewer-Scatter-plot"]');
      const items = sp.querySelectorAll('.d4-legend-item');
      items[0].dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true}));
      await new Promise(r => setTimeout(r, 500));
      return Array.from(sp.querySelectorAll('.d4-legend-item')).map(i => ({
        text: i.querySelector('.d4-legend-value')?.textContent,
        current: i.classList.contains('d4-legend-item-current'),
      }));
    });
    expect(result[0].current).toBe(true);
    expect(result[1].current).toBe(false);
    // Reset
    await page.evaluate(async () => {
      const sp = document.querySelector('[name="viewer-Scatter-plot"]');
      for (const item of sp.querySelectorAll('.d4-legend-item'))
        if (!item.classList.contains('d4-legend-item-current'))
          item.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, ctrlKey: true}));
      await new Promise(r => setTimeout(r, 500));
    });
  });

  // Step 11: Different Row Source values
  await softStep('Row Source affects legend (All, Filtered, Selected)', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.init((i) => {
        const val = df.col('Stereo Category').get(i);
        return val === 'R_ONE' || val === 'S_ACHIR';
      });
      df.filter.fireChanged();
      df.selection.init((i) => df.col('Stereo Category').get(i) === 'R_ONE');
      await new Promise(r => setTimeout(r, 500));
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot');
      sp.setOptions({rowSource: 'All'});
      await new Promise(r => setTimeout(r, 500));
      const allLegend = Array.from(sp.root.querySelectorAll('.d4-legend-value')).map(el => el.textContent);
      sp.setOptions({rowSource: 'Filtered'});
      await new Promise(r => setTimeout(r, 500));
      const filteredLegend = Array.from(sp.root.querySelectorAll('.d4-legend-value')).map(el => el.textContent);
      sp.setOptions({rowSource: 'Selected'});
      await new Promise(r => setTimeout(r, 500));
      const selectedLegend = Array.from(sp.root.querySelectorAll('.d4-legend-value')).map(el => el.textContent);
      sp.setOptions({rowSource: 'Filtered'});
      df.filter.setAll(true);
      df.selection.setAll(false);
      return {allLegend, filteredLegend, selectedLegend};
    });
    expect(result.allLegend.length).toBe(5);
    expect(result.filteredLegend).toEqual(['R_ONE', 'S_ACHIR']);
    expect(result.selectedLegend).toEqual(['R_ONE']);
  });

  // Section 2: Bar chart edge case
  await softStep('Bar chart edge case: stacked legend shows only filtered categories', async () => {
    const result = await page.evaluate(async () => {
      // Reset filters and selection first
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      df.selection.setAll(false);
      await new Promise(r => setTimeout(r, 500));

      const bc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Bar chart');
      bc.setOptions({
        valueColumnName: 'CAST Idea ID',
        splitColumnName: 'Stereo Category',
        stackColumnName: 'Primary scaffold name',
        includeNulls: false,
        legendVisibility: 'Always',
      });
      await new Promise(r => setTimeout(r, 2000));

      // Filter to 3 scaffold categories
      const keep = ['AMINOPYRROLIDINE', 'TRISUBSTITUTED', 'UNSUBSTDIAZA'];
      df.filter.init((i) => keep.includes(df.col('Primary scaffold name').get(i)));
      df.filter.fireChanged();
      await new Promise(r => setTimeout(r, 1500));

      // The bar chart legend shows the STACK column categories (Primary scaffold name)
      const legend = Array.from(bc.root.querySelectorAll('.d4-legend-value')).map(el => el.textContent);
      return {legend, filteredRows: df.filter.trueCount};
    });
    expect(result.legend.length).toBe(3);
    expect(result.legend).toContain('AMINOPYRROLIDINE');
    expect(result.legend).toContain('TRISUBSTITUTED');
    expect(result.legend).toContain('UNSUBSTDIAZA');
  });

  // Cleanup
  await softStep('Close All', async () => {
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1000);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
