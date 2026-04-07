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

test('Legend visibility and positioning', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  const page = context.pages()[0] || await context.newPage();

  // Phase 1: Navigate
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell, {timeout: 15000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch(e) {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(() => resolve(undefined), 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Add viewers and set legends
  await softStep('Add viewers and set color legends', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const viewerTypes = ['Scatter plot', 'Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'];
      for (const type of viewerTypes) {
        tv.addViewer(type);
        await new Promise(r => setTimeout(r, 500));
      }

      // Remove duplicate viewers if any
      const viewers = Array.from(tv.viewers);
      const seen = new Set();
      seen.add('Grid');
      for (const v of viewers) {
        if (v.type === 'Grid') continue;
        if (seen.has(v.type))
          v.close();
        else
          seen.add(v.type);
      }
      await new Promise(r => setTimeout(r, 500));

      // Set color/split/stack columns
      const remaining = Array.from(tv.viewers);
      for (const v of remaining) {
        if (v.type === 'Scatter plot') v.setOptions({colorColumnName: 'Stereo Category'});
        else if (v.type === 'Histogram') v.setOptions({splitColumnName: 'Stereo Category'});
        else if (v.type === 'Line chart') v.setOptions({splitColumnName: 'Stereo Category'});
        else if (v.type === 'Bar chart') v.setOptions({stackColumnName: 'Stereo Category'});
        else if (v.type === 'Pie chart') v.setOptions({categoryColumnName: 'Stereo Category'});
        else if (v.type === 'Trellis plot') v.setOptions({colorColumnName: 'Stereo Category'});
        else if (v.type === 'Box plot') v.setOptions({colorColumnName: 'Stereo Category'});
      }
      await new Promise(r => setTimeout(r, 1000));

      const count = Array.from(tv.viewers).filter(v => v.type !== 'Grid').length;
      return {count};
    });

    expect(result.count).toBe(7);
  });

  await softStep('Verify legends are visible', async () => {
    const result = await page.evaluate(() => {
      const tv = grok.shell.tv;
      const viewers = Array.from(tv.viewers).filter(v => v.type !== 'Grid');
      const withLegend = viewers.filter(v => {
        const legend = v.root.querySelector('.d4-legend');
        return legend && legend.offsetHeight > 0;
      });
      return {total: viewers.length, withLegend: withLegend.length};
    });

    // At least Scatter Plot, Histogram, Pie Chart, Trellis Plot should have legends
    expect(result.withLegend).toBeGreaterThanOrEqual(4);
  });

  await softStep('Verify changing color column updates legend', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const sp = Array.from(tv.viewers).find(v => v.type === 'Scatter plot');
      sp.setOptions({colorColumnName: 'Scaffold Names'});
      await new Promise(r => setTimeout(r, 1000));
      const legendTexts = Array.from(sp.root.querySelectorAll('.d4-legend-value')).map(el => el.textContent.trim());
      const hasScaffold = legendTexts.some(t => t.includes('AMINO') || t.includes('TRISUBSTITUTED'));

      sp.setOptions({colorColumnName: 'Stereo Category'});
      await new Promise(r => setTimeout(r, 500));
      const legendTextsAfter = Array.from(sp.root.querySelectorAll('.d4-legend-value')).map(el => el.textContent.trim());
      const hasStereo = legendTextsAfter.includes('R_ONE');

      return {hasScaffold, hasStereo};
    });

    expect(result.hasScaffold).toBe(true);
    expect(result.hasStereo).toBe(true);
  });

  await softStep('CTRL+click multi-select on legend', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const sp = Array.from(tv.viewers).find(v => v.type === 'Scatter plot');
      const legendItems = sp.root.querySelectorAll('.d4-legend-item');

      legendItems[0].click();
      await new Promise(r => setTimeout(r, 500));
      const firstCurrent = legendItems[0].classList.contains('d4-legend-item-current');

      legendItems[1].dispatchEvent(new MouseEvent('click', {bubbles: true, ctrlKey: true}));
      await new Promise(r => setTimeout(r, 500));
      const bothCurrent = legendItems[0].classList.contains('d4-legend-item-current') &&
                          legendItems[1].classList.contains('d4-legend-item-current');

      return {firstCurrent, bothCurrent};
    });

    expect(result.firstCurrent).toBe(true);
    expect(result.bothCurrent).toBe(true);
  });

  await softStep('Save and apply layout', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(() => resolve(undefined), 3000);
      });

      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const viewers = Array.from(tv2.viewers).filter(v => v.type !== 'Grid');
      const withLegend = viewers.filter(v => {
        const legend = v.root.querySelector('.d4-legend');
        return legend && legend.offsetHeight > 0;
      });

      await grok.dapi.layouts.delete(saved);
      return {viewerCount: viewers.length, legendCount: withLegend.length};
    });

    expect(result.viewerCount).toBe(7);
    expect(result.legendCount).toBeGreaterThanOrEqual(4);
  });

  await softStep('Set Visibility=Always, Position=Auto and verify', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const viewers = Array.from(tv.viewers).filter(v => v.type !== 'Grid');
      for (const v of viewers) {
        try {
          v.setOptions({legendVisibility: 'Always', legendPosition: 'Auto'});
        } catch (e) {}
      }
      await new Promise(r => setTimeout(r, 1000));

      const withLegend = viewers.filter(v => {
        const legend = v.root.querySelector('.d4-legend');
        return legend && legend.offsetHeight > 0;
      });
      return {total: viewers.length, withLegend: withLegend.length};
    });

    // With Visibility=Always, most viewers should show legends
    expect(result.withLegend).toBeGreaterThanOrEqual(5);
  });

  await softStep('Set Visibility=Auto and verify auto-hide on small size', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const sp = Array.from(tv.viewers).find(v => v.type === 'Scatter plot');
      sp.setOptions({legendVisibility: 'Auto'});
      await new Promise(r => setTimeout(r, 500));

      // Resize to small
      const container = document.querySelector('[name="viewer-Scatter-plot"]');
      const parent = container?.closest('.panel-base');
      if (!parent) return {error: 'no parent'};

      const origStyle = parent.style.cssText;
      parent.style.width = '100px';
      parent.style.height = '100px';
      window.dispatchEvent(new Event('resize'));
      await new Promise(r => setTimeout(r, 1000));

      const legend = sp.root.querySelector('.d4-legend');
      const hiddenWhenSmall = !legend || legend.offsetHeight === 0;

      // Restore
      parent.style.cssText = origStyle;
      window.dispatchEvent(new Event('resize'));
      await new Promise(r => setTimeout(r, 1000));

      return {hiddenWhenSmall};
    });

    expect(result.hiddenWhenSmall).toBe(true);
  });

  await softStep('Set corner positions and save layout', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const viewers = Array.from(tv.viewers).filter(v => v.type !== 'Grid');
      const positions = ['Top Right', 'Top Left', 'Bottom Right', 'Bottom Left', 'Top Right', 'Bottom Left', 'Top Right'];
      for (let i = 0; i < viewers.length; i++) {
        try {
          viewers[i].setOptions({legendPosition: positions[i], legendVisibility: 'Always'});
        } catch (e) {}
      }
      await new Promise(r => setTimeout(r, 1500));

      // Save layout
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      // Close and reopen
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(() => resolve(undefined), 3000);
      });

      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const restoredViewers = Array.from(tv2.viewers).filter(v => v.type !== 'Grid');
      const withLegend = restoredViewers.filter(v => {
        const legend = v.root.querySelector('.d4-legend');
        return legend && legend.offsetHeight > 0;
      });

      await grok.dapi.layouts.delete(saved);
      return {viewerCount: restoredViewers.length, legendCount: withLegend.length};
    });

    expect(result.viewerCount).toBe(7);
    expect(result.legendCount).toBeGreaterThanOrEqual(4);
  });

  await softStep('Close all', async () => {
    await page.evaluate(() => grok.shell.closeAll());
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
