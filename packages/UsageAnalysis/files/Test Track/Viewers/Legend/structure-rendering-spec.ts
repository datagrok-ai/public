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

async function openSPGI(page: any) {
  await page.evaluate(async (path: string) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(() => resolve(undefined), 3000);
    });
    // Ensure Molecule semType on Core and Structure columns
    try { df.columns.byName('Core').semType = 'Molecule'; } catch(_) {}
    try { df.columns.byName('Structure').semType = 'Molecule'; } catch(_) {}
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
}

test('Structure rendering in legends', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  const page = context.pages()[0] || await context.newPage();

  // Phase 1: Navigate
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell, {timeout: 15000});

  // Step 1-2: Open SPGI and add viewers
  await softStep('Open SPGI and add viewers', async () => {
    await openSPGI(page);

    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const viewerTypes = ['Scatter plot', 'Histogram', 'Line chart', 'Bar chart', 'Pie chart', 'Trellis plot', 'Box plot'];
      for (const v of viewerTypes) {
        tv.addViewer(v);
        await new Promise(r => setTimeout(r, 500));
      }
      await new Promise(r => setTimeout(r, 1000));
      return {count: tv.viewers.length};
    });
    expect(result.count).toBe(8);
  });

  // Step 3-4: Add structure legend (Core) to each viewer
  await softStep('Add structure legend (Core) to viewers', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      for (let i = 0; i < tv.viewers.length; i++) {
        const v = tv.viewers[i];
        try {
          if (v.type === 'Scatter plot') {
            v.props.markersColumnName = 'Core';
            v.props.colorColumnName = 'Core';
          }
          else if (['Histogram', 'Bar chart', 'Line chart'].includes(v.type))
            v.setOptions({split: 'Core'});
          else if (v.type === 'Pie chart')
            v.setOptions({category: 'Core'});
          else if (['Trellis plot', 'Box plot'].includes(v.type))
            v.setOptions({color: 'Core'});
        } catch(_) {}
      }
      await new Promise(r => setTimeout(r, 2000));

      let sp = null;
      for (let i = 0; i < tv.viewers.length; i++) {
        if (tv.viewers[i].type === 'Scatter plot') { sp = tv.viewers[i]; break; }
      }
      return {
        spColor: sp?.props.colorColumnName,
        spMarker: sp?.props.markersColumnName,
      };
    });
    expect(result.spColor).toBe('Core');
    expect(result.spMarker).toBe('Core');
  });

  // Step 5: Scatterplot: Marker=Core, Color=Series — check structure rendering
  await softStep('Scatterplot: Marker=Core, Color=Series — check structure rendering', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      let sp = null;
      for (let i = 0; i < tv.viewers.length; i++) {
        if (tv.viewers[i].type === 'Scatter plot') { sp = tv.viewers[i]; break; }
      }
      sp.props.colorColumnName = 'Series';
      await new Promise(r => setTimeout(r, 2000));

      const spContainer = document.querySelector('[name="viewer-Scatter-plot"]');
      const extraItems = spContainer?.querySelectorAll('.d4-legend-item-extra') || [];
      let structureIcons = 0;
      for (const item of Array.from(extraItems)) {
        if (item.querySelector('.grok-icon, .svg-icon, img, canvas'))
          structureIcons++;
      }

      return {
        color: sp.props.colorColumnName,
        marker: sp.props.markersColumnName,
        markerLegendCount: extraItems.length,
        structureIcons,
      };
    });
    expect(result.color).toBe('Series');
    expect(result.marker).toBe('Core');
    expect(result.markerLegendCount).toBeGreaterThan(0);
    expect(result.structureIcons).toBeGreaterThan(0);
  });

  // Step 6: Save and apply layout
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
      try { df.columns.byName('Core').semType = 'Molecule'; } catch(_) {}
      try { df.columns.byName('Structure').semType = 'Molecule'; } catch(_) {}

      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      let sp = null;
      for (let i = 0; i < tv2.viewers.length; i++) {
        if (tv2.viewers[i].type === 'Scatter plot') { sp = tv2.viewers[i]; break; }
      }

      const res = {
        viewerCount: tv2.viewers.length,
        spColor: sp?.props.colorColumnName,
        spMarker: sp?.props.markersColumnName,
      };

      await grok.dapi.layouts.delete(saved);
      return res;
    });
    expect(result.viewerCount).toBe(8);
    expect(result.spColor).toBe('Series');
    expect(result.spMarker).toBe('Core');
  });

  // Step 7: Save and open project — verify legend and structure rendering
  await softStep('Save and open project — verify legend and structure rendering', async () => {
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
      try { df.columns.byName('Core').semType = 'Molecule'; } catch(_) {}
      try { df.columns.byName('Structure').semType = 'Molecule'; } catch(_) {}

      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      let sp = null;
      for (let i = 0; i < tv2.viewers.length; i++) {
        if (tv2.viewers[i].type === 'Scatter plot') { sp = tv2.viewers[i]; break; }
      }
      const spContainer = document.querySelector('[name="viewer-Scatter-plot"]');
      const extraItems = spContainer?.querySelectorAll('.d4-legend-item-extra') || [];
      let structureIcons = 0;
      for (const item of Array.from(extraItems)) {
        if (item.querySelector('.grok-icon, .svg-icon, img, canvas'))
          structureIcons++;
      }

      const res = {
        viewerCount: tv2.viewers.length,
        spColor: sp?.props.colorColumnName,
        spMarker: sp?.props.markersColumnName,
        structureIcons,
      };

      await grok.dapi.layouts.delete(saved);
      grok.shell.closeAll();
      return res;
    });
    expect(result.viewerCount).toBe(8);
    expect(result.spColor).toBe('Series');
    expect(result.spMarker).toBe('Core');
    expect(result.structureIcons).toBeGreaterThan(0);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
