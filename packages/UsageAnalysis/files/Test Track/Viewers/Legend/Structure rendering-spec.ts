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

async function openDataset(page: any) {
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell, {timeout: 15000});

  await page.evaluate(async (path: string) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
}

test('Structure rendering in legends', async ({page}) => {
  stepErrors.length = 0;

  await softStep('1-2. Open SPGI and add viewers', async () => {
    await openDataset(page);

    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 500));
      tv.addViewer('Histogram');
      await new Promise(r => setTimeout(r, 500));
      tv.addViewer('Line chart');
      await new Promise(r => setTimeout(r, 500));
      tv.addViewer('Bar chart');
      await new Promise(r => setTimeout(r, 500));
      tv.addViewer('Pie chart');
      await new Promise(r => setTimeout(r, 500));
      tv.addViewer('Trellis plot');
      await new Promise(r => setTimeout(r, 500));
      tv.addViewer('Box plot');
      await new Promise(r => setTimeout(r, 1000));

      return {count: Array.from(tv.viewers).length};
    });
    expect(result.count).toBe(8);
  });

  await softStep('3-4. Add structure legend (Core) to viewers', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const viewers = Array.from(tv.viewers);
      const sp = viewers.find((v: any) => v.type === 'Scatter plot') as any;
      const hist = viewers.find((v: any) => v.type === 'Histogram') as any;
      const bc = viewers.find((v: any) => v.type === 'Bar chart') as any;
      const bp = viewers.find((v: any) => v.type === 'Box plot') as any;
      const lc = viewers.find((v: any) => v.type === 'Line chart') as any;
      const tp = viewers.find((v: any) => v.type === 'Trellis plot') as any;
      const pie = viewers.find((v: any) => v.type === 'Pie chart') as any;

      sp.props.markersColumnName = 'Core';
      sp.props.colorColumnName = 'Core';

      try { hist.props.splitColumnName = 'Core'; } catch(e) {}
      try { bc.props.splitColumnName = 'Core'; } catch(e) {}
      try { bp.props.colorColumnName = 'Core'; } catch(e) {}
      try { pie.props.splitColumnName = 'Core'; } catch(e) {}
      try { lc.props.splitColumnName = 'Core'; } catch(e) {}
      try { tp.props.colorColumnName = 'Core'; } catch(e) {}

      await new Promise(r => setTimeout(r, 2000));

      return {
        spColor: sp.props.colorColumnName,
        spMarker: sp.props.markersColumnName,
      };
    });
    expect(result.spColor).toBe('Core');
    expect(result.spMarker).toBe('Core');
  });

  await softStep('5. Scatterplot: Marker=Core, Color=Series — check structure rendering', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const viewers = Array.from(tv.viewers);
      const sp = viewers.find((v: any) => v.type === 'Scatter plot') as any;

      sp.props.colorColumnName = 'Series';
      await new Promise(r => setTimeout(r, 1500));

      const spContainer = document.querySelector('[name="viewer-Scatter-plot"]');
      const extraItems = spContainer!.querySelectorAll('.d4-legend-item-extra');
      let structureIcons = 0;
      for (const item of extraItems)
        if (item.querySelector('.grok-icon, .svg-icon, img, canvas'))
          structureIcons++;

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

  await softStep('6. Save and apply layout', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const viewers = Array.from(tv2.viewers);
      const sp = viewers.find((v: any) => v.type === 'Scatter plot') as any;

      const res = {
        viewerCount: viewers.length,
        viewerTypes: viewers.map((v: any) => v.type),
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

  await softStep('7. Save and open project — verify legend and structure rendering', async () => {
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
        setTimeout(resolve, 3000);
      });

      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const viewers = Array.from(tv2.viewers);
      const sp = viewers.find((v: any) => v.type === 'Scatter plot') as any;
      const spContainer = document.querySelector('[name="viewer-Scatter-plot"]');
      const extraItems = spContainer!.querySelectorAll('.d4-legend-item-extra');
      let structureIcons = 0;
      for (const item of extraItems)
        if (item.querySelector('.grok-icon, .svg-icon, img, canvas'))
          structureIcons++;

      const res = {
        viewerCount: viewers.length,
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
