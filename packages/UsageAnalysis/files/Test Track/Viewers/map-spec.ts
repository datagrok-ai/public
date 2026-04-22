import {test, expect, chromium} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

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

test('Map Viewer: add map, set color/size/renderType', async () => {
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

  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;

    // Use earthquakes.csv (MYmeteoritesTest.csv not found)
    const df = await grok.dapi.files.readCsv('System:DemoFiles/geo/earthquakes.csv');
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  });

  // Steps 1-2: Add Map viewer
  await softStep('Steps 1-2: Add Map viewer', async () => {
    const result = await page!.evaluate(async () => {
      const tv = grok.shell.tv;
      const mv = tv.addViewer('Map');
      await new Promise(r => setTimeout(r, 3000));
      return {type: mv?.type, viewerCount: Array.from(tv.viewers).length};
    });
    expect(result.type).toBe('Map');
  });

  // Steps 5-6: Set color and size
  await softStep('Steps 5-6: Set color and size columns', async () => {
    const result = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const mv = viewers.find((v: any) => v.type === 'Map');
      if (!mv) return {error: 'not found'};
      mv.setOptions({colorColumnName: 'mag', sizeColumnName: 'depth'});
      await new Promise(r => setTimeout(r, 1000));
      const opts = mv.getOptions();
      return {color: opts.look?.colorColumnName, size: opts.look?.sizeColumnName};
    });
    expect(result.color).toBe('mag');
    expect(result.size).toBe('depth');
  });

  // Steps 8, 13: Marker size and render type
  await softStep('Steps 8, 13: Marker min size and render type', async () => {
    const result = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const mv = viewers.find((v: any) => v.type === 'Map');
      if (!mv) return {error: 'not found'};

      mv.setOptions({markerMinSize: 5});
      mv.setOptions({renderType: 'heatmap'});
      await new Promise(r => setTimeout(r, 1000));
      mv.setOptions({renderType: 'markers'});
      await new Promise(r => setTimeout(r, 1000));

      const opts = mv.getOptions();
      return {minSize: opts.look?.markerMinSize, renderType: opts.look?.renderType};
    });
    expect(result.minSize).toBe(5);
    expect(result.renderType).toBe('markers');
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
