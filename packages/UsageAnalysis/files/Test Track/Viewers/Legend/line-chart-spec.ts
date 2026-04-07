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

test('Line chart legend', async () => {
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
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Step 1-2: Add line chart
  await softStep('Add line chart', async () => {
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-line-chart"]') as HTMLElement;
      icon.click();
    });
    await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 10000});
    const viewerCount = await page.evaluate(() => grok.shell.tv.viewers.length);
    expect(viewerCount).toBe(2);
  });

  // Step 3: Set Split to Primary Series Name and Series, verify diverse colors
  await softStep('Set Split and verify color diversity', async () => {
    await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart')
          viewers[i].props.splitColumnNames = ['Primary Series Name', 'Series'];
      }
    });
    await page.waitForTimeout(1500);
    const splitNames = await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart')
          return viewers[i].props.splitColumnNames;
      }
    });
    expect(splitNames).toEqual(['Primary Series Name', 'Series']);
  });

  // Step 4: Enable Multi Axis, verify per-line legend categories
  await softStep('Multi Axis legend categories', async () => {
    await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart')
          viewers[i].props.multiAxis = true;
      }
    });
    await page.waitForTimeout(1500);
    const multiAxis = await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart')
          return viewers[i].props.multiAxis;
      }
    });
    expect(multiAxis).toBe(true);
  });

  // Step 5: Save and apply layout (round 1)
  await softStep('Save and apply layout (round 1)', async () => {
    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id;
    });
    await page.waitForTimeout(1000);

    await page.evaluate(() => grok.shell.tv.resetLayout());
    await page.waitForTimeout(2000);
    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
    }, layoutId);
    await page.waitForTimeout(3000);

    const restored = await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart')
          return {splitNames: viewers[i].props.splitColumnNames, multiAxis: viewers[i].props.multiAxis};
      }
    });
    expect(restored?.splitNames).toEqual(['Primary Series Name', 'Series']);
    expect(restored?.multiAxis).toBe(true);

    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      await grok.dapi.layouts.delete(saved);
    }, layoutId);
  });

  // Step 6: Save and open project (round 1)
  await softStep('Save and open project (round 1)', async () => {
    await page.locator('button:has-text("SAVE")').click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(3000);

    const projectId = await page.evaluate(() => grok.shell.project.id);

    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1000);
    await page.evaluate(async (id) => {
      const p = await grok.dapi.projects.find(id);
      await p.open();
    }, projectId);
    await page.waitForTimeout(6000);

    const restored = await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart')
          return {splitNames: viewers[i].props.splitColumnNames, multiAxis: viewers[i].props.multiAxis};
      }
    });
    expect(restored?.splitNames).toEqual(['Primary Series Name', 'Series']);
    expect(restored?.multiAxis).toBe(true);
  });

  // Step 7: Configure two Y columns
  await softStep('Configure two Y columns', async () => {
    await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart') {
          viewers[i].props.multiAxis = false;
          viewers[i].props.splitColumnNames = [];
          viewers[i].props.yColumnNames = ['CAST Idea ID', 'Average Mass'];
        }
      }
    });
    await page.waitForTimeout(1000);
    const yColumns = await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart')
          return viewers[i].props.yColumnNames;
      }
    });
    expect(yColumns).toEqual(['CAST Idea ID', 'Average Mass']);
  });

  // Step 8: Change Y column via in-plot selector, verify legend updates
  await softStep('Change Y column via in-plot selector', async () => {
    await page.evaluate(() => {
      const lc = document.querySelector('[name="viewer-Line-chart"]');
      const combos = lc!.querySelectorAll('[name="div-column-combobox-"]');
      const colLabel = combos[1]?.querySelector('.d4-column-selector-column') as HTMLElement;
      const r = colLabel.getBoundingClientRect();
      colLabel.dispatchEvent(new MouseEvent('mousedown', {clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, bubbles: true, button: 0}));
    });

    await page.waitForSelector('.d4-column-selector-backdrop', {timeout: 5000});
    const searchInput = page.locator('.d4-column-selector-search-input');
    await searchInput.waitFor({timeout: 3000});
    await page.evaluate(() => {
      const input = document.querySelector('.d4-column-selector-search-input') as HTMLInputElement;
      input.focus();
      input.value = 'TPSA';
      input.dispatchEvent(new Event('input', {bubbles: true}));
    });
    await page.waitForTimeout(500);

    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(1000);

    const yColumns = await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart')
          return viewers[i].props.yColumnNames;
      }
    });
    expect(yColumns).toEqual(['CAST Idea ID', 'TPSA']);
  });

  // Step 9: Save and apply layout (round 2)
  await softStep('Save and apply layout (round 2)', async () => {
    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id;
    });
    await page.waitForTimeout(1000);

    await page.evaluate(() => grok.shell.tv.resetLayout());
    await page.waitForTimeout(2000);
    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
    }, layoutId);
    await page.waitForTimeout(3000);

    const restored = await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart')
          return {yColumns: viewers[i].props.yColumnNames};
      }
    });
    expect(restored?.yColumns).toEqual(['CAST Idea ID', 'TPSA']);

    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      await grok.dapi.layouts.delete(saved);
    }, layoutId);
  });

  // Step 10: Save and open project (round 2)
  await softStep('Save and open project (round 2)', async () => {
    await page.locator('button:has-text("SAVE")').click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(3000);

    const projectId = await page.evaluate(() => grok.shell.project.id);

    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(1000);
    await page.evaluate(async (id) => {
      const p = await grok.dapi.projects.find(id);
      await p.open();
    }, projectId);
    await page.waitForTimeout(6000);

    const restored = await page.evaluate(() => {
      const viewers = grok.shell.tv.viewers;
      for (let i = 0; i < viewers.length; i++) {
        if (viewers[i].type === 'Line chart')
          return {yColumns: viewers[i].props.yColumnNames};
      }
    });
    expect(restored?.yColumns).toEqual(['CAST Idea ID', 'TPSA']);

    // Cleanup
    await page.evaluate(async (id) => {
      try { const p = await grok.dapi.projects.find(id); await grok.dapi.projects.delete(p); } catch(e) {}
    }, projectId);
  });

  // Summary check
  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
