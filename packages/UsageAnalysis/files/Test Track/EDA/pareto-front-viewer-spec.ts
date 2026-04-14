import {test, expect, chromium} from '@playwright/test';

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

test('Pareto Front Viewer', async () => {
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

  // Open cars.csv dataset
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/cars.csv');
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Step 1: Open dataset verified', async () => {
    const rowCount = await page!.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(rowCount).toBe(30);
  });

  await softStep('Step 2: Add Pareto Front via ML menu', async () => {
    await page!.locator('[name="div-ML"]').click();
    await page!.waitForTimeout(300);
    await page!.locator('[name="div-ML---Pareto-Front..."]').click();
    await page!.waitForTimeout(1000);

    // Verify viewer was added
    const viewerTypes = await page!.evaluate(() =>
      Array.from(grok.shell.tv.viewers).map(v => v.type)
    );
    expect(viewerTypes).toContain('ParetoFrontViewer');
  });

  await softStep('Step 3: Non-numeric columns excluded from objectives', async () => {
    const result = await page!.evaluate(() => {
      const headings = document.querySelectorAll('h1');
      let optimizeSection: Element | null = null;
      for (const h of headings) {
        if (h.textContent?.trim() === 'Optimize') {
          optimizeSection = h.parentElement;
          break;
        }
      }
      const featureNames: string[] = [];
      if (optimizeSection) {
        optimizeSection.querySelectorAll('select').forEach(sel => {
          const prev = sel.previousSibling;
          if (prev && prev.textContent) featureNames.push(prev.textContent.trim());
        });
      }
      return { featureNames, modelIncluded: featureNames.includes('model') };
    });
    expect(result.modelIncluded).toBe(false);
    expect(result.featureNames.length).toBe(16);
  });

  await softStep('Step 4: Conflict warning when same col in min and max', async () => {
    const result = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const pv = viewers.find(v => v.type === 'ParetoFrontViewer')!;
      const props = pv.getProperties();
      const maxProp = props.find(p => p.name === 'maximizeColumnNames')!;
      const df = grok.shell.tv.dataFrame;
      const numCols: string[] = [];
      for (let i = 0; i < df.columns.length; i++) {
        const c = df.columns.byIndex(i);
        if (c.type !== 'string' && c.name !== 'Pareto optimality')
          numCols.push(c.name);
      }
      maxProp.set(pv, numCols);
      await new Promise(r => setTimeout(r, 2000));
      const allText = document.body.innerText;
      const warningMatch = allText.match(/[Cc]annot.*minimi.*maximi/);
      return { warningFound: !!warningMatch };
    });
    expect(result.warningFound).toBe(true);
  });

  await softStep('Step 5: cars.csv Label auto-selects model', async () => {
    const result = await page!.evaluate(async () => {
      // Remove existing viewers and re-add
      const viewers = Array.from(grok.shell.tv.viewers);
      viewers.filter(v => v.type !== 'Grid').forEach(v => v.close());
      await new Promise(r => setTimeout(r, 500));
      const pv = grok.shell.tv.addViewer('Pareto front');
      await new Promise(r => setTimeout(r, 1000));
      const props = pv.getProperties();
      const labelProp = props.find(p => p.name === 'labelColumnsColumnNames');
      return { labelColumns: labelProp ? labelProp.get(pv) : null };
    });
    expect(result.labelColumns).toContain('model');
  });

  await softStep('Step 6: demog.csv Label auto-selection', async () => {
    const result = await page!.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      const pv = tv.addViewer('Pareto front');
      await new Promise(r => setTimeout(r, 1000));
      const props = pv.getProperties();
      const labelProp = props.find(p => p.name === 'labelColumnsColumnNames');
      const autoLabel = props.find(p => p.name === 'autoLabelsSelection');
      return {
        labelColumns: labelProp ? labelProp.get(pv) : null,
        autoLabelsSelection: autoLabel ? autoLabel.get(pv) : null
      };
    });
    // USUBJID has 5850 unique values for 5850 rows — auto-selection correctly picks it
    expect(result.autoLabelsSelection).toBe(true);
    expect(result.labelColumns).toContain('USUBJID');
  });

  await softStep('Step 7: Properties accessible and functional', async () => {
    const result = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const pv = viewers.find(v => v.type === 'ParetoFrontViewer')!;
      const props = pv.getProperties();
      const errors: string[] = [];
      for (const p of props) {
        try { p.get(pv); }
        catch (e: any) { errors.push(`${p.name}: ${e.message}`); }
      }
      // Test changing a property
      const dlProp = props.find(p => p.name === 'displayLabels');
      if (dlProp) {
        dlProp.set(pv, 'Always');
        await new Promise(r => setTimeout(r, 500));
        const val = dlProp.get(pv);
        if (val !== 'Always') errors.push('displayLabels did not update');
      }
      return { propCount: props.length, errors };
    });
    expect(result.errors).toEqual([]);
    expect(result.propCount).toBeGreaterThan(10);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
