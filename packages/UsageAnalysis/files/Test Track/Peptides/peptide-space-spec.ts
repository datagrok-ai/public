import {test, expect, chromium} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/bio/peptides.csv';
const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Peptide Space — SAR and MCL', async () => {
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

  // Setup: Open peptides.csv
  await softStep('Setup: Open peptides.csv', async () => {
    const result = await page!.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach(d => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = false;

      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      const hasBio = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
        .some(c => c.semType === 'Macromolecule');
      if (hasBio) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise(r => setTimeout(r, 200));
        }
        await new Promise(r => setTimeout(r, 5000));
      }
      return {rows: df.rowCount, semType: df.col('AlignedSequence')?.semType};
    }, datasetPath);
    expect(result.rows).toBe(647);
    expect(result.semType).toBe('Macromolecule');
  });

  // Step 1: Launch SAR from Bio > Analyze > SAR
  await softStep('Step 1: Launch SAR via Bio menu', async () => {
    await page!.evaluate(async () => {
      const bio = document.querySelector('[name="div-Bio"]');
      bio.click();
      await new Promise(r => setTimeout(r, 500));
      const analyze = document.querySelector('[name="div-Bio---Analyze"]');
      analyze.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      analyze.dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const sar = document.querySelector('[name="div-Bio---Analyze---SAR..."]');
      sar.click();
      await new Promise(r => setTimeout(r, 2000));
    });
    // Dialog should appear
    const hasDialog = await page!.evaluate(() => !!document.querySelector('.d4-dialog'));
    expect(hasDialog).toBe(true);
  });

  // Step 2: Click OK and wait for calculation results
  await softStep('Step 2: Run SAR and wait for results', async () => {
    await page!.evaluate(async () => {
      const ok = document.querySelector('[name="button-OK"]');
      if (ok) ok.click();
    });
    // Wait for MCL viewer to appear (up to 30s)
    await page!.waitForFunction(() => {
      const viewers = Array.from(grok.shell.tv.viewers);
      return viewers.some(v => v.type === 'MCL');
    }, {timeout: 30000});

    const viewers = await page!.evaluate(() => {
      return Array.from(grok.shell.tv.viewers).map(v => v.type);
    });
    expect(viewers).toContain('MCL');
    expect(viewers).toContain('Most Potent Residues');
    expect(viewers).toContain('Sequence Variability Map');
  });

  // Step 3-5: Wait for MCL clustering to complete and verify results
  await softStep('Step 3-5: Verify MCL clustering results', async () => {
    // Wait for Cluster (MCL) column to appear (MCL runs async after viewer creation)
    await page!.waitForFunction(() => {
      return !!grok.shell.tv.dataFrame.col('Cluster (MCL)');
    }, {timeout: 30000});

    const result = await page!.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const clusterCol = df.col('Cluster (MCL)');
      const uniqueClusters = new Set();
      for (let i = 0; i < df.rowCount; i++) {
        const v = clusterCol.get(i);
        if (v !== null && v !== undefined) uniqueClusters.add(v);
      }
      return {clusters: uniqueClusters.size, colCount: df.columns.length};
    });
    expect(result.clusters).toBeGreaterThanOrEqual(1);
  });

  // Cleanup
  await page!.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
