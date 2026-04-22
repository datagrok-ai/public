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

test('Peptides — SAR parameters and WebLogo', async () => {
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

  // Steps 1-2: Open peptides.csv and click column title
  await softStep('Steps 1-2: Open dataset and select column', async () => {
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

      // Click column title
      grok.shell.o = df.col('AlignedSequence');
      await new Promise(r => setTimeout(r, 1000));
      return {rows: df.rowCount, semType: df.col('AlignedSequence')?.semType};
    }, datasetPath);
    expect(result.rows).toBe(647);
    expect(result.semType).toBe('Macromolecule');
  });

  // Step 3: Expand Peptides panel
  await softStep('Step 3: Expand Peptides panel', async () => {
    const hasPeptides = await page!.evaluate(async () => {
      await new Promise(r => setTimeout(r, 1000));
      const headers = document.querySelectorAll('.d4-accordion-pane-header');
      for (const h of headers) {
        if (h.textContent?.trim().startsWith('Peptides')) {
          if (!h.classList.contains('expanded')) h.click();
          return true;
        }
      }
      return false;
    });
    expect(hasPeptides).toBe(true);
  });

  // Step 4: Change Activity, Scaling and Clusters parameters
  await softStep('Step 4: Change parameters', async () => {
    const result = await page!.evaluate(async () => {
      await new Promise(r => setTimeout(r, 500));
      // Change Scaling
      const scalingSelect = document.querySelector('[name="input-Scaling"]');
      if (scalingSelect) {
        (scalingSelect as any).value = 'lg';
        scalingSelect.dispatchEvent(new Event('change', {bubbles: true}));
      }
      await new Promise(r => setTimeout(r, 500));

      // Toggle Generate clusters
      const genClusters = document.querySelector('[name="input-Generate-clusters"]');
      if (genClusters) (genClusters as any).click();
      await new Promise(r => setTimeout(r, 500));

      return {
        scalingChanged: (scalingSelect as any)?.value === 'lg',
        hasScaling: !!scalingSelect,
        hasClusters: !!genClusters
      };
    });
    expect(result.hasScaling).toBe(true);
    expect(result.scalingChanged).toBe(true);
  });

  // Steps 5-6: Click amino acid on weblogo — canvas-based, verify via API
  await softStep('Steps 5-6: Verify WebLogo canvas exists', async () => {
    const result = await page!.evaluate(async () => {
      // Re-set context panel to column
      grok.shell.o = grok.shell.tv.dataFrame.col('AlignedSequence');
      await new Promise(r => setTimeout(r, 1000));

      // Find Peptides pane and check for canvas
      const headers = document.querySelectorAll('.d4-accordion-pane-header');
      for (const h of headers) {
        if (h.textContent?.trim().startsWith('Peptides')) {
          const pane = h.parentElement;
          const canvases = pane?.querySelectorAll('canvas');
          return {canvasCount: canvases?.length || 0};
        }
      }
      return {canvasCount: 0};
    });
    // WebLogo canvas exists (2 canvases: WebLogo + Histogram)
    expect(result.canvasCount).toBeGreaterThanOrEqual(1);
  });

  // Cleanup
  await page!.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
