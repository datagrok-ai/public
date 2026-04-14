import {test, expect, chromium} from '@playwright/test';

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

test('SAR — Launch and verify viewers', async () => {
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

  // Steps 1-4: Open peptides, select column, launch SAR
  await softStep('Steps 1-4: Open peptides and launch SAR', async () => {
    await page!.evaluate(async (path) => {
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

      // Select column and expand Peptides panel
      grok.shell.o = df.col('AlignedSequence');
      await new Promise(r => setTimeout(r, 1500));

      const headers = document.querySelectorAll('.d4-accordion-pane-header');
      for (const h of headers) {
        if (h.textContent?.trim().startsWith('Peptides') && !h.classList.contains('expanded'))
          h.click();
      }
      await new Promise(r => setTimeout(r, 1000));

      // Click Launch SAR
      const launchBtn = document.querySelector('[name="button-Launch-SAR"]');
      if (launchBtn) launchBtn.click();
    }, datasetPath);

    // Wait for SAR viewers
    await page!.waitForFunction(() => {
      const viewers = Array.from(grok.shell.tv.viewers);
      return viewers.some(v => v.type === 'Sequence Variability Map');
    }, {timeout: 30000});
  });

  // Step 5: Verify viewers appeared
  await softStep('Step 5: Verify SAR viewers', async () => {
    await page!.waitForTimeout(3000);
    const viewers = await page!.evaluate(() => {
      return Array.from(grok.shell.tv.viewers).map(v => v.type);
    });
    expect(viewers).toContain('Sequence Variability Map');
    expect(viewers).toContain('Most Potent Residues');
    expect(viewers).toContain('MCL');
  });

  // Steps 6-9: Change parameters and re-launch SAR
  await softStep('Steps 6-9: Change parameters and re-launch', async () => {
    await page!.evaluate(async () => {
      grok.shell.o = grok.shell.tv.dataFrame.col('AlignedSequence');
      await new Promise(r => setTimeout(r, 1500));

      const headers = document.querySelectorAll('.d4-accordion-pane-header');
      for (const h of headers) {
        if (h.textContent?.trim().startsWith('Peptides') && !h.classList.contains('expanded'))
          h.click();
      }
      await new Promise(r => setTimeout(r, 1000));

      const scalingSelect = document.querySelector('[name="input-Scaling"]');
      if (scalingSelect) {
        (scalingSelect as any).value = '-lg';
        scalingSelect.dispatchEvent(new Event('change', {bubbles: true}));
      }
      await new Promise(r => setTimeout(r, 500));

      const launchBtn = document.querySelector('[name="button-Launch-SAR"]');
      if (launchBtn) launchBtn.click();
    });

    // Wait for Logo Summary Table to appear (indicates second SAR completed)
    await page!.waitForFunction(() => {
      const viewers = Array.from(grok.shell.tv.viewers);
      return viewers.some(v => v.type === 'Logo Summary Table');
    }, {timeout: 30000});
  });

  // Step 10: Switch between Mutation Cliffs and Invariant Map
  await softStep('Step 10: Toggle Mutation Cliffs / Invariant Map', async () => {
    // Switch to Invariant Map
    await page!.evaluate(async () => {
      const svms = document.querySelectorAll('[name="viewer-Sequence-Variability-Map"]');
      // Use the latest (second) SVM viewer
      const svm = svms[svms.length - 1];
      const labels = svm.querySelectorAll('.ui-input-bool');
      for (const label of labels) {
        if (label.textContent?.includes('Invariant Map')) {
          const input = label.querySelector('input');
          if (input) input.click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 1000));
    });

    // Switch back to Mutation Cliffs
    await page!.evaluate(async () => {
      const svms = document.querySelectorAll('[name="viewer-Sequence-Variability-Map"]');
      const svm = svms[svms.length - 1];
      const labels = svm.querySelectorAll('.ui-input-bool');
      for (const label of labels) {
        if (label.textContent?.includes('Mutation Cliffs')) {
          const input = label.querySelector('input');
          if (input) input.click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 1000));
    });
  });

  // Cleanup
  await page!.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
