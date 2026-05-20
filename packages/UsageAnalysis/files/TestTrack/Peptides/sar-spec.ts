import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

test('SAR — Launch and verify viewers', async ({page}) => {
  // Multiple waitForFunctions (30s each) for SAR/MCL compute — won't fit in
  // the playwright default 60s per-test budget.
  test.setTimeout(300_000);
  await loginToDatagrok(page);

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
    const setup = await page!.evaluate(async () => {
      grok.shell.o = grok.shell.tv.dataFrame.col('AlignedSequence');
      await new Promise(r => setTimeout(r, 1500));

      const headers = document.querySelectorAll('.d4-accordion-pane-header');
      let expanded = false;
      for (const h of headers) {
        if (h.textContent?.trim().startsWith('Peptides') && !h.classList.contains('expanded')) {
          h.click();
          expanded = true;
        }
      }
      await new Promise(r => setTimeout(r, 1000));

      const scalingSelect = document.querySelector('[name="input-Scaling"]');
      const scalingFound = !!scalingSelect;
      if (scalingSelect) {
        (scalingSelect as any).value = '-lg';
        scalingSelect.dispatchEvent(new Event('change', {bubbles: true}));
      }
      await new Promise(r => setTimeout(r, 500));

      const launchBtn = document.querySelector('[name="button-Launch-SAR"]');
      const launchFound = !!launchBtn;
      if (launchBtn) launchBtn.click();
      return {expanded, scalingFound, launchFound};
    });
    // Surface the real failure if the SAR panel UI shape changed and our
    // selectors silently no-opped — previously this would just time out on
    // the Logo Summary Table wait below, hiding the actual cause.
    expect(setup.scalingFound, '[name="input-Scaling"] not found in SAR panel — UI shape may have changed').toBe(true);
    expect(setup.launchFound, '[name="button-Launch-SAR"] not found — UI shape may have changed').toBe(true);

    // Wait for the second SAR launch to materialize. Verified live against dev
    // 2026-05-21 (MCP recon): re-launching SAR with a different scaling does NOT
    // produce a 'Logo Summary Table' viewer on this build — it simply adds a
    // second `Sequence Variability Map` (and second `Most Potent Residues`) to
    // the workspace. The reliable invariant is therefore "SVM count went from 1
    // to 2" — not the name of a viewer the platform doesn't create.
    await page!.waitForFunction(() => {
      const svmCount = (Array.from(grok.shell.tv.viewers) as any[])
        .filter((v: any) => v.type === 'Sequence Variability Map').length;
      return svmCount >= 2;
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
