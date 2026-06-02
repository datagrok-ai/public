/* ---
sub_features_covered: [peptides.workflow.sar-dialog, peptides.workflow.analyze-ui, peptides.workflow.start-analysis, peptides.model.add-sequence-space, peptides.viewers.cluster-max-activity, peptides.viewers.logo-summary-table, peptides.widgets.settings-dialog, peptides.compute.calculate-cluster-statistics]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
test.use(specTestOptions);
const datasetPath = 'System:DemoFiles/bio/peptides.csv';
test('Peptide Space — top-menu SAR launch with sequence-space + MCL clustering', async ({page}) => {
  test.setTimeout(300_000);
  await loginToDatagrok(page);
  await softStep('Setup: open the peptides Macromolecule table', async () => {
    const result = await page.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = false;
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 4000);
      });
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 4000));
      try { await grok.functions.call('Peptides:initPeptides'); }
      catch (e) { console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', String(e)); }
      return {
        rows: df.rowCount,
        semType: df.col('AlignedSequence')?.semType ?? null,
      };
    }, datasetPath);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
    expect(result.rows, 'peptides.csv should load 647 rows').toBe(647);
    expect(result.semType, 'AlignedSequence must be a Macromolecule column').toBe('Macromolecule');
  });
  await softStep('Scenario 1 (step 1): launch SAR via Bio | Analyze | SAR... top menu', async () => {
    const opened = await page.evaluate(async () => {
      const bio = document.querySelector('[name="div-Bio"]') as HTMLElement | null;
      const bioVisible = bio ? bio.offsetParent !== null : false;
      if (bio) bio.click();
      await new Promise((r) => setTimeout(r, 700));
      const analyze = document.querySelector('[name="div-Bio---Analyze"]') as HTMLElement | null;
      if (analyze) {
        analyze.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        analyze.dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
      }
      await new Promise((r) => setTimeout(r, 700));
      const sar = document.querySelector('[name="div-Bio---Analyze---SAR..."]') as HTMLElement | null;
      if (sar) sar.click();
      await new Promise((r) => setTimeout(r, 2500));
      const dlg = document.querySelector('[name="dialog-Analyze-Peptides"]');
      return {
        bioFound: !!bio,
        bioVisible,
        analyzeFound: !!analyze,
        sarFound: !!sar,
        dialogFound: !!dlg,
      };
    });
    expect(opened.bioFound, '[name="div-Bio"] top-menu entry not found').toBe(true);
    expect(opened.bioVisible,
      '[name="div-Bio"] present but not visible (offsetParent null — menu in overflow at this viewport)').toBe(true);
    expect(opened.analyzeFound, '[name="div-Bio---Analyze"] submenu entry not found').toBe(true);
    expect(opened.sarFound, '[name="div-Bio---Analyze---SAR..."] menu item not found').toBe(true);
    expect(opened.dialogFound, '[name="dialog-Analyze-Peptides"] config dialog did not open').toBe(true);
  });
  await softStep('Scenario 1 (step 2): accept default config, run SAR, verify viewers', async () => {
    await page.evaluate(async () => {
      const dlg = document.querySelector('[name="dialog-Analyze-Peptides"]');
      const ok = (dlg?.querySelector('[name="button-OK"]')
        ?? document.querySelector('[name="button-OK"]')) as HTMLElement | null;
      if (ok) ok.click();
    });
    await page.waitForFunction(() => {
      return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
    }, {timeout: 90000});
    await page.waitForTimeout(8000);
    const viewers = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      return Array.from(tv.viewers).map((v) => v.type);
    });
    expect(viewers, 'Sequence Variability Map must attach after SAR launch').toContain('Sequence Variability Map');
    expect(viewers, 'Most Potent Residues must attach after SAR launch').toContain('Most Potent Residues');
    expect(viewers, 'MCL clustering viewer must attach after SAR launch').toContain('MCL');
    if (!viewers.includes('Logo Summary Table'))
      console.log('[note] Logo Summary Table not in default top-menu SAR attach set on this build (cluster-/settings-dependent).');
  });
  await softStep('Scenario 2 (step 3): open the SAR settings dialog via the wrench', async () => {
    const opened = await page.evaluate(async () => {
      const wrench = document.querySelector(
        'i.grok-icon.fa-wrench[aria-label="Peptides analysis settings"]') as HTMLElement | null;
      const wrenchFound = !!wrench;
      if (wrench) wrench.click();
      await new Promise((r) => setTimeout(r, 2000));
      const dlg = document.querySelector('[name="dialog-Peptides-settings"]');
      const panes = dlg
        ? Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
          .map((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim())
        : [];
      return {wrenchFound, dlgFound: !!dlg, panes};
    });
    expect(opened.wrenchFound, 'Peptides analysis settings wrench not found on the SAR toolbar').toBe(true);
    expect(opened.dlgFound, '[name="dialog-Peptides-settings"] did not open').toBe(true);
    expect(opened.panes, 'settings dialog must expose the MCL accordion pane').toContain('MCL');
  });
  await softStep('Scenario 2 (step 4): change a representative MCL parameter + OK', async () => {
    const changed = await page.evaluate(async () => {
      const dlg = document.querySelector('[name="dialog-Peptides-settings"]')!;
      const tvBefore = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']);
      const inflationBefore = tvBefore
        ? (tvBefore.dataFrame.temp['peptidesModel'] as any)?._settings?.mclSettings?.inflation
        : null;
      const mclPane = Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
        .find((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim() === 'MCL');
      if (mclPane) {
        const h = mclPane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        if (h && !mclPane.classList.contains('expanded')) h.click();
      }
      await new Promise((r) => setTimeout(r, 800));
      const host = dlg.querySelector('[name="input-host-Inflation-Factor"]');
      const inner = host ? (host.querySelector('input') as HTMLInputElement | null) : null;
      const inputFound = !!inner;
      const prevVal = inner ? inner.value : null;
      if (inner) {
        inner.focus();
        inner.value = '2.5';
        inner.dispatchEvent(new Event('input', {bubbles: true}));
        inner.dispatchEvent(new Event('change', {bubbles: true}));
        inner.blur();
      }
      await new Promise((r) => setTimeout(r, 600));
      const newVal = inner ? inner.value : null;
      const ok = dlg.querySelector('[name="button-OK"]') as HTMLElement | null;
      if (ok) ok.click();
      await new Promise((r) => setTimeout(r, 2000));
      const dlgClosed = !document.querySelector('[name="dialog-Peptides-settings"]');
      return {inputFound, prevVal, newVal, inflationBefore, dlgClosed};
    });
    expect(changed.inputFound, '[name="input-host-Inflation-Factor"] not found in the MCL pane').toBe(true);
    expect(changed.newVal, 'MCL Inflation Factor did not accept the new value').toBe('2.5');
    expect(changed.inflationBefore, 'pre-change MCL inflation should be the default 1.4').toBe(1.4);
    expect(changed.dlgClosed, 'settings dialog did not close after OK').toBe(true);
  });
  await softStep('Scenario 2 (step 5): verify the MCL Viewer re-renders after the settings change', async () => {
    await page.waitForTimeout(12000);
    const state = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const viewers = Array.from(tv.viewers).map((v) => v.type);
      const model = tv.dataFrame.temp['peptidesModel'] as any;
      const inflationAfter = model?._settings?.mclSettings?.inflation ?? null;
      const mclDom = document.querySelector('[name="viewer-MCL"]');
      let mclHasRender = false;
      if (mclDom) {
        const cs = Array.from(mclDom.querySelectorAll('canvas')) as HTMLCanvasElement[];
        let best: HTMLCanvasElement | null = null, area = 0;
        for (const c of cs) { const r = c.getBoundingClientRect(); if (r.width * r.height > area) { area = r.width * r.height; best = c; } }
        if (best && best.width > 0 && best.height > 0) {
          try {
            const ctx = best.getContext('2d')!;
            const data = ctx.getImageData(0, 0, best.width, best.height).data;
            let nonBg = 0;
            for (let i = 0; i < data.length; i += 41)
              if (data[i + 3] > 0 && (data[i] < 250 || data[i + 1] < 250 || data[i + 2] < 250)) nonBg++;
            mclHasRender = nonBg > 50;
          } catch (e) { mclHasRender = false; }
        }
      }
      const lastError = grok.shell.lastError ? String(grok.shell.lastError) : '';
      return {viewers, modelPresent: !!model, inflationAfter, mclHasRender, lastError};
    });
    expect(state.inflationAfter, 'MCL Inflation Factor change did not propagate to the model').toBe(2.5);
    expect(state.modelPresent, 'PeptidesModel cache lost after the MCL settings change').toBe(true);
    expect(state.viewers, 'MCL viewer must persist after the settings change').toContain('MCL');
    expect(state.mclHasRender, 'MCL viewer did not re-render visible content after the settings change').toBe(true);
    expect(/setTrue|null/.test(state.lastError),
      `GROK-19145 invariant: post-OK MCL compute produced a null-receiver error: ${state.lastError}`).toBe(false);
  });
  // Cleanup.
  await page.evaluate(() => grok.shell.closeAll());
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
