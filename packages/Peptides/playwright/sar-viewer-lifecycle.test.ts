/* ---
sub_features_covered: [peptides.model.add-cluster-max-activity, peptides.model.add-dendrogram, peptides.model.add-logo-summary-table, peptides.model.add-monomer-position, peptides.model.add-most-potent-residues, peptides.model.viewer-type, peptides.widgets.settings-dialog, peptides.workflow.start-analysis]
--- */
// SAR viewer lifecycle — top-menu Bio | Analyze | SAR... launch, verify model.add-* surfaces realize,
// and round-trip the Settings-dialog Viewers-pane toggles (Dendrogram + Active peptide selection).
// Cluster max activity has no UI entry point and is driven through the model; Dendrogram goes
// through its Settings toggle, which is the only path that proves the viewer actually mounts.

import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import {waitForViewers} from './helpers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

const VIEWER_TYPE = {
  SEQUENCE_VARIABILITY_MAP: 'Sequence Variability Map',
  MOST_POTENT_RESIDUES: 'Most Potent Residues',
  LOGO_SUMMARY_TABLE: 'Logo Summary Table',
  DENDROGRAM: 'Dendrogram',
  CLUSTER_MAX_ACTIVITY: 'Active peptide selection',
  MCL: 'MCL',
};

const DEFAULT_ATTACH_VIEWERS = [
  VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP,
  VIEWER_TYPE.MOST_POTENT_RESIDUES,
  VIEWER_TYPE.MCL,
  VIEWER_TYPE.LOGO_SUMMARY_TABLE,
];

/** Open the Settings dialog from the ribbon wrench and expand its Viewers pane. */
async function openViewersPane(page: Page): Promise<void> {
  await page.evaluate(() => {
    const wrench = document.querySelector(
      'i.grok-icon.fa-wrench[aria-label="Peptides analysis settings"]') as HTMLElement | null;
    if (wrench) wrench.click();
  });
  await page.locator('[name="dialog-Peptides-settings"]').waitFor({timeout: 8000});
  await page.evaluate(async () => {
    const dlg = document.querySelector('[name="dialog-Peptides-settings"]')!;
    // d4 panes signal collapse via the content's display, not a class.
    const pane = Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
      .find((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim() === 'Viewers');
    const content = pane?.querySelector('.d4-accordion-pane-content') as HTMLElement | null;
    if (content && getComputedStyle(content).display === 'none') {
      (pane!.querySelector('.d4-accordion-pane-header') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 600));
    }
  });
}

test('SAR viewer lifecycle — model.add-* family + VIEWER_TYPE discriminator + Settings dialog Viewers-pane round-trip', async ({page}) => {
  test.setTimeout(360_000);
  await loginToDatagrok(page);

  await softStep('Setup: open peptides dataset, prewarm Peptides:initPeptides', async () => {
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
      // Macromolecule dataset: wait for grid canvas + Bio package settle.
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 4000));

      // GROK-17557 prewarm — also loads the Bio TreeHelper for the model.add-dendrogram path.
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

  await softStep('Scenario 1 (steps 1-3): launch SAR via Bio | Analyze | SAR... top menu', async () => {
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

    await page.evaluate(async () => {
      const dlg = document.querySelector('[name="dialog-Analyze-Peptides"]');
      const ok = (dlg?.querySelector('[name="button-OK"]')
        ?? document.querySelector('[name="button-OK"]')) as HTMLElement | null;
      if (ok) ok.click();
    });
    await page.waitForFunction(() => {
      return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
    }, {timeout: 90000});
  });

  await softStep('Scenario 1 (step 4): verify deterministic default-attach (SVM + MPR + MCL + LST)', async () => {
    const viewers = await waitForViewers(page, DEFAULT_ATTACH_VIEWERS);
    for (const type of DEFAULT_ATTACH_VIEWERS)
      expect(viewers, `${type} must attach after Launch SAR`).toContain(type);

    const state = await page.evaluate((VT) => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = tv.dataFrame.temp['peptidesModel'];
      const viewers = Array.from(tv.viewers).map((v) => v.type);
      return {
        viewers,
        modelPresent: !!model,
        findSvm: !!model.findViewer(VT.SEQUENCE_VARIABILITY_MAP),
        findMpr: !!model.findViewer(VT.MOST_POTENT_RESIDUES),
        findMcl: !!model.findViewer(VT.MCL),
      };
    }, VIEWER_TYPE);
    expect(state.modelPresent, 'peptides.model PeptidesModel singleton must attach after Launch SAR').toBe(true);
    expect(state.viewers, 'Sequence Variability Map (peptides.model.add-monomer-position) must attach').toContain(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP);
    expect(state.findSvm, 'findViewer(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP) must return non-null').toBe(true);
    expect(state.viewers, 'Most Potent Residues (peptides.model.add-most-potent-residues) must attach').toContain(VIEWER_TYPE.MOST_POTENT_RESIDUES);
    expect(state.findMpr, 'findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) must return non-null').toBe(true);
    expect(state.viewers, 'MCL clustering viewer must attach').toContain(VIEWER_TYPE.MCL);
    expect(state.findMcl, 'findViewer(VIEWER_TYPE.MCL) must return non-null').toBe(true);
  });

  // Cluster max activity has no launch-dialog input, so it is driven through the model directly.
  await softStep('Scenario 1 (step 4): exercise peptides.model.add-cluster-max-activity', async () => {
    const result = await page.evaluate(async (VT) => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = tv.dataFrame.temp['peptidesModel'];
      await model.addClusterMaxActivityViewer();
      await new Promise((r) => setTimeout(r, 6000));
      return {
        viewers: Array.from(tv.viewers).map((v) => v.type),
        findCma: !!model.findViewer(VT.CLUSTER_MAX_ACTIVITY),
      };
    }, VIEWER_TYPE);

    expect(result.viewers, 'Active peptide selection viewer (peptides.model.add-cluster-max-activity) must attach')
      .toContain(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY);
    expect(result.findCma, 'findViewer(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY) must return non-null').toBe(true);
  });

  // Dendrogram ships OFF by default. Driving it through the Settings toggle is the path a user
  // has; addDendrogram() swallows its own failures into the logger (model.ts#L1094), so calling
  // it directly can only ever report "did not throw" and proves nothing about the viewer.
  await softStep('Scenario 1 (step 5): turn Dendrogram ON via the Settings dialog', async () => {
    await openViewersPane(page);
    const cb = page.locator('[name="dialog-Peptides-settings"] [name="input-Dendrogram"]');
    await cb.waitFor({timeout: 8000});
    expect(await cb.isEnabled(),
      'Dendrogram toggle is disabled — TreeHelper did not load (settings.ts#L139)').toBe(true);
    // Real click: a dispatched MouseEvent flips the checkbox without the platform ever
    // seeing the input, so the setting is applied but no viewer is ever added.
    if (!await cb.isChecked())
      await cb.click();
    expect(await cb.isChecked(), 'Dendrogram checkbox should be ON before OK').toBe(true);
    await page.locator('[name="dialog-Peptides-settings"] [name="button-OK"]').click();

    await page.waitForFunction(() =>
      !document.querySelector('[name="dialog-Peptides-settings"]'), null, {timeout: 8000});

    // The dendrogram is not a viewer: Dendrogram's hierarchicalClusteringUI attaches it as a
    // GridNeighbor beside the grid and parks it in grid.temp, so tv.viewers and findViewer
    // never see it however long you wait for them.
    const attached = await page.waitForFunction(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      return !!tv.grid?.temp?.['__dendrogram_neighbor_temp__'];
    }, null, {timeout: 120_000}).then(() => true).catch(() => false);
    expect(attached,
      'Dendrogram did not attach to the grid once its Settings toggle was turned ON').toBe(true);
  });

  await softStep('Scenario 2 (steps 1-2): open Settings dialog, confirm Viewers-pane structure', async () => {
    const opened = await page.evaluate(async () => {
      const wrench = document.querySelector(
        'i.grok-icon.fa-wrench[aria-label="Peptides analysis settings"]') as HTMLElement | null;
      const wrenchFound = !!wrench;
      if (wrench) wrench.click();
      await new Promise((r) => setTimeout(r, 2000));
      const dlg = document.querySelector('[name="dialog-Peptides-settings"]');
      const dlgFound = !!dlg;
      const panes = dlg
        ? Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
          .map((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim())
        : [];
      // d4 panes signal collapse via .d4-accordion-pane-content display, not a class; click header only if hidden.
      const viewersPane = dlg
        ? Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
          .find((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim() === 'Viewers')
        : null;
      if (viewersPane) {
        const content = viewersPane.querySelector('.d4-accordion-pane-content') as HTMLElement | null;
        const isCollapsed = !!content && getComputedStyle(content).display === 'none';
        if (isCollapsed) {
          const h = viewersPane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
          if (h) h.click();
          await new Promise((r) => setTimeout(r, 600));
        }
      }
      const toggles = viewersPane
        ? Array.from(viewersPane.querySelectorAll('input[type="checkbox"]'))
          .map((el: any) => ({name: el.getAttribute('name'), disabled: el.disabled}))
        : [];
      return {wrenchFound, dlgFound, panes, toggles};
    });
    expect(opened.wrenchFound, 'Peptides analysis settings wrench not found on the SAR toolbar').toBe(true);
    expect(opened.dlgFound, '[name="dialog-Peptides-settings"] did not open').toBe(true);
    expect(opened.panes, 'Settings dialog must expose the General pane').toContain('General');
    expect(opened.panes, 'Settings dialog must expose the Viewers pane').toContain('Viewers');
    expect(opened.panes, 'Settings dialog must expose the MCL pane').toContain('MCL');
    const toggleNames = opened.toggles.map((t: any) => t.name);
    expect(toggleNames, 'Viewers pane must expose [name="input-Dendrogram"]')
      .toContain('input-Dendrogram');
    expect(toggleNames, 'Viewers pane must expose [name="input-Sequence-space"]')
      .toContain('input-Sequence-space');
    expect(toggleNames, 'Viewers pane must expose [name="input-Active-peptide-selection"]')
      .toContain('input-Active-peptide-selection');
    // Cancel; the dialog is re-opened per step to drive a clean settingsChanged dispatch.
    await page.evaluate(() => {
      const cancel = document.querySelector('[name="dialog-Peptides-settings"] [name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    await page.waitForTimeout(500);
  });

  // Toggle CMA OFF: expand pane only if collapsed (no .expanded class), drive cb/OK via composed MouseEvents.
  await softStep('Scenario 2 (steps 3-4): toggle Cluster max activity OFF via the Settings dialog, apply', async () => {
    const before = await page.evaluate((VT) => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = tv.dataFrame.temp['peptidesModel'];
      return {
        cmaAttached: !!model.findViewer(VT.CLUSTER_MAX_ACTIVITY),
        viewers: Array.from(tv.viewers).map((v) => v.type),
      };
    }, VIEWER_TYPE);
    // Sanity: CMA is currently attached (from Scenario 1's direct invocation).
    expect(before.cmaAttached, 'precondition: CMA viewer should be attached entering Scenario 2 step 3').toBe(true);

    // Re-open the Settings dialog via the wrench.
    await page.evaluate(async () => {
      const wrench = document.querySelector(
        'i.grok-icon.fa-wrench[aria-label="Peptides analysis settings"]') as HTMLElement | null;
      if (wrench) wrench.click();
    });
    await page.locator('[name="dialog-Peptides-settings"]').waitFor({timeout: 8000});

    const ok = await page.evaluate(async () => {
      const dlg = document.querySelector('[name="dialog-Peptides-settings"]');
      if (!dlg) return {error: 'dialog not found'};
      // Expand Viewers pane only if its content is actually collapsed.
      const viewersPane = Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
        .find((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim() === 'Viewers');
      if (viewersPane) {
        const content = viewersPane.querySelector('.d4-accordion-pane-content') as HTMLElement | null;
        const isCollapsed = !!content && getComputedStyle(content).display === 'none';
        if (isCollapsed) {
          const h = viewersPane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
          if (h) h.click();
          await new Promise((r) => setTimeout(r, 600));
        }
      }
      const cb = dlg.querySelector('[name="input-Active-peptide-selection"]') as HTMLInputElement | null;
      if (!cb) return {error: 'cb not found'};
      const initialChecked = cb.checked;
      // Entry state is FALSE (Scenario 1 leaves showClusterMaxActivity undefined); click twice to fire OFF.
      if (initialChecked === false) {
        cb.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, composed: true, view: window}));
        await new Promise((r) => setTimeout(r, 300));
      }
      cb.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, composed: true, view: window}));
      await new Promise((r) => setTimeout(r, 300));
      const cbStateBeforeOk = cb.checked;
      const okBtn = dlg.querySelector('[name="button-OK"]') as HTMLElement | null;
      if (!okBtn) return {error: 'OK btn not found'};
      okBtn.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, composed: true, view: window}));
      return {initialChecked, cbStateBeforeOk};
    });
    expect(ok && (ok as any).error, '[name="input-Active-peptide-selection"] driving setup failed').toBeFalsy();
    expect((ok as any).cbStateBeforeOk,
      'Active-peptide-selection checkbox should be OFF (false) before OK').toBe(false);

    await page.waitForFunction(() =>
      !document.querySelector('[name="dialog-Peptides-settings"]'), null, {timeout: 8000});

    let cmaAttachedSettled = true;
    for (let i = 0; i < 12; i++) {
      cmaAttachedSettled = await page.evaluate((VT) => {
        const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
        const model = tv.dataFrame.temp['peptidesModel'];
        return !!model.findViewer(VT.CLUSTER_MAX_ACTIVITY);
      }, VIEWER_TYPE);
      if (!cmaAttachedSettled) break;
      await page.waitForTimeout(500);
    }

    const after = await page.evaluate((VT) => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = tv.dataFrame.temp['peptidesModel'];
      return {
        cmaAttached: !!model.findViewer(VT.CLUSTER_MAX_ACTIVITY),
        viewers: Array.from(tv.viewers).map((v) => v.type),
        findSvm: !!model.findViewer(VT.SEQUENCE_VARIABILITY_MAP),
        findMpr: !!model.findViewer(VT.MOST_POTENT_RESIDUES),
        showCMASetting: model.settings?.showClusterMaxActivity ?? null,
      };
    }, VIEWER_TYPE);
    expect(after.showCMASetting,
      'settings.showClusterMaxActivity should be false after the click-twice OFF round-trip').toBe(false);
    expect(after.cmaAttached,
      'findViewer(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY) should be null after the closeViewer dispatch').toBe(false);
    expect(after.findSvm, 'SVM viewer must persist after CMA toggle-off').toBe(true);
    expect(after.findMpr, 'MPR viewer must persist after CMA toggle-off').toBe(true);
  });

  // Re-add path: toggle CMA ON from a non-empty analysisView (single click; entry state is FALSE).
  await softStep('Scenario 2 (steps 5-7): toggle Cluster max activity ON again, verify re-add path', async () => {
    await page.evaluate(async () => {
      const wrench = document.querySelector(
        'i.grok-icon.fa-wrench[aria-label="Peptides analysis settings"]') as HTMLElement | null;
      if (wrench) wrench.click();
    });
    await page.locator('[name="dialog-Peptides-settings"]').waitFor({timeout: 8000});

    const ok = await page.evaluate(async () => {
      const dlg = document.querySelector('[name="dialog-Peptides-settings"]');
      if (!dlg) return {error: 'dialog not found'};
      // Expand Viewers pane only if its content is actually collapsed.
      const viewersPane = Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
        .find((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim() === 'Viewers');
      if (viewersPane) {
        const content = viewersPane.querySelector('.d4-accordion-pane-content') as HTMLElement | null;
        const isCollapsed = !!content && getComputedStyle(content).display === 'none';
        if (isCollapsed) {
          const h = viewersPane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
          if (h) h.click();
          await new Promise((r) => setTimeout(r, 600));
        }
      }
      const cb = dlg.querySelector('[name="input-Active-peptide-selection"]') as HTMLInputElement | null;
      if (!cb) return {error: 'cb not found'};
      const initialChecked = cb.checked;
      if (initialChecked === false) {
        cb.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, composed: true, view: window}));
        await new Promise((r) => setTimeout(r, 300));
      }
      const cbStateBeforeOk = cb.checked;
      const okBtn = dlg.querySelector('[name="button-OK"]') as HTMLElement | null;
      if (!okBtn) return {error: 'OK btn not found'};
      okBtn.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, composed: true, view: window}));
      return {initialChecked, cbStateBeforeOk};
    });
    expect(ok && (ok as any).error, 'Active-peptide-selection driving setup failed').toBeFalsy();
    expect((ok as any).cbStateBeforeOk, 'Active-peptide-selection checkbox should be ON before OK').toBe(true);

    await page.waitForFunction(() =>
      !document.querySelector('[name="dialog-Peptides-settings"]'), null, {timeout: 8000});

    // Poll for re-add: addClusterMaxActivityViewer is async (awaits dockManager.dock).
    let cmaAttachedSettled = false;
    for (let i = 0; i < 16; i++) {
      cmaAttachedSettled = await page.evaluate((VT) => {
        const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
        const model = tv.dataFrame.temp['peptidesModel'];
        return !!model.findViewer(VT.CLUSTER_MAX_ACTIVITY);
      }, VIEWER_TYPE);
      if (cmaAttachedSettled) break;
      await page.waitForTimeout(500);
    }

    const after = await page.evaluate((VT) => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = tv.dataFrame.temp['peptidesModel'];
      return {
        cmaAttached: !!model.findViewer(VT.CLUSTER_MAX_ACTIVITY),
        viewers: Array.from(tv.viewers).map((v) => v.type),
        findSvm: !!model.findViewer(VT.SEQUENCE_VARIABILITY_MAP),
        findMpr: !!model.findViewer(VT.MOST_POTENT_RESIDUES),
        findMcl: !!model.findViewer(VT.MCL),
        showCMASetting: model.settings?.showClusterMaxActivity ?? null,
      };
    }, VIEWER_TYPE);
    expect(after.showCMASetting,
      'settings.showClusterMaxActivity should be truthy after the click-once ON round-trip').toBeTruthy();
    expect(after.cmaAttached,
      'findViewer(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY) should be non-null after the addClusterMaxActivityViewer re-add').toBe(true);
    expect(after.findSvm, 'SVM viewer must remain attached across the round-trip').toBe(true);
    expect(after.findMpr, 'MPR viewer must remain attached across the round-trip').toBe(true);
    expect(after.findMcl, 'MCL viewer must remain attached across the round-trip').toBe(true);
  });

  await softStep('Scenario 2 (step 8): no null-receiver crash across the toggle-off + toggle-on round-trip', async () => {
    const lastError = await page.evaluate(() =>
      grok.shell.lastError ? String(grok.shell.lastError) : null);
    const fatal = lastError && /setTrue|fire.*on (null|undefined)|Cannot read .* (null|undefined)|method not found.*null/i.test(lastError);
    expect(fatal,
      `Scenario 2 step 8 invariant: toggle-off + toggle-on produced a null-receiver / fatal error: ${lastError}`)
      .toBeFalsy();
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
