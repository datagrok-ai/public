import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const datasetPath = 'System:DemoFiles/bio/peptides.csv';
function isNullReceiverCrash(lastError: string): boolean {
  return /setTrue|method not found.*null|Cannot read propert.*null|reading 'getRawData'/.test(lastError);
}
test('Collaborative selection — WebLogo header click propagates through fireBitsetChanged', async ({page}) => {
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
      return {rows: df.rowCount, semType: df.col('AlignedSequence')?.semType ?? null};
    }, datasetPath);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
    expect(result.rows, 'peptides.csv should load 647 rows').toBe(647);
    expect(result.semType, 'AlignedSequence must be a Macromolecule column').toBe('Macromolecule');
  });
  await softStep('Setup step 2: launch SAR from the top menu (default config)', async () => {
    const bio = page.locator('[name="div-Bio"]');
    await bio.waitFor({state: 'attached', timeout: 30_000});
    const bioVisible = await bio.evaluate((el) => (el as HTMLElement).offsetParent !== null);
    expect(bioVisible,
      '[name="div-Bio"] present but not visible (offsetParent null — menu in overflow at this viewport)')
      .toBe(true);
    await bio.click();
    await page.waitForTimeout(700);
    const opened = await page.evaluate(async () => {
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
      let okClicked = false;
      if (dlg) {
        const ok = (dlg.querySelector('[name="button-OK"]')
          ?? document.querySelector('[name="button-OK"]')) as HTMLElement | null;
        if (ok) { ok.click(); okClicked = true; }
      }
      return {analyzeFound: !!analyze, sarFound: !!sar, dialogFound: !!dlg, okClicked};
    });
    expect(opened.analyzeFound, '[name="div-Bio---Analyze"] submenu entry not found').toBe(true);
    expect(opened.sarFound, '[name="div-Bio---Analyze---SAR..."] menu item not found').toBe(true);
    expect(opened.dialogFound, '[name="dialog-Analyze-Peptides"] config dialog did not open').toBe(true);
    expect(opened.okClicked, 'Analyze Peptides dialog OK was not clicked').toBe(true);
    await page.waitForFunction(() => {
      return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
    }, {timeout: 90000});
    await page.waitForTimeout(8000);
  });
  await softStep('Setup step 3: confirm SAR layout + empty pre-exercise selection', async () => {
    await page.evaluate(() => {
      (window as any).__paneHasContent = async (paneName: string) => {
        const deadline = Date.now() + 8000;
        let last: any = {found: false, hasContent: false};
        while (Date.now() < deadline) {
          const pane = document.querySelector(`[name="${paneName}"]`);
          if (pane) {
            const header = pane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
            if (header && !pane.classList.contains('expanded')) header.click();
            const txt = pane.textContent || '';
            const hasContent =
              !!pane.querySelector('table') ||
              !!pane.querySelector('.d4-grid') ||
              !!pane.querySelector('.d4-histogram') ||
              !!pane.querySelector('.d4-viewer') ||
              /Count|Mean activity|Mean difference/.test(txt);
            last = {found: true, hasContent, childCount: pane.querySelectorAll('*').length};
            if (hasContent) return last;
          }
          await new Promise((r) => setTimeout(r, 300));
        }
        return last;
      };
    });
    const state = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const model = tv.dataFrame.temp['peptidesModel'] as any;
      df.selection.setAll(false);
      const viewers = Array.from(tv.viewers).map((v) => v.type);
      const monomerCols: string[] = [];
      for (let i = 0; i < df.columns.length; i++) {
        const c = df.columns.byIndex(i);
        if (c.semType === 'Monomer') monomerCols.push(c.name);
      }
      return {
        viewers,
        colHeaderHeight: tv.grid?.props?.colHeaderHeight ?? 0,
        monomerColCount: monomerCols.length,
        selBefore: df.selection.trueCount,
        modelPresent: !!model,
      };
    });
    expect(state.modelPresent, 'PeptidesModel did not attach after SAR launch').toBe(true);
    expect(state.viewers, 'Sequence Variability Map (MonomerPosition) must attach').toContain('Sequence Variability Map');
    expect(state.viewers, 'Most Potent Residues must attach').toContain('Most Potent Residues');
    expect(state.viewers, 'MCL clustering viewer must attach').toContain('MCL');
    expect(state.colHeaderHeight,
      'grid colHeaderHeight did not grow for the WebLogo column-headers').toBeGreaterThan(40);
    expect(state.monomerColCount, 'no per-position Monomer columns were split out by SAR').toBeGreaterThan(0);
    expect(state.selBefore, 'selection should be empty before exercising the backbone').toBe(0);
  });
  await softStep('Scenario 1 (steps 1-8): WebLogo header pick broadcasts to all surfaces', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const model = tv.dataFrame.temp['peptidesModel'] as any;
      const mps = model.monomerPositionStats;
      const fresh: Record<string, string[]> = {};
      for (const pos of Object.keys(mps).filter((k) => k !== 'general')) fresh[pos] = [];
      model.webLogoSelection = fresh;
      df.selection.setAll(false);
      let pick: {position: string; monomer: string; count: number} | null = null;
      for (const pos of Object.keys(mps).filter((k) => k !== 'general')) {
        const monomers = Object.keys(mps[pos]).filter((m) => m !== 'general');
        for (const mon of monomers) {
          const cnt = mps[pos][mon]?.count ?? 0;
          if (cnt > 0 && cnt < df.rowCount) { pick = {position: pos, monomer: mon, count: cnt}; break; }
        }
        if (pick) break;
      }
      if (!pick) return {pickFound: false};
      let threw: string | null = null;
      try {
        model.modifyWebLogoSelection(
          {positionOrClusterType: pick.position, monomerOrCluster: pick.monomer},
          {shiftPressed: false, ctrlPressed: false}, true);
        model.fireBitsetChanged('WebLogo');
      } catch (e) { threw = String(e); }
      await new Promise((r) => setTimeout(r, 1200));
      const selAfter = df.selection.trueCount;
      const combined = model.getCombinedSelection();
      const combinedCount = combined?.trueCount ?? null;
      const mapForPos = Array.isArray(model.webLogoSelection[pick.position])
        ? model.webLogoSelection[pick.position].slice() : null;
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const mpr = document.querySelector('[name="viewer-Most-Potent-Residues"]');
      const svmHasCanvas = svm ? !!svm.querySelector('canvas') : false;
      const mprPresent = !!mpr;
      const widgets: Record<string, any> = {};
      for (const paneName of ['pane-Distribution', 'pane-Selection']) {
        widgets[paneName] = await (window as any).__paneHasContent(paneName);
      }
      const lastError = grok.shell.lastError ? String(grok.shell.lastError) : '';
      return {
        pickFound: true, pick, threw, selAfter, combinedCount, mapForPos,
        svmHasCanvas, mprPresent, widgets, lastError,
      };
    });
    expect(result.pickFound, 'no populated (monomer, position) pick was found in the WebLogo stats').toBe(true);
    expect(result.threw, `WebLogo selection handler / fireBitsetChanged threw: ${result.threw}`).toBeNull();
    // Step 3: DataFrame.selection populated (exactly the picked-monomer rows).
    expect(result.selAfter, 'WebLogo header pick did not populate df.selection (backbone did not fire)')
      .toBeGreaterThan(0);
    expect(result.selAfter, 'df.selection trueCount should equal the picked monomer count')
      .toBe(result.pick!.count);
    // get-selection-bitset projection is consistent with the broadcast BitSet.
    expect(result.combinedCount, 'getCombinedSelection (get-selection-bitset) did not match df.selection')
      .toBe(result.selAfter);
    // The unified Selection map carries exactly the picked monomer at the picked position.
    expect(result.mapForPos, 'webLogoSelection did not record the picked monomer at the picked position')
      .toContain(result.pick!.monomer);
    // Steps 5-6: cross-viewer surfaces survive the broadcast.
    expect(result.svmHasCanvas, 'Sequence Variability Map lost its canvas after the broadcast').toBe(true);
    expect(result.mprPresent, 'Most Potent Residues viewer disappeared after the broadcast').toBe(true);
    // Steps 7-8: the Distribution + Selection widget surfaces carry rendered content.
    expect(result.widgets['pane-Distribution']?.found, 'Distribution widget pane not found').toBe(true);
    expect(result.widgets['pane-Distribution']?.hasContent,
      'Distribution widget pane has no rendered content for the selected subset').toBe(true);
    expect(result.widgets['pane-Selection']?.found, 'Selection widget pane not found').toBe(true);
    expect(result.widgets['pane-Selection']?.hasContent,
      'Selection widget pane has no rendered content for the selected subset').toBe(true);
    // GROK-14298 invariant: no null-receiver crash in the cross-viewer broadcast path.
    expect(isNullReceiverCrash(result.lastError!),
      `GROK-14298 invariant: the first broadcast produced a null-receiver crash: ${result.lastError}`)
      .toBe(false);
  });
  // Scenario 2 — Shift-click additive (union) then Ctrl-click toggle-off, asserting cross-surface
  // consistency and no null-receiver crash on either re-broadcast (GROK-14298 listener-mutation mode).
  await softStep('Scenario 2 (steps 1-8): Shift additive then Ctrl toggle-off re-broadcast', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const model = tv.dataFrame.temp['peptidesModel'] as any;
      const mps = model.monomerPositionStats;
      // Re-establish the Scenario-1 single pick deterministically (first partial-count item).
      const fresh: Record<string, string[]> = {};
      for (const pos of Object.keys(mps).filter((k) => k !== 'general')) fresh[pos] = [];
      model.webLogoSelection = fresh;
      df.selection.setAll(false);
      const partials: {position: string; monomer: string; count: number}[] = [];
      for (const pos of Object.keys(mps).filter((k) => k !== 'general')) {
        const monomers = Object.keys(mps[pos]).filter((m) => m !== 'general');
        for (const mon of monomers) {
          const cnt = mps[pos][mon]?.count ?? 0;
          if (cnt > 0 && cnt < df.rowCount) { partials.push({position: pos, monomer: mon, count: cnt}); break; }
        }
        if (partials.length >= 2) break;
      }
      if (partials.length < 2) return {twoPicksFound: false};
      const [first, second] = partials;
      // Step 1 (re-establish Scenario-1 single pick).
      model.modifyWebLogoSelection(
        {positionOrClusterType: first.position, monomerOrCluster: first.monomer},
        {shiftPressed: false, ctrlPressed: false}, true);
      model.fireBitsetChanged('WebLogo');
      await new Promise((r) => setTimeout(r, 800));
      const selSingle = df.selection.trueCount;
      // Step 2: Shift-click a second distinct (monomer, position) — additive (union).
      let shiftThrew: string | null = null;
      try {
        model.modifyWebLogoSelection(
          {positionOrClusterType: second.position, monomerOrCluster: second.monomer},
          {shiftPressed: true, ctrlPressed: false}, true);
        model.fireBitsetChanged('WebLogo');
      } catch (e) { shiftThrew = String(e); }
      await new Promise((r) => setTimeout(r, 1000));
      const selAdditive = df.selection.trueCount;
      const mapAfterShift = {
        first: (model.webLogoSelection[first.position] ?? []).slice(),
        second: (model.webLogoSelection[second.position] ?? []).slice(),
      };
      // Step 4: SVM reflects the additive selection (viewer survives the re-broadcast).
      const svmAfterShift = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const svmAliveAfterShift = svmAfterShift ? !!svmAfterShift.querySelector('canvas') : false;
      // Steps 5-6: widget surfaces re-render against the new selection (cold-stable predicate).
      const widgetsAfterShift: Record<string, boolean> = {};
      for (const paneName of ['pane-Distribution', 'pane-Selection']) {
        widgetsAfterShift[paneName] = (await (window as any).__paneHasContent(paneName)).hasContent;
      }
      const lastErrorAfterShift = grok.shell.lastError ? String(grok.shell.lastError) : '';
      // Step 7: Ctrl-click the second pick — toggle it off; re-broadcast the reduced selection.
      let ctrlThrew: string | null = null;
      try {
        model.modifyWebLogoSelection(
          {positionOrClusterType: second.position, monomerOrCluster: second.monomer},
          {shiftPressed: false, ctrlPressed: true}, true);
        model.fireBitsetChanged('WebLogo');
      } catch (e) { ctrlThrew = String(e); }
      await new Promise((r) => setTimeout(r, 1000));
      const selReduced = df.selection.trueCount;
      const mapAfterCtrl = {
        first: (model.webLogoSelection[first.position] ?? []).slice(),
        second: (model.webLogoSelection[second.position] ?? []).slice(),
      };
      const combinedReduced = model.getCombinedSelection()?.trueCount ?? null;
      // Step 8: cross-surfaces reflect the reduced selection.
      const svmAfterCtrl = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const svmAliveAfterCtrl = svmAfterCtrl ? !!svmAfterCtrl.querySelector('canvas') : false;
      const widgetsAfterCtrl: Record<string, boolean> = {};
      for (const paneName of ['pane-Distribution', 'pane-Selection']) {
        widgetsAfterCtrl[paneName] = (await (window as any).__paneHasContent(paneName)).hasContent;
      }
      const lastErrorAfterCtrl = grok.shell.lastError ? String(grok.shell.lastError) : '';
      return {
        twoPicksFound: true, first, second,
        selSingle, selAdditive, selReduced, combinedReduced,
        mapAfterShift, mapAfterCtrl,
        shiftThrew, ctrlThrew,
        svmAliveAfterShift, svmAliveAfterCtrl,
        widgetsAfterShift, widgetsAfterCtrl,
        lastErrorAfterShift, lastErrorAfterCtrl,
      };
    });
    expect(result.twoPicksFound, 'need two distinct partial-count WebLogo picks for the modifier flow').toBe(true);
    expect(result.shiftThrew, `Shift-add WebLogo handler threw: ${result.shiftThrew}`).toBeNull();
    expect(result.selAdditive, 'Shift-click did not grow the selection (additive union expected, not replacement)')
      .toBeGreaterThan(result.selSingle!);
    // The unified Selection map carries BOTH picks after the Shift-add.
    expect(result.mapAfterShift!.first, 'first pick lost from the Selection map after Shift-add')
      .toContain(result.first!.monomer);
    expect(result.mapAfterShift!.second, 'second pick not recorded in the Selection map after Shift-add')
      .toContain(result.second!.monomer);
    // Step 4: SVM survives the re-broadcast. Steps 5-6: widget surfaces re-render.
    expect(result.svmAliveAfterShift, 'Sequence Variability Map lost its canvas after the Shift re-broadcast').toBe(true);
    expect(result.widgetsAfterShift!['pane-Distribution'],
      'Distribution widget did not re-render after the Shift re-broadcast').toBe(true);
    expect(result.widgetsAfterShift!['pane-Selection'],
      'Selection widget did not re-render after the Shift re-broadcast').toBe(true);
    // No null-receiver crash on the Shift re-broadcast.
    expect(isNullReceiverCrash(result.lastErrorAfterShift!),
      `GROK-14298 invariant: the Shift re-broadcast produced a null-receiver crash: ${result.lastErrorAfterShift}`)
      .toBe(false);
    // Step 7: Ctrl-toggle-off did not throw; selection returns to the single-pick state.
    expect(result.ctrlThrew, `Ctrl-toggle WebLogo handler threw: ${result.ctrlThrew}`).toBeNull();
    expect(result.selReduced, 'Ctrl-toggle-off did not return the selection to the single-pick count')
      .toBe(result.selSingle);
    // get-selection-bitset projection matches the reduced selection.
    expect(result.combinedReduced, 'getCombinedSelection did not match the reduced df.selection')
      .toBe(result.selReduced);
    // The second pick is removed from the unified Selection map; the first remains.
    expect(result.mapAfterCtrl!.second, 'Ctrl-toggle did not remove the second pick from the Selection map')
      .not.toContain(result.second!.monomer);
    expect(result.mapAfterCtrl!.first, 'Ctrl-toggle erroneously dropped the first (untouched) pick')
      .toContain(result.first!.monomer);
    // Step 8: cross-surfaces reflect the reduced selection.
    expect(result.svmAliveAfterCtrl, 'Sequence Variability Map lost its canvas after the Ctrl re-broadcast').toBe(true);
    expect(result.widgetsAfterCtrl!['pane-Distribution'],
      'Distribution widget did not re-render after the Ctrl re-broadcast').toBe(true);
    expect(result.widgetsAfterCtrl!['pane-Selection'],
      'Selection widget did not re-render after the Ctrl re-broadcast').toBe(true);
    // No null-receiver crash on the Ctrl re-broadcast (the listener-mutation crash mode).
    expect(isNullReceiverCrash(result.lastErrorAfterCtrl!),
      `GROK-14298 invariant: the Ctrl re-broadcast produced a null-receiver crash: ${result.lastErrorAfterCtrl}`)
      .toBe(false);
  });
  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
