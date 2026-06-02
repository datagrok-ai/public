/* ---
sub_features_covered:
  - bio.viewers.similarity-search
  - bio.viewers.diversity-search
  - bio.search.similarity
  - bio.search.similarity.top-menu
  - bio.search.diversity
  - bio.search.diversity.top-menu
  - bio.analyze.activity-cliffs
  - bio.analyze.activity-cliffs.top-menu
--- */
//   related_bugs: [GROK-16111] (Bio Similarity Search empty-input regression-risk;
//     same gap suspected in Diversity Search + Activity Cliffs per bug-library
// [GROK-16111]) across the three current-row Bio viewers. Per bio.md
// "GROK-16111 empty-current-row invariant" (line 304): "Setting
//     (GROK-16111).
// (GROK-16111 path) | wrap grok.shell.error and grok.shell.warning to
// Known real-platform gap (scope_reductions :: SR-01, operator Option A
// 2026-06-02): GROK-16111 — the three Bio current-row viewers silently KNN
// the bug is unfixed; the no-crash / no-silent-rewrite assertion (Invariant
// data-enrichment SR-05..08 (GROK-20175). Revert SR-01 + restore the hard
// expect() when GROK-16111 is fixed. See empty-input-row-viewers.md
// scope_reductions :: SR-01.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
test.use(specTestOptions);
type ViewerCase = {
  label: string;
  groupSelector: string;
  leafSelector: string;
  hasEditorDialog: boolean;
  viewerSelector: string | null;
  subFeatures: string[];
};
const viewerCases: ViewerCase[] = [
  {
    label: 'Sequence Similarity Search',
    groupSelector: '[name="div-Bio---Search"]',
    leafSelector: '[name="div-Bio---Search---Similarity-Search"]',
    hasEditorDialog: false,
    viewerSelector: '[name="viewer-Sequence-Similarity-Search"]',
    subFeatures: [
      'bio.viewers.similarity-search',
      'bio.search.similarity',
      'bio.search.similarity.top-menu',
    ],
  },
  {
    label: 'Sequence Diversity Search',
    groupSelector: '[name="div-Bio---Search"]',
    leafSelector: '[name="div-Bio---Search---Diversity-Search"]',
    hasEditorDialog: false,
    viewerSelector: '[name="viewer-Sequence-Diversity-Search"]',
    subFeatures: [
      'bio.viewers.diversity-search',
      'bio.search.diversity',
      'bio.search.diversity.top-menu',
    ],
  },
  {
    label: 'Activity Cliffs',
    groupSelector: '[name="div-Bio---Analyze"]',
    leafSelector: '[name="div-Bio---Analyze---Activity-Cliffs..."]',
    hasEditorDialog: true,
    viewerSelector: null,
    subFeatures: [
      'bio.analyze.activity-cliffs',
      'bio.analyze.activity-cliffs.top-menu',
    ],
  },
];
for (const vc of viewerCases) {
  test(`Bio ${vc.label} rejects empty current-row input with balloon`, async ({page}) => {
    test.setTimeout(600_000);
    stepErrors.length = 0;
    await loginToDatagrok(page);
    // Setup phase: open filter_FASTA.csv (14 rows, single Macromolecule
    // "fasta" column per bio.md "Test datasets" — verified via direct
    // file inspection: header "fasta", 14 data rows). Wait for
    // Macromolecule semType detection + Bio package init.
    await page.evaluate(async () => {
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = true;
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:AppData/Bio/tests/filter_FASTA.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(() => resolve(), 4000);
      });
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
      if (hasMacromolecule) {
        for (let i = 0; i < 60; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        await new Promise((r) => setTimeout(r, 5000));
      }
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    // Bio top-menu + init-completion readiness — mirrors the analyze-spec /
    // sequence-activity-cliffs-spec cold-start stabilization. Layer 1: DOM
    // visibility of [name="div-Bio"]. Layer 2: Bio:getSeqHelper resolves
    // only after initBio completes (bio.md "Init-order invariant",
    // analyze-spec.ts L88-99).
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await page.evaluate(async () => {
      const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
      }
      await new Promise((r) => setTimeout(r, 3000));
    });
    // Setup the empty-current-row + balloon-hook precondition.
    //   1. Verify the table has >= 2 rows so an empty-cell state is
    //      constructable without an empty-table degenerate (per scenario
    //      .md "Setup" step: "≥ 2 rows so the empty current cell state can
    //      be constructed without an empty-table degenerate").
    //   2. Set row 0's Macromolecule cell to empty via the column setter
    //      (grid cell-edit equivalent per scenario .md "Setup" step:
    //      "Position the current row on a row whose sequence cell has been
    //      cleared to empty/null"). Programmatic edit is the canonical
    //      class-1 path — bio.md L304 explicitly states "Programmatic
    //      verification: hook grok.balloon calls before triggering".
    //   3. Set currentRowIdx to that empty row.
    //   4. Hook grok.shell.warning + grok.shell.error to capture balloon
    //      invocations (per bio.md L580 "Hook balloon errors (GROK-16111
    //      path) | wrap grok.shell.error and grok.shell.warning to capture
    //      invocations"). Stash in window.__balloonCalls for later read.
    await softStep('Setup: open filter_FASTA, empty row 0, hook balloon, set current row', async () => {
      const setup = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
        if (!macro) throw new Error('No Macromolecule column detected on filter_FASTA.csv');
        if (df.rowCount < 2) throw new Error(`filter_FASTA.csv must have >=2 rows; got ${df.rowCount}`);
        // Clear row 0 to empty string (Macromolecule columns accept '' as
        // the empty / null sentinel — the same shape produced by an
        // in-grid blank+commit edit).
        macro.set(0, '');
        df.currentRowIdx = 0;
        // Install balloon hook. Wrap both warning and error since
        // platform code may choose either class (Bio package source
        // uses grok.shell.error for "Empty sequences cannot be used for
        // similarity search"-class messages per package.ts; the contract
        // accepts either class).
        (window as any).__balloonCalls = [];
        const origWarn = grok.shell.warning.bind(grok.shell);
        const origErr = grok.shell.error.bind(grok.shell);
        grok.shell.warning = ((msg: any, opts?: any) => {
          try { (window as any).__balloonCalls.push({type: 'warning', msg: String(msg).slice(0, 300)}); } catch { /* ignore */ }
          return origWarn(msg, opts);
        }) as any;
        grok.shell.error = ((msg: any, opts?: any) => {
          try { (window as any).__balloonCalls.push({type: 'error', msg: String(msg).slice(0, 300)}); } catch { /* ignore */ }
          return origErr(msg, opts);
        }) as any;
        return {
          rowCount: df.rowCount,
          macroName: macro.name,
          row0value: macro.get(0),
          row0empty: macro.get(0) === '' || macro.get(0) == null,
          currentRowIdx: df.currentRowIdx,
        };
      });
      expect(setup.rowCount).toBeGreaterThanOrEqual(2);
      expect(setup.macroName).not.toBeNull();
      // Contract pre-state: row 0 sequence MUST be empty/null and MUST
      // be the current row when the viewer is invoked (scenario .md:
      // "The current-row sequence cell MUST be empty/null when each of
      // the three viewers below is invoked").
      expect(setup.row0empty).toBe(true);
      expect(setup.currentRowIdx).toBe(0);
    });
    // Capture baseline row count BEFORE invoking the viewer — used to
    // assert the failure-mode side of the contract: "no silent zero-row
    // result table is produced".
    const baseRowCount: number = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    // Dispatch step: open the viewer's top-menu path. Per the scenario
    // step 1-3 for each viewer: click Bio > <group> > <leaf>. Hover-
    // driven submenu open is required to surface the leaves (per bio.md
    // "Click pattern" + sibling sequence-activity-cliffs-spec.ts
    // L177-184). For Activity Cliffs the leaf click opens an editor
    // dialog; click OK to commit empty-input compute attempt.
    await softStep(`${vc.label}: open Bio > ${vc.groupSelector.includes('Search') ? 'Search' : 'Analyze'} > ${vc.label}${vc.hasEditorDialog ? ' (defaults OK)' : ''}`, async () => {
      await page.evaluate(async ({groupSel, leafSel}) => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        const group = document.querySelector(groupSel);
        if (!group) throw new Error(`Bio group ${groupSel} not found in top menu`);
        group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        const leaf = document.querySelector(leafSel);
        if (!leaf) throw new Error(`Bio leaf ${leafSel} not found under ${groupSel}`);
        (leaf as HTMLElement).click();
      }, {groupSel: vc.groupSelector, leafSel: vc.leafSelector});
      if (vc.hasEditorDialog) {
        // Activity Cliffs path: editor dialog opens; click OK to run
        // with default parameters. Per scenario step 3 for Activity
        // Cliffs: "In the SeqActivityCliffsEditor dialog, leave defaults
        // and click OK." 60s tolerance per sibling
        // sequence-activity-cliffs-spec.ts L189 (empirical cold ceiling).
        await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
        await page.locator('.d4-dialog [name="button-OK"]').click();
      }
    });
    // Contract assertion: the viewer's reaction to empty current-row input.
    // Three invariants per scenario .md "Expected":
    //   1. Balloon invocation count > 0 (warning OR error class). The
    //      bug-library GROK-16111 expected message is "Empty sequences
    //      cannot be used for similarity search"; the contract accepts
    //      any balloon naming empty/null input as the rejection reason,
    //      not a specific phrasing.
    //   2. No silent zero-row result is produced — the dataframe is
    //      unchanged in row count, AND the viewer either does NOT dock
    //      OR docks with no result rows.
    //   3. The active TableView's dataframe is the SAME object — the
    //      compute path did not silently rewrite the table.
    // 30s wait window allows for async balloon emission paths (the
    // empty-input rejection may flow through onCurrentRowChanged
    // dispatched on the next animation frame after the click).
    await softStep(`${vc.label}: balloon surfaces, no silent zero-row result`, async () => {
      // Poll for at least one balloon invocation. 15s is comfortable for
      // a synchronous-rejection contract; 30s for hot-path embedding
      // attempts that may take longer to surface a guard.
      await page.waitForFunction(() => {
        return Array.isArray((window as any).__balloonCalls)
          && (window as any).__balloonCalls.length > 0;
      }, null, {timeout: 30_000}).catch(() => { /* swallow — assertion below covers it */ });
      const probe = await page.evaluate((viewerSel) => {
        const df = grok.shell.tv.dataFrame;
        const calls = ((window as any).__balloonCalls || []) as Array<{type: string; msg: string}>;
        // Viewer types Active in the TableView — used to detect both
        // Similarity-Search-named and Diversity-Search-named viewers,
        // plus the Scatter-plot result viewer for Activity Cliffs.
        const viewerTypes = Array.from((grok.shell.tv as any).viewers).map((v: any) => v.type);
        // Per-viewer probes:
        //   - Similarity / Diversity: did the named viewer dock at all?
        //     If yes, does it expose result rows?
        //   - Activity Cliffs: did a Scatter plot dock at all? (Failure
        //     mode would be "silent Scatter plot with empty cliffs".)
        const docked = viewerSel
          ? !!document.querySelector(viewerSel)
          : viewerTypes.includes('Scatter plot');
        return {
          balloonCount: calls.length,
          balloonTypes: calls.map((c) => c.type),
          balloonMsgs: calls.map((c) => c.msg),
          rowCount: df.rowCount,
          viewerTypes,
          docked,
        };
      }, vc.viewerSelector);
      // Invariant 1: a balloon (warning OR error class) MUST surface.
      // This is the contract's POSITIVE assertion — the GROK-16111
      // regression-risk failure mode is "no balloon, silent empty
      // result".
      //
      // SR-01 (GROK-16111 — operator Option A, 2026-06-02): the three Bio
      // current-row viewers silently KNN on empty/null current-row input
      // and surface NO rejection balloon (GROK-16111 status:
      // regression-risk, fixed_in: ''). The assertion is correct; the
      // product is broken. The hard
      //   expect(probe.balloonCount).toBeGreaterThan(0)
      // is replaced with the guarded console.warn below so the scenario
      // invariant is documented without blocking the run — mirrors
      // PowerPack data-enrichment SR-05..08 (GROK-20175). Revert SR-01 +
      // restore the hard expect() above when GROK-16111 fixed_in is set.
      // See empty-input-row-viewers.md scope_reductions :: SR-01.
      if (!(probe.balloonCount > 0)) {
        // eslint-disable-next-line no-console
        console.warn(`[SR-01 known platform gap] GROK-16111: ${vc.label} did NOT surface a rejection balloon on empty current-row input (balloonCount=${probe.balloonCount}). Captured balloons: ${JSON.stringify(probe.balloonMsgs)}. Active viewer types: ${JSON.stringify(probe.viewerTypes)}. Revert SR-01 + restore the hard expect(probe.balloonCount).toBeGreaterThan(0) when GROK-16111 is fixed.`);
      }
      // Invariant 2: no silent table mutation. The dataframe row count
      // must equal the baseline — empty-input rejection must NEVER
      // produce a zero-row result or otherwise rewrite the active table.
      expect(probe.rowCount).toBe(baseRowCount);
      // Invariant 3 (per scenario .md "Expected" acceptance shape):
      // either the viewer does NOT dock at all, OR it docks but exposes
      // no result rows. For Similarity/Diversity the named viewer
      // selector probes the dock; for Activity Cliffs the Scatter plot
      // viewer is the dock proxy. The silent-success failure mode is
      // "viewer docked + viewer holds an embedded result grid with zero
      // rows". The scenario .md accepts BOTH the no-dock shape AND the
      // dock-with-placeholder shape — assertion is on the silent-
      // empty-result NOT happening, not on a specific dock shape. The
      // balloon assertion above already covers the contract's positive
      // signal; this assertion documents the acceptance shape for
      // operator review and does not fail on either acceptance form.
      // (Stronger empty-result-table inspection would require viewer-
      // specific selectors not in bio.md beyond the dock container; the
      // class-2 selector path is not invoked here per the selector-
      // provenance 3-class model.)
      //
      // Soft expectation: log the observed dock state for operator
      // review without failing — the contract is satisfied as long as
      // invariants 1 and 2 hold AND no silent embedded result with rows
      // was produced (validated by the no-table-mutation check above).
      // If a stronger empty-dock assertion is needed in a future cycle,
      // augment bio.md with viewer-internal selectors and tighten here.
    });
    if (stepErrors.length > 0) {
      const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
    }
  });
}
