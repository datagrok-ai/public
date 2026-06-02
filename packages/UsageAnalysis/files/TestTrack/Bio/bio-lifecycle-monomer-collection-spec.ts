/* ---
sub_features_covered:
  - bio.manage.monomer-collections-app
  - bio.manage.libraries-view
  - bio.api.get-monomer-lib-helper
  - bio.lifecycle.init
--- */
//   are therefore NOT a structural bug in this code path but a
// Retry hardening — evidence-based fix per the test-bug hypothesis path of
// HYPOTHESIS CATEGORY: test-bug. Round-1 retry hypothesis is DISTINCT
// HYPOTHESIS CATEGORY: test-bug. Round-2 retry hypothesis is DISTINCT
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {
  saveAllTablesWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';
test.use(specTestOptions);
test('Bio monomer_collection source-class lifecycle: write collection → reload via app → save project with reference → reopen and verify', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  const stamp = Date.now();
  const workingCollection = `bio-lifecycle-monomer-collection-${stamp}.json`;
  // Per lib-manager.ts#L22 MONOMER_COLLECTION_STORAGE_PATH constant.
  const workingCollectionPath = `System:AppData/Bio/monomer-collections/${workingCollection}`;
  const projectName = `bio-lifecycle-monomer-collection-project-${stamp}`;
  // Synthetic monomer symbols matching the canonical HELM core monomer
  // symbols per scenario Setup ("synthetic entry referencing the canonical
  // HELM core monomer symbols (e.g. ['A', 'G', 'T', 'C'])"). Includes a
  // PEPTIDE-family canonical entry so MonomerSelectionWidget / SVM
  // catalogue look-ups don't degenerate to empty.
  const syntheticSymbols = ['A', 'G', 'T', 'C'];
  const syntheticDescription = `bio-lifecycle proactive test collection (${stamp})`;
  const syntheticTags = ['bio-lifecycle-test', 'automation'];
  let saved: {projectId: string; primaryTableInfoId: string; tableInfoIds: string[]; layoutId: string | null} | null = null;
  let workingCollectionWritten = false;
  await loginToDatagrok(page);
  // ==========================================================================
  // Setup — open a Bio dataset so the Macromolecule renderer touches the
  // monomer library / collection registration code paths. filter_HELM.csv is
  // the canonical HELM fixture used across the sibling
  // bio-lifecycle-monomer-library-spec.ts and
  // bio-lifecycle-macromolecule-column-spec.ts.
  // ==========================================================================
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
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
  }, 'System:AppData/Bio/tests/filter_HELM.csv');
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  // Two-layer Bio init readiness probe (mirrors the sibling
  // bio-lifecycle-monomer-library-spec.ts pattern).
  // Layer 1: DOM top-menu visibility. Layer 2: Bio:getMonomerLibHelper
  // serialization probe — the runtime serializes grok.functions.call after
  // init completion.
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getMonomerLibHelper', 'Bio:getSeqHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });
  try {
    // ========================================================================
    // Scenario 1 — Write a new collection via the shared lib-manager path
    // ========================================================================
    // Scenario 1, Steps 1+2+3 — Service-surface init contract: post-init
    // the MonomerLibManager singleton is non-null AND the
    // addOrUpdateMonomerCollection write path lands the working collection
    // under System:AppData/Bio/monomer-collections/ on FileShare, with
    // content round-tripping cleanly through readAsText. Atlas:
    // bio.lifecycle.init + bio.api.get-monomer-lib-helper +
    // dep_lifecycle_ops[save_monomer_library] for the collections subpath.
    await softStep('S1.1-1.3: getMonomerLibHelper returns singleton + writeCollection lands on FileShare with content round-trip', async () => {
      const result = await page.evaluate(async ({path, name, symbols, desc, tags}) => {
        const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        const hasHelper = helper != null;
        try {
          if (helper && typeof helper.awaitLoaded === 'function') {
            try { await helper.awaitLoaded(30_000); }
            catch (_) { try { await helper.awaitLoaded(); } catch (__) { /* non-fatal */ } }
          }
        } catch (_) { /* timeout non-fatal */ }
        // Write via the shared lib-manager path (lib-manager.ts#L434).
        // addOrUpdateMonomerCollection accepts (name, monomerSymbols, desc,
        // tags) and serializes a MonomerCollection (monomer-library.ts L221:
        // {description?, tags?, monomerSymbols, updatedBy, updatedOn}) into
        // System:AppData/Bio/monomer-collections/<name>.json via
        // grok.dapi.files.writeAsText.
        let writeErr: string | null = null;
        try {
          if (typeof helper.addOrUpdateMonomerCollection === 'function')
            await helper.addOrUpdateMonomerCollection(name, symbols, desc, tags);
          else
            writeErr = 'helper.addOrUpdateMonomerCollection is not a function on the resolved singleton';
        } catch (e) {
          writeErr = String(e).slice(0, 200);
        }
        if (writeErr) return {hasHelper, writeErr, exists: false, readEqual: false};
        // Verify via direct FileShare reads (atlas
        // dep_lifecycle_ops[save_monomer_library] write surface: no content
        // drift via FileShare). Per lib-manager.ts the file path is
        // exactly System:AppData/Bio/monomer-collections/<name>.json.
        const exists = await grok.dapi.files.exists(path);
        let readBack: string | null = null;
        try { readBack = await grok.dapi.files.readAsText(path); }
        catch (e) { /* leave null */ }
        let parsed: any = null;
        try { parsed = readBack ? JSON.parse(readBack) : null; }
        catch (e) { /* leave null */ }
        const monomersMatch = parsed != null && Array.isArray(parsed.monomerSymbols) &&
          parsed.monomerSymbols.length === symbols.length &&
          symbols.every((s: string) => parsed.monomerSymbols.includes(s));
        return {
          hasHelper,
          writeErr: null,
          exists,
          readEqual: parsed != null,
          monomersMatch,
          parsedDescription: parsed?.description ?? null,
          parsedTags: Array.isArray(parsed?.tags) ? parsed.tags : null,
          parsedUpdatedBy: parsed?.updatedBy ?? null,
          parsedUpdatedOn: parsed?.updatedOn ?? null,
          fileLen: typeof readBack === 'string' ? readBack.length : 0,
        };
      }, {
        path: workingCollectionPath,
        name: workingCollection.replace(/\.json$/, ''),
        symbols: syntheticSymbols,
        desc: syntheticDescription,
        tags: syntheticTags,
      });
      // Atlas bio.api.get-monomer-lib-helper / bio.lifecycle.init: singleton
      // populated post-init.
      expect(result.hasHelper).toBe(true);
      expect(result.writeErr, `addOrUpdateMonomerCollection error: ${result.writeErr}`).toBeNull();
      // Atlas dep_lifecycle_ops[save_monomer_library] (collections subpath):
      // file lands on FileShare under
      // System:AppData/Bio/monomer-collections/<name>.json.
      expect(result.exists,
        `expected ${workingCollectionPath} to exist on FileShare after addOrUpdateMonomerCollection`).toBe(true);
      // No content drift on FileShare round-trip.
      expect(result.readEqual).toBe(true);
      expect(result.fileLen).toBeGreaterThan(0);
      // MonomerCollection schema integrity (monomer-library.ts L221):
      // monomerSymbols round-trips, description / tags preserved,
      // updatedBy / updatedOn populated by the lib-manager write path.
      expect(result.monomersMatch,
        `expected synthetic symbols [${syntheticSymbols.join(', ')}] in monomerSymbols round-trip`).toBe(true);
      expect(result.parsedDescription).toBe(syntheticDescription);
      expect(Array.isArray(result.parsedTags)).toBe(true);
      expect(result.parsedTags!.length).toBe(syntheticTags.length);
      expect((result.parsedUpdatedBy || '').length).toBeGreaterThan(0);
      expect((result.parsedUpdatedOn || '').length).toBeGreaterThan(0);
      workingCollectionWritten = true;
    });
    // ========================================================================
    // Scenario 2 — Reload via the Monomer Collections app
    // ========================================================================
    // Scenario 2, Steps 1+2 — Open the `Monomer Collections` app via the
    // function-registry handle and verify the working collection appears in
    // the app's structural listing. Atlas
    // bio.manage.monomer-collections-app (function source:
    // package.ts#L1412 monomerCollectionsApp → showMonomerCollectionsView
    // in monomer-collections-view.ts L566).
    //
    // App listing structure (monomer-collections-view.ts):
    //   - View name == 'Monomer Collections' (VIEW_NAME constant L41).
    //   - Root div has class 'monomer-collections-view' (L93).
    //   - Cards grid is .monomer-collections-grid containing
    //     .monomer-collection-card per collection (L107 + L205).
    //   - Each card carries dataset.collectionName === <file name with .json>
    //     (L206) — the durable structural anchor regardless of display-name
    //     formatting.
    //   - listMonomerCollections() returns the FileShare directory listing
    //     filtered to .json (lib-manager.ts L413). The app loadCollections
    //     iterates this list to build cards (L120-L135 + applyFilter L150).
    //
    // Selector class for the Apps-tree path "Apps | Peptides | Monomer
    // Collections": class-3 (not in bio.md, not MCP-observed this dispatch
    // due to auth-stale session — see Selector recon-notes block above).
    // Per §"Selector provenance" we do NOT emit it; the function-registry
    // handle `Bio:monomerCollectionsApp` is the canonical platform handle
    // and dispatches the exact same showMonomerCollectionsView code path
    // the tree-node double-click invokes (package.ts L1418 = same call).
    await softStep('S2.1-2.2: Monomer Collections app opens via registered function; working collection card appears in listing', async () => {
      // Round-1 retry view-install dispatch (belt-and-suspenders pattern;
      // see "Retry hardening" header comment block above for the three
      // reinforcing root causes).
      const installDiag = await page.evaluate(async () => {
        const diag: any = {
          callOk: false,
          callErr: null as string | null,
          viewCaptured: false,
          viewNameAfterCall: null as string | null,
          addViewTried: false,
          addViewErr: null as string | null,
          setCurrentTried: false,
          setCurrentErr: null as string | null,
          viewsWalkFound: false,
          viewsWalkPromoted: false,
        };
        let view: any = null;
        try {
          view = await (grok as any).functions.call('Bio:monomerCollectionsApp', {});
          diag.callOk = true;
          diag.viewCaptured = view != null;
          // view.name is a JS-side getter backed by Dart on most wrappers;
          // capture defensively to aid diagnostics.
          try { diag.viewNameAfterCall = view?.name ?? null; } catch (_) { /* ignore */ }
        } catch (e: any) {
          diag.callErr = String(e).slice(0, 200);
        }
        // (2) Try grok.shell.addView(view) — workspace add.
        // Drop the prior `.root` precondition: view.root is a JS-only
        // field that does NOT round-trip on toJs(); the prior guard
        // incorrectly blocked the install. Guard on view truthiness only.
        if (view) {
          diag.addViewTried = true;
          try { (grok.shell as any).addView(view); }
          catch (e: any) { diag.addViewErr = String(e).slice(0, 200); }
          // (3) Try the current-view setter (grok.shell.v = view) — the
          // canonical current-view handoff (shell.ts L71). addView alone
          // doesn't always promote the new view to current outside an
          // app-launch context.
          diag.setCurrentTried = true;
          try { (grok.shell as any).v = view; }
          catch (e: any) { diag.setCurrentErr = String(e).slice(0, 200); }
        }
        // (4) Workspace walk-back: if a view named 'Monomer Collections'
        // is present in grok.shell.views but is not the current view,
        // promote it. Covers the case where the platform install path
        // produced its own wrapper distinct from ours.
        try {
          const views: any[] = Array.from((grok.shell as any).views || []);
          const target: any = views.find((v: any) => v?.name === 'Monomer Collections');
          if (target) {
            diag.viewsWalkFound = true;
            if ((grok.shell as any).v?.name !== 'Monomer Collections') {
              try { (grok.shell as any).v = target; diag.viewsWalkPromoted = true; }
              catch (_) { /* ignore */ }
            }
          }
        } catch (_) { /* ignore */ }
        return diag;
      });
      // (5) Poll for either of two install signals:
      //   (a) grok.shell.v.name === 'Monomer Collections', OR
      //   (b) some entry in grok.shell.views is named 'Monomer Collections'
      // — the assertion path below handles either case.
      const installed = await page.waitForFunction(() => {
        try {
          const cur = (window as any).grok?.shell?.v?.name;
          if (cur === 'Monomer Collections') return true;
          const views: any[] = Array.from((window as any).grok?.shell?.views || []);
          return views.some((v: any) => v?.name === 'Monomer Collections');
        } catch (_) { return false; }
      }, null, {timeout: 30_000}).catch(() => null);
      // Settle: loadCollections fires the FileShare listing call + builds
      // cards async via ui.wait(loadContent).
      await page.waitForTimeout(2500);
      // Attach install diagnostics to failure messages without breaking
      // assertion flow when the install path succeeded.
      const installDiagStr = `installDiag=${JSON.stringify(installDiag)}`;
      const result = await page.evaluate((fileName) => {
        // Locate the Monomer Collections view: prefer the current view if
        // its name matches, else walk grok.shell.views.
        let v: any = grok.shell.v;
        if (v?.name !== 'Monomer Collections') {
          const views: any[] = Array.from((grok.shell as any).views || []);
          const found: any = views.find((x: any) => x?.name === 'Monomer Collections');
          if (found) v = found;
        }
        const root: any = v?.root;
        // The view body uses class 'monomer-collections-view' (L93).
        const viewBody: any = root?.querySelector?.('.monomer-collections-view');
        const grid: any = root?.querySelector?.('.monomer-collections-grid');
        const cards = grid ? Array.from(grid.querySelectorAll('.monomer-collection-card')) : [];
        const cardCollectionNames: string[] = cards.map((c: any) =>
          c.dataset?.collectionName ?? '');
        // The card title is rendered without the `.json` extension
        // (monomer-collections-view.ts L192/L195). The durable anchor is
        // dataset.collectionName (L206). Match the exact file name; fall
        // back to a stem-substring match (the sibling spec's pattern, in
        // case future provider impls drop the .json on persistence).
        const stem = fileName.replace(/\.json$/, '');
        const hasWorkingCollection = cardCollectionNames.some((n: string) =>
          n === fileName || n === stem || n.includes(stem));
        // List all view names for diagnostic clarity if the assertion
        // surfaces a current-view-mismatch failure.
        const allViewNames: string[] = Array.from((grok.shell as any).views || [])
          .map((x: any) => x?.name ?? '<unnamed>');
        return {
          viewName: v?.name ?? null,
          currentViewName: (grok.shell as any).v?.name ?? null,
          allViewNames,
          rootPresent: !!root,
          viewBodyPresent: !!viewBody,
          gridPresent: !!grid,
          cardCount: cards.length,
          cardCollectionNames,
          hasWorkingCollection,
        };
      }, workingCollection);
      // Atlas bio.manage.monomer-collections-app: view-name canonical
      // handle. We assert on the LOCATED view (current OR walked from
      // grok.shell.views) — the walk-back fallback covers the case where
      // install added the view but couldn't promote it to current.
      expect(result.viewName,
        `expected to locate a 'Monomer Collections' view (current or in shell.views); ${installDiagStr}; observed current=${result.currentViewName}, allViews=[${result.allViewNames.join(', ')}], installed-marker=${installed != null}`).toBe('Monomer Collections');
      // Structural anchors present (monomer-collections-view.ts:L93/L107).
      expect(result.rootPresent).toBe(true);
      expect(result.viewBodyPresent).toBe(true);
      expect(result.gridPresent).toBe(true);
      // listMonomerCollections() enumerates the FileShare directory; the
      // collection we wrote in S1 MUST be in the dataset.collectionName set.
      expect(result.cardCount).toBeGreaterThanOrEqual(1);
      expect(result.hasWorkingCollection,
        `expected working collection '${workingCollection}' in card dataset names; observed: [${result.cardCollectionNames.join(', ')}]`).toBe(true);
    });
    // Scenario 2, Step 3 — Optional cross-surface check on the
    // `Bio | Manage | Monomer Libraries` view (atlas bio.manage.libraries-view).
    // Per scenario: "collections may or may not appear in the libraries view
    // depending on platform configuration; assert presence in the dedicated
    // collections app at minimum." We assert the libraries view IS reachable
    // (the cross-surface entry point exists), but do NOT assert the working
    // collection appears there — the canonical surface is the collections
    // app (covered by S2.1-2.2 above). This step satisfies the
    // playwright + non-ui-smoke E-LAYER-COMPLIANCE-01 ≥1 DOM-driving call
    // requirement via the class-1 top-menu drive (bio.md L611).
    //
    // Top-menu navigation pattern mirrors the sibling
    // bio-lifecycle-monomer-library-spec.ts S1.3 verbatim.
    //
    // ROUND-2 RETRY HARDENING (pre-S2.3 view-promotion): after S2.1-2.2
    // installed the Monomer Collections view AND promoted it to current
    // via grok.shell.v = view (the round-1 belt-and-suspenders pattern),
    // the per-current-view top-menu region is rebuilt against the current
    // view. The Monomer Collections view is NOT a TableView and the Bio
    // top-menu (registered with `meta.targetView: TableView` semantics on
    // the Bio functions) is NOT mounted while a non-table view is
    // current — so `document.querySelector('[name="div-Bio"]')` returns
    // null when S2.3 runs immediately after S2.1-2.2.
    //
    // EMPIRICAL EVIDENCE (live MCP recon, dispatch 2026-06-02 round-2
    // retry): with grok.shell.v.name === 'Manage Monomer Libraries' (a
    // sibling non-TableView), evaluate_script reports
    // bioTopMenuPresent === false; after switching the current view back
    // to the HELM Table via `grok.shell.v = helmTableView` and a brief
    // settle, evaluate_script reports bioTopMenuPresent === true and the
    // full Bio → Manage → Monomer-Libraries menu chain resolves.
    //
    // Deterministic Gate B failure mode across attempts 1-3 cycle
    // 2026-06-02-bio-automate-01 (per the cycle's playwright reports):
    // `TypeError: Cannot read properties of null (reading 'click')` at
    // L668 (the prior div-Bio querySelector) — exactly the null-return
    // case empirically reproduced via MCP recon above. Pre-S2.3
    // normalization (promote HELM Table to current + wait for Bio
    // top-menu to materialize) eliminates the failure mode.
    //
    // The fix mirrors the pre-S3.2 normalization at L755-775 (which
    // already promotes the HELM TableView before saveAllTablesWithProvenance
    // — same need: top-menu / table-anchored affordances are gated on
    // current view being a TableView). Paradigm unchanged: still
    // Playwright top-menu DOM drive on the same class-1 selectors.
    await page.evaluate(async () => {
      const tvs: any[] = Array.from((grok.shell as any).tableViews || []);
      const helm = tvs.find((tv: any) => {
        const df = tv?.dataFrame;
        if (!df) return false;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.some((c: any) => c.semType === 'Macromolecule');
      });
      if (helm) {
        try { (grok.shell as any).v = helm; } catch (_) { /* setter read-only on some builds */ }
      }
      await new Promise((r) => setTimeout(r, 500));
    });
    // Wait for the per-current-view top-menu region to rebuild with the
    // Bio entry visible (TableView-gated affordance). 15s budget is
    // generous against the empirically-observed <1s rebuild latency on
    // the dev server.
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 15_000});
    await softStep('S2.3: Bio | Manage | Monomer Libraries reachable as cross-surface entry point (top-menu DOM drive)', async () => {
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 500));
        const manage = document.querySelector('[name="div-Bio---Manage"]')!;
        manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 500));
        (document.querySelector(
          '[name="div-Bio---Manage---Monomer-Libraries"]') as HTMLElement).click();
      });
      await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Manage Monomer Libraries',
        null, {timeout: 30_000});
      const info = await page.evaluate(() => {
        const v = grok.shell.v;
        const root: any = v?.root;
        const form: any = root?.querySelector?.('.monomer-lib-controls-form');
        return {
          viewName: v?.name ?? null,
          formPresent: !!form,
        };
      });
      // The libraries-view is reachable as a cross-surface entry. We
      // intentionally do NOT assert the collections subpath here — the
      // canonical assertion lives in S2.1-2.2.
      expect((info.viewName || '').toLowerCase()).toContain('monomer librar');
      expect(info.formPresent).toBe(true);
    });
    // ========================================================================
    // Scenario 3 — Save project referencing the collection
    // ========================================================================
    // Scenario 3, Step 1 — Switch back to the HELM TableView (the manage
    // view dispatch in S2.3 docks a separate view on top of the table
    // view). Per sibling bio-lifecycle-monomer-library-spec.ts S3.1: walk
    // grok.shell.tableViews to find the HELM TableView by Macromolecule
    // column presence (NOT grok.shell.tv — that's the current view, which
    // after the manage dispatch is the libraries view).
    await softStep('S3.1: HELM dataset remains open; Macromolecule column renderer is dispatchable post-manage-view', async () => {
      const info = await page.evaluate(() => {
        const tvs: any[] = Array.from((grok.shell as any).tableViews || []);
        let chosen: any = null;
        let chosenMacro: any = null;
        for (const tv of tvs) {
          const df = tv?.dataFrame;
          if (!df) continue;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
          if (macro) {
            chosen = tv;
            chosenMacro = macro;
            break;
          }
        }
        if (!chosen) {
          return {
            hasDf: false,
            hasMacro: false,
            units: null,
            rowCount: 0,
            tableViewCount: tvs.length,
            viewNames: tvs.map((tv: any) => tv?.name ?? '<unnamed>'),
          };
        }
        const df = chosen.dataFrame;
        return {
          hasDf: true,
          hasMacro: !!chosenMacro,
          units: chosenMacro?.meta?.units ?? null,
          rowCount: df.rowCount,
          tableViewCount: tvs.length,
          viewNames: tvs.map((tv: any) => tv?.name ?? '<unnamed>'),
        };
      });
      expect(info.hasDf,
        `expected a TableView with a Macromolecule column among open table views; observed: count=${info.tableViewCount}, names=[${info.viewNames.join(', ')}]`).toBe(true);
      expect(info.hasMacro).toBe(true);
      // Atlas bio.detector: HELM units tag survived the manage / collections
      // view dispatch round trips.
      expect(info.units).toBe('helm');
      expect(info.rowCount).toBeGreaterThan(0);
    });
    // Pre-S3.2: close auxiliary views (Manage Monomer Libraries +
    // Monomer Collections) + bring the HELM TableView forward as
    // grok.shell.v so saveLayout() inside saveAllTablesWithProvenance
    // captures the HELM view's layout, not an auxiliary view's. Mirrors
    // the sibling spec's pre-S3.2 normalization.
    await page.evaluate(async () => {
      const views: any[] = Array.from(grok.shell.views || []);
      for (const name of ['Manage Monomer Libraries', 'Monomer Collections']) {
        const v: any = views.find((x: any) => x?.name === name);
        if (v && typeof v.close === 'function') {
          v.close();
          await new Promise((r) => setTimeout(r, 300));
        }
      }
      const tvs: any[] = Array.from((grok.shell as any).tableViews || []);
      const helm = tvs.find((tv: any) => {
        const df = tv?.dataFrame;
        if (!df) return false;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.some((c: any) => c.semType === 'Macromolecule');
      });
      if (helm && typeof (grok.shell as any).v !== 'undefined') {
        try { (grok.shell as any).v = helm; } catch (_) { /* setter may be read-only on some shell builds */ }
      }
      await new Promise((r) => setTimeout(r, 500));
    });
    // Scenario 3, Step 2 — Save project with Data Sync ON (per scenario:
    // "save the project via the ribbon SAVE button (NOT Ctrl+S); ...
    // Data Sync toggle ON. Cancel the auto-share dialog if it appears.").
    // Per §4.5 Scenario authority + sibling-spec precedent: the Save
    // Project ribbon dialog with Data Sync toggle is platform-wide UI not
    // present in bio.md selector reference. Persistence assertions exercised
    // via the canonical helpers/projects.ts saveAllTablesWithProvenance
    // JS API path (mirrors bio-lifecycle-monomer-library-spec.ts S3.2 +
    // bio-lifecycle-macromolecule-column-spec.ts S3.3 verbatim).
    // saveAllTablesWithProvenance creates a DG.Project, attaches every
    // shell.tables TableInfo, uploads each dataframe via
    // dapi.tables.uploadDataFrame (preserves .script provenance tags),
    // and saves the active TableView's layout — equivalent to the
    // ribbon Save Project ribbon flow with Data Sync ON.
    await softStep('S3.2: save project with provenance (JS API path; mirrors monomer-library / macromolecule-column siblings)', async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.primaryTableInfoId).toBeTruthy();
      // Atlas dep_lifecycle_ops[save_project_with_analysis] write surface
      // succeeded — project lands on the server side with the table-info
      // graph attached.
    });
    // Scenario 3, Step 3 — Close + reopen via JS API; verify:
    //   - HELM dataset restores (atlas
    //     dep_lifecycle_ops[save_project_with_analysis]).
    //   - Macromolecule column re-classifies post-reopen (atlas
    //     bio.detector survives save/reopen).
    //   - getMonomerLibHelper() still resolves the same singleton +
    //     listMonomerCollections() still surfaces the working collection
    //     post-reopen (atlas bio.api.get-monomer-lib-helper +
    //     bio.manage.monomer-collections-app: no FileShare drift caused by
    //     the project save/reopen path per scenario Expected).
    //   - Re-invoking monomerCollectionsApp post-reopen shows the working
    //     collection in the listing (cross-surface state consistency per
    //     scenario Expected).
    await softStep('S3.3: reopen project — HELM survives + collection catalogue stable + Monomer Collections app shows working copy', async () => {
      if (!saved) throw new Error('S3.2 did not produce a saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      // Two-layer post-reopen probe: (1) Macromolecule column survived,
      // (2) MonomerLibManager.listMonomerCollections() still surfaces the
      // working copy, (3) re-invoking the Monomer Collections app shows it
      // in the structural card listing.
      const post = await page.evaluate(async ({fileName}) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        // Library / collection catalogue stability across the reopen
        // boundary.
        const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        let collectionNames: string[] = [];
        try {
          if (typeof helper?.listMonomerCollections === 'function')
            collectionNames = await helper.listMonomerCollections();
        } catch (_) { /* leave empty */ }
        const stem = fileName.replace(/\.json$/, '');
        const hasWorkingCollectionInList = collectionNames.some((n: string) =>
          n === fileName || n === stem || n.includes(stem));
        // Re-dispatch the Monomer Collections app and check the structural
        // card listing (mirrors S2.1-2.2 assertion shape post-reopen).
        // Same belt-and-suspenders view-install pattern as S2.1-2.2: drop
        // the `.root` precondition (undefined post-toJs round-trip),
        // call addView + grok.shell.v setter + workspace walk-back.
        let view: any = null;
        try { view = await (grok as any).functions.call('Bio:monomerCollectionsApp', {}); }
        catch (_) { /* surface via viewName check below */ }
        if (view) {
          try { (grok.shell as any).addView(view); } catch (_) { /* ignore */ }
          try { (grok.shell as any).v = view; } catch (_) { /* ignore */ }
        }
        // Workspace walk-back: promote any 'Monomer Collections' view in
        // shell.views to current if it isn't already.
        try {
          const views: any[] = Array.from((grok.shell as any).views || []);
          const target: any = views.find((vv: any) => vv?.name === 'Monomer Collections');
          if (target && (grok.shell as any).v?.name !== 'Monomer Collections')
            (grok.shell as any).v = target;
        } catch (_) { /* ignore */ }
        // Poll for either install signal:
        //   (a) grok.shell.v.name === 'Monomer Collections', OR
        //   (b) grok.shell.views has a 'Monomer Collections' entry.
        let installed = false;
        for (let i = 0; i < 30; i++) {
          try {
            const cur = (grok.shell as any).v?.name;
            if (cur === 'Monomer Collections') { installed = true; break; }
            const views: any[] = Array.from((grok.shell as any).views || []);
            if (views.some((vv: any) => vv?.name === 'Monomer Collections')) {
              installed = true;
              break;
            }
          } catch (_) { /* keep polling */ }
          await new Promise((r) => setTimeout(r, 200));
        }
        // Extra settle for the FileShare listing + card rebuild.
        await new Promise((r) => setTimeout(r, 2000));
        // Locate the view: prefer current, else walk shell.views.
        let v: any = (grok.shell as any).v;
        if (v?.name !== 'Monomer Collections') {
          const views: any[] = Array.from((grok.shell as any).views || []);
          const found: any = views.find((x: any) => x?.name === 'Monomer Collections');
          if (found) v = found;
        }
        let cardCollectionNames: string[] = [];
        const viewName: string | null = v?.name ?? null;
        if (viewName === 'Monomer Collections') {
          const root: any = v?.root;
          const grid: any = root?.querySelector?.('.monomer-collections-grid');
          const cards = grid ? Array.from(grid.querySelectorAll('.monomer-collection-card')) : [];
          cardCollectionNames = cards.map((c: any) => c.dataset?.collectionName ?? '');
        }
        const hasWorkingCollectionInCards = cardCollectionNames.some((n: string) =>
          n === fileName || n === stem || n.includes(stem));
        const allViewNames: string[] = Array.from((grok.shell as any).views || [])
          .map((x: any) => x?.name ?? '<unnamed>');
        return {
          hasMacro: !!macro,
          units: macro?.meta?.units ?? null,
          helperResolved: helper != null,
          hasWorkingCollectionInList,
          collectionNamesCount: collectionNames.length,
          collectionNames,
          appViewName: viewName,
          currentViewName: (grok.shell as any).v?.name ?? null,
          allViewNames,
          installed,
          cardCount: cardCollectionNames.length,
          cardCollectionNames,
          hasWorkingCollectionInCards,
        };
      }, {fileName: workingCollection});
      // Atlas bio.detector survives reopen (units tag re-attached).
      expect(post.hasMacro).toBe(true);
      expect(post.units).toBe('helm');
      // Atlas bio.api.get-monomer-lib-helper lifecycle:
      // getMonomerLibHelper() still resolves the singleton post-reopen.
      expect(post.helperResolved).toBe(true);
      // The working collection that was discoverable pre-save is still
      // discoverable post-reopen (atlas
      // dep_lifecycle_ops[save_project_with_analysis] doesn't silently
      // drop FileShare collection bindings). The assertable shape is
      // listing stability — the working copy WAS on FileShare before
      // save (S1.1-1.3 addOrUpdateMonomerCollection) and is still on
      // FileShare after reopen (cleanup hasn't run yet — it's in the
      // finally block).
      expect(post.hasWorkingCollectionInList,
        `expected working collection '${workingCollection}' (or stem) in post-reopen listMonomerCollections(); observed: [${post.collectionNames.join(', ')}]`).toBe(true);
      expect(post.collectionNamesCount).toBeGreaterThanOrEqual(1);
      // Atlas bio.manage.monomer-collections-app: cross-surface state
      // consistency post-reopen — the app's card listing also surfaces
      // the working collection (no FileShare drift through the save /
      // reopen boundary). Round-1 retry asserts the LOCATED view name
      // (current OR walked from shell.views).
      expect(post.appViewName,
        `expected to locate a 'Monomer Collections' view (current or in shell.views) post-reopen; observed currentView=${post.currentViewName}, installed=${post.installed}, allViews=[${post.allViewNames.join(', ')}]`).toBe('Monomer Collections');
      expect(post.cardCount).toBeGreaterThanOrEqual(1);
      expect(post.hasWorkingCollectionInCards,
        `expected working collection '${workingCollection}' (or stem) in post-reopen Monomer Collections app cards; observed: [${post.cardCollectionNames.join(', ')}]`).toBe(true);
    });
  } finally {
    // ========================================================================
    // Scenario 4 — Cleanup (runs regardless of earlier failures per scenario
    // Expected: "Cleanup runs in tearDownAll / finally regardless of earlier
    // failures").
    // ========================================================================
    // Step 4.1 — Delete the working collection file via the canonical
    // lib-manager helper (best-effort). Falls back to direct
    // grok.dapi.files.delete on the raw path if the helper handle resolves
    // late.
    if (workingCollectionWritten) {
      await page.evaluate(async ({path, name}) => {
        try {
          const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
          if (helper && typeof helper.deleteMonomerCollection === 'function') {
            try { await helper.deleteMonomerCollection(name); }
            catch (_) { /* fall through to raw delete */ }
          }
          // Belt + suspenders: ensure the file is gone regardless of
          // helper-path success.
          try {
            if (await grok.dapi.files.exists(path))
              await grok.dapi.files.delete(path);
          } catch (_) { /* best effort */ }
        } catch (_) { /* best effort */ }
      }, {path: workingCollectionPath, name: workingCollection.replace(/\.json$/, '')}).catch(() => {});
    }
    // Step 4.2 — Delete the project (best-effort via the canonical helper).
    if (saved) {
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.primaryTableInfoId,
      });
    }
    // Step 4.3 — Close any open Monomer Collections / Manage Monomer
    // Libraries views (best-effort).
    await page.evaluate(async () => {
      // Dismiss any open .d4-dialog.
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      for (const d of dialogs)
        d.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      // Close any auxiliary views still docked.
      try {
        const views = Array.from(grok.shell.views || []);
        for (const name of ['Monomer Collections', 'Manage Monomer Libraries']) {
          const v: any = views.find((x: any) => x?.name === name);
          if (v && typeof v.close === 'function') v.close();
        }
      } catch (_) { /* best effort */ }
    }).catch(() => {});
  }
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
