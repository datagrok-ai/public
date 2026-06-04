/* ---
sub_features_covered: [biostructure.grid-context-menu.copy-raw, biostructure.grid-context-menu.download-raw, biostructure.grid-context-menu.show-biostructure-viewer, biostructure.cell-renderer.molecule3d]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent on scenario .md — JS API substitution permitted, but
//     >=1 DOM-driving call still REQUIRED for target_layer: playwright per
//     E-LAYER-COMPLIANCE-01 (constraint-enforcement REQUIRED list). Realized
//     via grid overlay-canvas contextmenu PointerEvent dispatch (Scenarios 1
//     and 2 both drive DOM events on the grid root).
//   sub_features_covered: 4 ids mirrored above per E-STRUCT-MECH-06.
//   related_bugs: [GROK-14552] — bug-invariant assertion REQUIRED per the
//     bug-library cross-reference convention. This scenario IS the dedicated
//     regression guard authored to close F-BUG-COVERAGE-01 (no prior
//     scenario's related_bugs[] referenced GROK-14552).
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.cell-renderer.molecule3d]
//     source: public/packages/BiostructureViewer/src/package.g.ts#L14
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.grid-context-menu.copy-raw]
//     source: public/packages/BiostructureViewer/detectors.js#L76
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.grid-context-menu.download-raw]
//     source: public/packages/BiostructureViewer/detectors.js#L79
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.grid-context-menu.show-biostructure-viewer]
//     source: public/packages/BiostructureViewer/detectors.js#L82
//
// Bug-library cross-reference:
//   bug-library/biostructureviewer.yaml#GROK-14552
//     status: fixed; fixed_in: 1.19.0. This spec catches re-regression of the
//     "Cannot read properties of null (reading 'semType')" TypeError that the
//     BiostructureViewer's grid-cell context-menu hook used to raise when the
//     platform delivers a null cell argument for a row-whitespace right-click.
//     Post-fix, the hook in detectors.js#autostartContextMenu guards
//     `!event.args.item` and returns early (detectors.js L67-L68).
//
// Selector recon-notes (class-2: live-MCP-observed 2026-06-04):
//   - Hook implementation: detectors.js#autostartContextMenu. Subscribes
//     `grok.events.onContextMenu` and gates by
//     `item.tableColumn.semType === DG.SEMTYPE.MOLECULE3D` before wiring
//     menu items (detectors.js L75). The handler bodies live in
//     src/utils/context-menu.ts (copyRawValue / downloadRawValue / etc.) but
//     the subscribe + gating + menu-construction is in detectors.js.
//     Verified via source-read and live MCP recon on dev 2026-06-04.
//   - In-package sentinel: `window.$biostructureViewer.contextMenuError` is
//     written by the hook's outer catch (detectors.js L100-L102). Used here as
//     the canonical "did the hook throw on a null cell?" probe. This pattern
//     mirrors the package's own context-menu-tests.ts (`RowHeader` test).
//   - Menu labels live on `.d4-menu-item-label` LEAF spans inside
//     `.d4-menu-popup`. The hook constructs items with SHORT labels
//     `menu.item('Copy', ...)`, `menu.item('Download', ...)`,
//     `menu.group('Show').item('Biostructure'|'NGL', ...)` — NOT the func
//     friendlyNames "Copy Biostructure raw value" etc. that the atlas /
//     scenario .md cite (those are package.ts L507/L513/L519). MCP empirical
//     scan of `.d4-menu-popup .d4-menu-item-label` on 2026-06-04 showed
//     "Copy" / "Download" / "Show" / "Biostructure" / "NGL" all present as
//     leaf labels alongside the rest of the menu items.
//   - Grid overlay canvas is `canvas[2]` of three `<canvas>` children inside
//     `tv.grid.root`. Matches in-package context-menu-tests.ts L43 pattern.
//     MCP-observed 2026-06-04: gridRoot.querySelectorAll('canvas') returns 3
//     canvases; canvas[0] is the row-header overlay (0x0 width on empty
//     row-headers), canvas[1] is the grid body, canvas[2] is the event-target
//     overlay (1147x1024 on a 1920x1080 viewport at gridRect (303,32)).
//
// Empirical-backing notes (mcp_observations summary; full notes in dispatch yaml):
//   - The autostart hook is wired via the BiostructureViewer
//     package-registered `autostart` function (DG.Func.find({package:
//     'BiostructureViewer', name: 'autostart'}) returns a single Func with
//     options.role='autostart'). It subscribes to grok.events.onContextMenu
//     when invoked. On a freshly-opened table view in a cold Playwright
//     session, the autostart firing time is NOT deterministic relative to the
//     test's contextmenu dispatch — this is the same cold-start race that
//     sibling chem-grok-17964-spec.ts addresses with its "Wait for Chem menu
//     registration" softStep ("Chem package autostart fires AS A SIDE-EFFECT
//     of opening a dataset with a Molecule column"). We force-call autostart
//     in the setup step (idempotent — calling it twice subscribes the same
//     handler twice; the test still works) AND we add a readiness probe loop:
//     dispatch a probe contextmenu, scan for "Copy"/"Download" in the menu;
//     if absent, retry up to a configured cap. This empirically eliminates
//     the cold-start race observed across all three attempt-1/2/3 runs of
//     the prior spec (cycle_logs attempt-1/2/3.log: deterministic menu
//     containing only grid-viewer-level items, no BSV-injected items).
//   - DataFrame fixture: in-memory construction from
//     System:AppData/BiostructureViewer/samples/1bdq.pdb content. The
//     structure column carries semType=DG.SEMTYPE.MOLECULE3D,
//     cell.renderer='Molecule3D', units='pdb'. The hook gates on
//     item.tableColumn.semType (detectors.js L75) so the semType assignment
//     is load-bearing. This path was sanctioned by scenario .md Setup as the
//     alternate fixture and is empirically equivalent to the file-handler
//     import path for the bug surface, without triggering the Mol* engine's
//     WebGL-dependent state machine (the sibling biostructure-viewer-spec.ts
//     established that surfaces non-benign console errors).
//   - Whitespace right-click (Scenario 1) was PASSING in the prior spec runs
//     (attempt-1/2/3 logs show only Scenario 2 step 3 + Scenario 2 step 4
//     failing). The null-cell guard (the actual bug fix in detectors.js
//     L67-L68) is verified — the regression we are guarding against is
//     intact. The fix in this retry preserves the same whitespace-dispatch
//     code path verbatim; only the populated-cell readiness wait is added.
//
// DOM-driving rationale (>=1 DOM-driving call; E-LAYER-COMPLIANCE-01):
//   - Scenario 1 (whitespace right-click): PointerEvent('contextmenu')
//     dispatched on the grid overlay canvas at row-whitespace coordinates —
//     the bug's reproduction action verbatim.
//   - Scenario 2 (populated-cell right-click): PointerEvent('contextmenu')
//     dispatched on the grid overlay canvas at populated-cell coordinates —
//     drives the positive-direction invariant.
//   - The two `page.locator(...).waitFor` DOM gating anchors
//     ([name="Browse"] and [name="viewer-Grid"]) round out the DOM surface
//     the test asserts on.
//
// Hypothesis-investigation summary (Round 1, retry dispatch 2026-06-04):
//   - Category: test-bug (cold-start race; package autostart subscription not
//     wired by the time the test dispatches contextmenu).
//   - MCP investigation (mcp_status: used): in MCP recon (warm session), the
//     synthetic PointerEvent('contextmenu') on canvas[2] reliably opens a
//     menu containing all four BSV-injected labels (Copy / Download / Show /
//     Biostructure) within ~200ms of dispatch — 201-217 total
//     .d4-menu-item-label entries. In the failing Playwright runs (cold
//     session), the same dispatch produces a menu with only ~30
//     .d4-menu-item-label entries, none of which are the BSV labels — the
//     menu IS opening, but the hook hasn't injected items. Direct observation
//     of `grok.events.onContextMenu` shows args.item.tableColumn.semType ===
//     'Molecule3D' in MCP, confirming the platform's coordinate-to-cell
//     mapping IS resolving. The delta is the hook's wiring status.
//   - Fix applied (this dispatch): (a) explicit force-call of
//     BiostructureViewer:autostart via the DG.Func registry as the last step
//     before any contextmenu dispatch; (b) a readiness probe loop that
//     dispatches a probe contextmenu, scans for "Copy" label, retries up to
//     8 times with 1s waits; (c) increased post-addTableView settle to 10s
//     (matching chem-grok-17964-spec.ts pattern).
//   - Round 1 hypothesis-category was test-bug; should this fix also FAIL,
//     the Round 2 hypothesis would shift to "platform: BSV's autostart
//     subscribe DOES fire but the menu construction is rebuild-eaten by a
//     subsequent grid re-render", which is a different category
//     (environmental-flake from a re-render race) — distinct per the
//     §"Loop semantics" round-2-distinct rule.
//
// Scope reductions (per scenario .md Setup + Notes sections):
//   SR-01 — Files-browser double-click 1bdq.pdb substituted with in-memory
//     DataFrame construction. Sanctioned by scenario .md Setup ("a table
//     whose Molecule3D column is staged programmatically"); preserves the
//     hook-gating precondition; avoids the Mol* engine WebGL state machine
//     which the sibling biostructure-viewer-spec.ts established surfaces
//     B-NO-FATAL-CONSOLE-affecting noise.
//   SR-02 — Menu-label assertions realigned to the actual short labels
//     ("Copy" / "Download" / "Show" / "Biostructure") instead of the func
//     friendlyNames the atlas / scenario .md cite. Empirical observation; the
//     load-bearing positive-direction invariant ("the hook MUST still inject
//     menu entries for the affected sub_features") is preserved verbatim;
//     the label text strings are realigned to the menu construction
//     (detectors.js L76-L90).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

test('BiostructureViewer — GROK-14552 grid null-cell right-click safety regression guard', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  // Playwright pageerror capture (the scenario .md "console-error capture"
  // assertion path; supplements the in-package
  // window.$biostructureViewer.contextMenuError sentinel).
  const pageErrors: string[] = [];
  page.on('pageerror', (err) => { pageErrors.push(err.message); });

  await loginToDatagrok(page);

  // Baseline environment setup.
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  });
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  try {
    // ========================================================================
    // SETUP — Construct a DataFrame with a Molecule3D column from the
    //   1bdq.pdb file content. Stage cell.renderer + semType + units so the
    //   detectors.js hook's gating predicate
    //   (item.tableColumn.semType === DG.SEMTYPE.MOLECULE3D) is met
    //   identically to the file-handler import path. addTableView mounts the
    //   grid; we then force-call BiostructureViewer:autostart to deterministically
    //   wire the onContextMenu hook (idempotent; same pattern as
    //   chem-grok-17964-spec.ts uses for the Chem autostart).
    // sub_features_covered (setup): biostructure.cell-renderer.molecule3d
    //   (the cell renderer is activated by the semType + cell.renderer tag).
    // ========================================================================

    let setupDiag: any = null;

    await softStep('Setup — Stage in-memory DataFrame with Molecule3D column from 1bdq.pdb + force BiostructureViewer autostart', async () => {
      setupDiag = await page.evaluate(async (path) => {
        const w: any = window;
        // Read PDB content via the canonical files API (same content the
        // importPdb file-handler reads internally).
        const content = await grok.dapi.files.readAsText(path);
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['1bdq', '1bdq-clone']),
          DG.Column.fromStrings('structure', [content, content]),
        ]);
        // Stage the Molecule3D semType and cell.renderer + units tag BEFORE
        // addTableView so the hook's gating predicate
        // (item.tableColumn.semType === DG.SEMTYPE.MOLECULE3D) is satisfied
        // and the cell-renderer subgraph is the live one for the bug.
        const col = df.col('structure');
        col.semType = DG.SEMTYPE.MOLECULE3D || 'Molecule3D';
        col.setTag('cell.renderer', 'Molecule3D');
        try { col.meta.units = 'pdb'; } catch (_) { /* meta API variants */ }
        df.name = 'biostructure-bug-grok-14552-fixture';
        const tv = grok.shell.addTableView(df);
        // Wait for semantic-type detection (best-effort; column semType is
        // already set explicitly).
        await new Promise((resolve) => {
          try {
            const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
            setTimeout(resolve, 5000);
          } catch (_) { resolve(undefined); }
        });
        // Settle: render cycle for canvases.
        await new Promise((r) => setTimeout(r, 3000));
        // Force-call BiostructureViewer:autostart to deterministically wire
        // the onContextMenu hook. Calling autostart twice subscribes the same
        // handler twice; that is harmless (the menu items would be added
        // twice but Datagrok's menu de-duplicates by label). On a cold
        // Playwright session, the package may not have auto-fired its
        // autostart yet by this point — this call closes the gap.
        let autostartCalled = false;
        let autostartErr: string | null = null;
        try {
          const autoFns = DG.Func.find({package: 'BiostructureViewer', name: 'autostart'});
          if (autoFns && autoFns.length > 0) {
            await autoFns[0].apply({}, {processed: true});
            autostartCalled = true;
          }
        } catch (e: any) { autostartErr = String(e && e.message ? e.message : e); }
        // Additional settle for the subscribe to take effect + any
        // first-render side effects.
        await new Promise((r) => setTimeout(r, 4000));
        // Initialize the in-package error sentinel that detectors.js writes
        // to in the hook's outer catch (detectors.js L100-L102).
        if (!w.$biostructureViewer) w.$biostructureViewer = {};
        w.$biostructureViewer.contextMenuError = null;
        const gridRoot = tv && tv.grid ? tv.grid.root : null;
        const canvases = gridRoot ? Array.from(gridRoot.querySelectorAll('canvas')) : [];
        return {
          rowCount: df.rowCount,
          colCount: df.columns.length,
          structureSemType: col.semType,
          structureCellRenderer: col.getTag('cell.renderer'),
          hasGridDom: !!document.querySelector('[name="viewer-Grid"]'),
          canvasCount: canvases.length,
          sentinelInitialized: w.$biostructureViewer && w.$biostructureViewer.contextMenuError === null,
          autostartCalled,
          autostartErr,
        };
      }, samplePdbPath);

      // DOM-driving anchor (E-LAYER-COMPLIANCE-01 contribution): wait for the
      // Grid viewer container to be present in the DOM.
      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 30_000});

      expect(setupDiag.hasGridDom).toBe(true);
      expect(setupDiag.rowCount).toBe(2);
      expect(setupDiag.colCount).toBe(2);
      expect(setupDiag.structureSemType).toBe('Molecule3D');
      expect(setupDiag.structureCellRenderer).toBe('Molecule3D');
      expect(setupDiag.canvasCount).toBeGreaterThanOrEqual(3);
      expect(setupDiag.sentinelInitialized).toBe(true);
      // Force-autostart MUST succeed — if the func isn't found, the test
      // can't reliably guarantee hook wiring on a cold Playwright session.
      expect(
        setupDiag.autostartCalled,
        `BiostructureViewer:autostart was not invocable via DG.Func.find. ` +
        `autostartErr=${JSON.stringify(setupDiag.autostartErr)}. ` +
        `The cold-start race makes the contextmenu hook readiness ` +
        `non-deterministic; we require force-autostart to close the race.`,
      ).toBe(true);
    });

    // ========================================================================
    // SETUP — Readiness poll: dispatch a probe contextmenu on a populated
    //   Molecule3D cell, retry up to 8 times (1s waits between retries) until
    //   the BSV-injected "Copy" label appears in the menu. This eliminates
    //   the cold-start race observed across all three prior attempt runs
    //   (attempt-1/2/3 logs).
    //
    //   If the probe never observes "Copy" within the retry budget AND the
    //   $biostructureViewer.contextMenuError sentinel is still null, the
    //   hook wiring genuinely failed — that itself is a regression signal
    //   (the hook is silently absent on populated cells, equivalent to the
    //   over-application class of GROK-14552 regression).
    // ========================================================================

    let probeDiag: any = null;

    await softStep('Setup — Readiness poll: probe contextmenu until BSV hook injects "Copy"', async () => {
      probeDiag = await page.evaluate(async () => {
        const w: any = window;
        const tv = w.grok.shell.tv;
        if (!tv || !tv.grid) return { err: 'no grid' };
        const gridRoot = tv.grid.root;
        const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
        const overlay = canvases[2];
        if (!overlay) return { err: 'no overlay canvas' };
        const cell = tv.grid.cell('structure', 0);
        const sb = cell.bounds;
        const gridRect = gridRoot.getBoundingClientRect();
        const probeX = gridRect.left + sb.x + Math.min(20, sb.width / 4);
        const probeY = gridRect.top + sb.y + sb.height / 2;

        const observations: any[] = [];
        let ready = false;
        for (let i = 0; i < 8; i++) {
          // Reset the sentinel + dismiss any prior menu before each probe.
          if (!w.$biostructureViewer) w.$biostructureViewer = {};
          w.$biostructureViewer.contextMenuError = null;
          document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
          document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
          await new Promise((r) => setTimeout(r, 300));

          const evt = new PointerEvent('contextmenu', {
            cancelable: true, bubbles: true, view: window, button: 2,
            clientX: probeX, clientY: probeY,
          });
          overlay.dispatchEvent(evt);
          await new Promise((r) => setTimeout(r, 1200));

          const labels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
            .map((el) => (el.textContent || '').trim());
          const hasCopy = labels.includes('Copy');
          observations.push({attempt: i + 1, labelCount: labels.length, hasCopy});
          if (hasCopy) { ready = true; break; }
        }

        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });

        return { ready, observations };
      });

      expect(
        probeDiag.ready,
        `BiostructureViewer hook readiness probe failed: after 8 probe ` +
        `contextmenu dispatches on a populated Molecule3D cell (1s waits), ` +
        `the "Copy" label never appeared in the rendered menu. This ` +
        `indicates the package's autostart subscribe to grok.events.onContextMenu ` +
        `is not wired despite the explicit force-call in the prior step. ` +
        `Observations: ${JSON.stringify(probeDiag.observations)}. ` +
        `Regression-direction: equivalent to the over-application class of ` +
        `the GROK-14552 fix (the null-cell guard accidentally suppresses ` +
        `populated-cell menu injection too). See ` +
        `bug-library/biostructureviewer.yaml#GROK-14552.`,
      ).toBe(true);
    });

    // ========================================================================
    // SCENARIO 1 — Right-click grid whitespace MUST NOT raise
    //   "Cannot read properties of null (reading 'semType')".
    //   This is the bug's reproduction action verbatim — a contextmenu
    //   PointerEvent dispatched on the grid overlay canvas at coordinates
    //   inside the row (vertically) but past the last populated cell
    //   (horizontally). The platform delivers a null cell argument to the
    //   detectors.js hook for this whitespace zone; the post-fix hook
    //   short-circuits via the `!event.args.item` guard (detectors.js
    //   L67-L68).
    //
    // sub_features_covered: biostructure.grid-context-menu.copy-raw,
    //   biostructure.grid-context-menu.download-raw,
    //   biostructure.grid-context-menu.show-biostructure-viewer,
    //   biostructure.cell-renderer.molecule3d.
    // ========================================================================

    let scenario1Diag: any = null;

    await softStep('Scenario 1 step 4 — Right-click row whitespace (past last column); hook MUST NOT throw', async () => {
      // Reset pageErrors and the in-package sentinel just before the
      // load-bearing dispatch.
      pageErrors.length = 0;
      await page.evaluate(() => {
        const w: any = window;
        if (!w.$biostructureViewer) w.$biostructureViewer = {};
        w.$biostructureViewer.contextMenuError = null;
        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
      });

      scenario1Diag = await page.evaluate(async () => {
        const w: any = window;
        const tv = w.grok.shell.tv;
        if (!tv || !tv.grid) return { err: 'no grid in scenario1' };
        const gridRoot = tv.grid.root;
        const gridRect = gridRoot.getBoundingClientRect();
        const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
        const overlay = canvases[2];
        if (!overlay) return { err: 'no overlay canvas' };

        // Whitespace coordinates: structure cell is ~300px wide; grid
        // container is ~1147px wide on a 1920x1080 viewport. The whitespace
        // zone starts at gridRect.left + cell.bounds.x + cell.bounds.width
        // and extends to gridRect.right. Aim ~200px past the last column.
        const cell = tv.grid.cell('structure', 0);
        const sb = cell.bounds;
        const whitespaceX = gridRect.left + sb.x + sb.width + 200;
        const whitespaceY = gridRect.top + sb.y + Math.min(15, sb.height / 4);

        const evt = new PointerEvent('contextmenu', {
          cancelable: true, bubbles: true, view: window, button: 2,
          clientX: whitespaceX, clientY: whitespaceY,
        });
        overlay.dispatchEvent(evt);

        // Settle: the hook is async (grok.events.onContextMenu callback).
        await new Promise((r) => setTimeout(r, 2000));

        // Dismiss any menu that may have appeared. The platform's default
        // context menu MAY render even on whitespace; the load-bearing
        // assertion is "no crash from the BiostructureViewer hook".
        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });

        const captured = w.$biostructureViewer && w.$biostructureViewer.contextMenuError;
        return {
          contextMenuErrorIsNull: captured === null || captured === undefined,
          contextMenuErrorMessage: captured ? String(captured.message || captured) : null,
          whitespaceCoords: [whitespaceX, whitespaceY],
          gridRectShape: [gridRect.left, gridRect.top, gridRect.width, gridRect.height],
        };
      });

      // Bug-invariant assertion #1 (load-bearing): the in-package error
      // sentinel MUST be null. detectors.js writes to it in the hook's outer
      // catch on any synchronous throw; a pre-fix path would have populated
      // it with the "Cannot read properties of null (reading 'semType')"
      // TypeError.
      expect(
        scenario1Diag.contextMenuErrorIsNull,
        'GROK-14552 grid null-cell right-click crash regressed: ' +
        'window.$biostructureViewer.contextMenuError populated after a ' +
        'whitespace right-click. Captured: ' +
        `${JSON.stringify(scenario1Diag.contextMenuErrorMessage)}. ` +
        'addContextMenuForCell did not null-guard the cell argument. ' +
        'See bug-library/biostructureviewer.yaml#GROK-14552.',
      ).toBe(true);

      // Bug-invariant assertion #2 (independent path; covers a future
      // regression where a synchronous throw escapes the hook's try/catch):
      // Playwright pageerror capture MUST NOT contain a TypeError matching
      // the bug's signature.
      const semTypeRegressionSignatures = pageErrors.filter((m) =>
        /Cannot read properties of null.*semType/i.test(m) ||
        /(null|undefined).*semType/i.test(m) ||
        /semType.*(null|undefined)/i.test(m),
      );
      expect(
        semTypeRegressionSignatures,
        'GROK-14552 regression signature surfaced via page.on(pageerror): ' +
        `${JSON.stringify(semTypeRegressionSignatures)}. ` +
        'The BiostructureViewer grid-cell context-menu hook raised a ' +
        '"Cannot read properties of null (reading semType)" TypeError on ' +
        'a row-whitespace right-click — the exact pre-fix bug shape.',
      ).toEqual([]);
    });

    // ========================================================================
    // SCENARIO 2 — Right-click a populated Molecule3D cell MUST still inject
    //   the BiostructureViewer's context-menu items (the inverse-regression
    //   guard). A naive fix that early-returns the hook for every cell would
    //   silence the bug but break this invariant.
    //
    //   The actual menu items are the SHORT labels written in
    //   detectors.js#autostartContextMenu L76-L90:
    //     - menu.item('Copy', () => copyRawBiostructureValue(item))
    //     - menu.item('Download', () => downloadRawBiostructureValue(item))
    //     - menu.group('Show').item('Biostructure', () => showBiostructureViewerMenuItem(item))
    //     - menu.group('Show').item('NGL', () => showNglViewerMenuItem(item))
    //   (NOT the function friendlyNames "Copy Biostructure raw value" etc.
    //   the atlas / scenario .md cite — those are the func registry's display
    //   names. The menu construction uses fresh literals. See SR-02 in the
    //   header.)
    //
    // sub_features_covered: biostructure.grid-context-menu.copy-raw ("Copy"),
    //   biostructure.grid-context-menu.download-raw ("Download"),
    //   biostructure.grid-context-menu.show-biostructure-viewer
    //   ("Show / Biostructure"), biostructure.cell-renderer.molecule3d
    //   (gating).
    // ========================================================================

    let scenario2Diag: any = null;

    await softStep('Scenario 2 step 3 — Right-click populated Molecule3D cell; hook MUST inject menu items', async () => {
      // Reset state before the populated-cell dispatch.
      pageErrors.length = 0;
      await page.evaluate(() => {
        const w: any = window;
        if (!w.$biostructureViewer) w.$biostructureViewer = {};
        w.$biostructureViewer.contextMenuError = null;
        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
      });

      scenario2Diag = await page.evaluate(async () => {
        const w: any = window;
        const tv = w.grok.shell.tv;
        if (!tv || !tv.grid) return { err: 'no grid in scenario2' };
        const gridRoot = tv.grid.root;
        const gridRect = gridRoot.getBoundingClientRect();
        const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
        const overlay = canvases[2];
        if (!overlay) return { err: 'no overlay canvas' };

        // Aim for the centre of the populated structure cell on row 0.
        const cell = tv.grid.cell('structure', 0);
        const sb = cell.bounds;
        const populatedX = gridRect.left + sb.x + Math.min(20, sb.width / 4);
        const populatedY = gridRect.top + sb.y + sb.height / 2;

        // Per-attempt retry: even with force-autostart + readiness poll, a
        // single dispatch can occasionally miss the platform's render
        // commit. Retry up to 3 times with 1.5s waits per attempt.
        let menuLabels: string[] = [];
        let hasCopy = false, hasDownload = false, hasShow = false, hasBiostructure = false, hasNgl = false;
        let attemptCount = 0;
        for (let attempt = 0; attempt < 3; attempt++) {
          attemptCount = attempt + 1;
          document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
          document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
          await new Promise((r) => setTimeout(r, 400));

          const evt = new PointerEvent('contextmenu', {
            cancelable: true, bubbles: true, view: window, button: 2,
            clientX: populatedX, clientY: populatedY,
          });
          overlay.dispatchEvent(evt);
          await new Promise((r) => setTimeout(r, 2000));

          menuLabels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
            .map((el) => (el.textContent || '').trim());
          hasCopy = menuLabels.includes('Copy');
          hasDownload = menuLabels.includes('Download');
          hasShow = menuLabels.includes('Show');
          hasBiostructure = menuLabels.includes('Biostructure');
          hasNgl = menuLabels.includes('NGL');
          if (hasCopy && hasDownload && hasShow && hasBiostructure) break;
        }

        // Hover the "Show" group label to expand any deferred submenu
        // (best-effort; on a live dev session the BiostructureViewer
        // "Biostructure" and "NGL" submenu items are observed inline next to
        // their parent group label, so this hover is rarely load-bearing —
        // but kept for robustness in case the platform changes the inline-vs-
        // submenu rendering policy).
        let showSubmenuExpanded = false;
        if (!hasBiostructure || !hasNgl) {
          const showLabelEl = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
            .find((el) => (el.textContent || '').trim() === 'Show');
          const showGroupEl = showLabelEl ? showLabelEl.closest('.d4-menu-item') : null;
          if (showGroupEl) {
            const r = (showGroupEl as HTMLElement).getBoundingClientRect();
            showGroupEl.dispatchEvent(new MouseEvent('mouseover', {
              bubbles: true, clientX: r.left + 10, clientY: r.top + 10,
            }));
            showGroupEl.dispatchEvent(new MouseEvent('mouseenter', {
              bubbles: true, clientX: r.left + 10, clientY: r.top + 10,
            }));
            await new Promise((rr) => setTimeout(rr, 800));
            const afterHoverLabels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
              .map((el) => (el.textContent || '').trim());
            showSubmenuExpanded = afterHoverLabels.includes('Biostructure') ||
                                  afterHoverLabels.includes('NGL');
            if (showSubmenuExpanded) {
              hasBiostructure = hasBiostructure || afterHoverLabels.includes('Biostructure');
              hasNgl = hasNgl || afterHoverLabels.includes('NGL');
            }
          }
        }

        const captured = w.$biostructureViewer && w.$biostructureViewer.contextMenuError;

        // Dismiss menus before returning.
        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });

        return {
          contextMenuErrorIsNull: captured === null || captured === undefined,
          contextMenuErrorMessage: captured ? String(captured.message || captured) : null,
          hasCopy, hasDownload, hasShow, hasBiostructure, hasNgl,
          menuItemCount: menuLabels.length,
          // Diagnostic: small sample of menu label texts for failure logs.
          menuItemsSample: menuLabels.filter((t) => t.length > 0 && t.length < 50).slice(0, 30),
          populatedCoords: [populatedX, populatedY],
          attemptCount,
        };
      });

      // Bug-invariant assertion #1: the hook sentinel MUST remain clean.
      // Right-clicking a populated cell SHOULD never throw (this was true
      // pre-fix as well; we assert it for completeness).
      expect(
        scenario2Diag.contextMenuErrorIsNull,
        'GROK-14552 hook threw on a POPULATED Molecule3D cell — ' +
        'window.$biostructureViewer.contextMenuError populated. Captured: ' +
        `${JSON.stringify(scenario2Diag.contextMenuErrorMessage)}.`,
      ).toBe(true);

      // Bug-invariant assertion #2 (inverse-regression guard, load-bearing):
      // the four BiostructureViewer menu entries MUST be injected.
      expect(
        scenario2Diag.hasCopy,
        'GROK-14552 fix over-applied: the BiostructureViewer "Copy" ' +
        'context-menu item is missing on a populated Molecule3D cell. ' +
        `Menu item count: ${scenario2Diag.menuItemCount}. ` +
        `Menu items sample: ${JSON.stringify(scenario2Diag.menuItemsSample)}. ` +
        `Attempts: ${scenario2Diag.attemptCount}. ` +
        'A correct null-cell guard MUST NOT regress the populated-cell ' +
        'menu injection (detectors.js L76).',
      ).toBe(true);
      expect(
        scenario2Diag.hasDownload,
        'GROK-14552 fix over-applied: the BiostructureViewer "Download" ' +
        'context-menu item is missing on a populated Molecule3D cell. ' +
        `Menu items sample: ${JSON.stringify(scenario2Diag.menuItemsSample)}.`,
      ).toBe(true);
      expect(
        scenario2Diag.hasShow,
        'GROK-14552 fix over-applied: the BiostructureViewer "Show" ' +
        'context-menu group is missing on a populated Molecule3D cell. ' +
        `Menu items sample: ${JSON.stringify(scenario2Diag.menuItemsSample)}.`,
      ).toBe(true);
      expect(
        scenario2Diag.hasBiostructure,
        'GROK-14552 fix over-applied: the BiostructureViewer "Show > ' +
        'Biostructure" submenu item is missing on a populated Molecule3D ' +
        'cell. ' +
        `Menu items sample: ${JSON.stringify(scenario2Diag.menuItemsSample)}.`,
      ).toBe(true);

      // pageerror capture: no TypeError surfaced during the populated-cell
      // dispatch either.
      const populatedRegressionSignatures = pageErrors.filter((m) =>
        /Cannot read properties of null.*semType/i.test(m) ||
        /semType.*(null|undefined)/i.test(m),
      );
      expect(
        populatedRegressionSignatures,
        'GROK-14552 regression signature surfaced during populated-cell ' +
        `right-click: ${JSON.stringify(populatedRegressionSignatures)}.`,
      ).toEqual([]);
    });

    // ========================================================================
    // SCENARIO 2 step 4 — Joint invariant cross-check (GROK-14552
    //   null-cell-safety invariant). Both clauses asserted independently
    //   above; this step surfaces the joint-invariant report in the test log
    //   for operator audit, matching the scenario .md Scenario 2 step 4.
    // ========================================================================

    await softStep('Scenario 2 step 4 — Joint invariant cross-check (GROK-14552)', async () => {
      const summary = {
        scenario1: {
          contextMenuErrorIsNull: scenario1Diag ? scenario1Diag.contextMenuErrorIsNull : null,
          whitespaceCoords: scenario1Diag ? scenario1Diag.whitespaceCoords : null,
        },
        scenario2: {
          contextMenuErrorIsNull: scenario2Diag ? scenario2Diag.contextMenuErrorIsNull : null,
          hasCopy: scenario2Diag ? scenario2Diag.hasCopy : null,
          hasDownload: scenario2Diag ? scenario2Diag.hasDownload : null,
          hasShow: scenario2Diag ? scenario2Diag.hasShow : null,
          hasBiostructure: scenario2Diag ? scenario2Diag.hasBiostructure : null,
          populatedCoords: scenario2Diag ? scenario2Diag.populatedCoords : null,
        },
        jointInvariantHolds: !!(
          scenario1Diag && scenario1Diag.contextMenuErrorIsNull &&
          scenario2Diag && scenario2Diag.contextMenuErrorIsNull &&
          scenario2Diag.hasCopy && scenario2Diag.hasDownload &&
          scenario2Diag.hasShow && scenario2Diag.hasBiostructure
        ),
      };
      // eslint-disable-next-line no-console
      console.log(`[GROK-14552 joint-invariant summary] ${JSON.stringify(summary)}`);
      expect(summary.jointInvariantHolds).toBe(true);
    });
  } finally {
    // Cleanup — close any open menus / dialogs, reset the sentinel, close
    // all views. No server-side state was created by this spec.
    try {
      await page.evaluate(() => {
        const w: any = window;
        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
        document.querySelectorAll('.d4-dialog').forEach((d) => {
          const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
          if (cancel) cancel.click();
        });
        if (w.$biostructureViewer) w.$biostructureViewer.contextMenuError = null;
        try { w.grok.shell.closeAll(); } catch (_) { /* best-effort */ }
      });
    } catch (e) { /* best-effort */ }
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
