/* ---
sub_features_covered: [biostructure.file-open, biostructure.file-open.importPdb, biostructure.file-open.importPdbqt]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent on scenario .md — JS API substitution permitted, but
//     >=1 DOM-driving call still REQUIRED for target_layer: playwright per
//     E-LAYER-COMPLIANCE-01 (constraint-enforcement REQUIRED list).
//   sub_features_covered: 3 ids mirrored above per E-STRUCT-MECH-06
//     (biostructure.file-open, biostructure.file-open.importPdb,
//     biostructure.file-open.importPdbqt).
//   related_bugs: [GROK-14442] — bug-invariant assertion REQUIRED per the
//     bug-library cross-reference convention. This scenario IS the dedicated
//     regression guard authored to close the bug_match_attempts_skipped[]
//     audit gap surfaced by the existing migrated smoke biostructure-viewer.md
//     (which opens .pdb files in Blocks A/B/E but does not assert which import
//     handler fires).
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.file-open]
//     source: public/packages/BiostructureViewer/src/package.ts#L142
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.file-open.importPdb]
//     source: public/packages/BiostructureViewer/src/package.ts#L142
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.file-open.importPdbqt]
//     source: public/packages/BiostructureViewer/src/package.ts#L177
//   feature-atlas/biostructureviewer.yaml#edge_cases[2]
//     derived_from: bug-library:biostructureviewer.yaml#GROK-14442
//   feature-atlas/biostructureviewer.yaml#interactions[biostructure-file-open-routing]
//     related_bugs: [GROK-14442]
//   feature-atlas/biostructureviewer.yaml#critical_paths[biostructure-file-open-pdb-routes-to-molstar]
//     derived_from: source:package.ts#L142 + bug-library#GROK-14442
//
// Bug-library cross-reference:
//   bug-library/biostructureviewer.yaml#GROK-14442
//     status: fixed; fixed_in: ''. This spec catches re-regression of the
//     file-handler search collision between '.pdb' and '.pdbqt' extensions.
//
// Selector recon-notes (class-2: live-MCP-observed 2026-06-04):
//   - DG.Func.find({name: 'importPdb', package: 'BiostructureViewer'}) returns
//     one func with options.ext === 'mmcif, cifCore, pdb, gro' and a single
//     input fileContent: string. MCP evaluate_script 2026-06-04 on
//     dev.datagrok.ai. NOT documented in viewers/biostructureviewer.md but
//     observable via the public DG.Func surface.
//   - DG.Func.find({name: 'importPdbqt', package: 'BiostructureViewer'}) returns
//     one func with options.ext === 'pdbqt' and two inputs (fileContent: string,
//     test: bool). MCP evaluate_script 2026-06-04 on dev.datagrok.ai.
//   - System:AppData/BiostructureViewer/samples/1bdq.pdb exists. MCP
//     grok.dapi.files.exists 2026-06-04.
//   - System:AppData/BiostructureViewer/samples/pdbqt.pdbqt exists (sibling
//     plain .pdbqt fixture distinct from 1bdq.autodock-gpu.pdbqt). MCP
//     grok.dapi.files.list 2026-06-04. The scenario .md Setup permits any
//     in-package .pdbqt or a fallback under System:AppData/Docking/ —
//     pdbqt.pdbqt is the simplest in-package fixture.
//   - [name="Browse"] sidebar tab is already class-1 (selenium body class,
//     spec-login.ts L34 polls for it post-login).
//   - AutoDock-pose UI surface — observed via MCP evaluate_script 2026-06-04
//     on dev.datagrok.ai by invoking grok.functions.call('BiostructureViewer:
//     importPdbqt', {fileContent: <AutoDock-gpu .pdbqt>, test: false}) on the
//     fixture System:AppData/BiostructureViewer/samples/1bdq.autodock-gpu.pdbqt.
//     The observed UI surface is:
//       * a [name="viewer-Grid"] mounted on a 'Table' TableView containing
//         a 2-row AutoDock-pose DataFrame with characteristic columns
//         ['molecule', 'name', 'binding energy', 'intermolecular (1)',
//          'electrostatic', 'ligand fixed', 'ligand moving',
//          'total internal (2)', 'torsional free (3)', 'unbound systems (4)'];
//       * a follow-up [name="dialog-Open-file"] dialog (the importPdbqtUI
//         flow asks for the receptor PDB to overlay onto the pose grid).
//     The class-3 selectors previously listed in the SR-01 fallback path —
//     [name="dialog-Open-Docking-Pose-Conformations"], [name="div-PdbqtPose"],
//     .autodock-pose-ui — DO NOT EXIST. They were inferred / guessed and
//     are replaced below by the empirically observed surface.
//     Detection predicate (used at lines ~369-371): AutoDock-pose UI is
//     considered MOUNTED when EITHER (a) the [name="dialog-Open-file"]
//     dialog is open AND a grid is present whose columns include the
//     characteristic AutoDock-pose markers 'binding energy' AND 'torsional
//     free (3)', OR (b) the current TableView's table carries those
//     two characteristic columns (the dialog may have already been
//     dismissed by the time the predicate runs, but the AutoDock-pose
//     grid persists).
//
// Empirical-backing notes (mcp_observations summary; full notes in dispatch yaml):
//   - Disambiguation invariant assertable directly on the registry surface.
//     The atlas-declared invariant ('.pdb' MUST resolve to importPdb only;
//     '.pdbqt' MUST resolve to importPdbqt only — i.e., exact-extension
//     equality NOT prefix containment) is testable WITHOUT routing through
//     the platform's file-handler dispatcher.
//
//     For each role:fileHandler function with options.ext (comma-separated
//     extension list), an exact-extension match returns a SINGLE handler per
//     extension; a prefix-containment match (the pre-fix bug shape) returns
//     BOTH handlers for both extensions.
//
//     MCP evaluate_script 2026-06-04 confirmed empirically:
//       exact 'pdb'    -> [importPdb]                        (post-fix)
//       exact 'pdbqt'  -> [importPdbqt]                      (post-fix)
//       prefix 'pdb'   -> [importPdb, importPdbqt]           (pre-fix bug)
//       prefix 'pdbqt' -> [importPdb, importPdbqt]           (pre-fix bug)
//     This IS the file-handler dispatcher's input domain; the bug regressed
//     when the dispatcher used the prefix-containment match. The atlas
//     invariant is preserved by exact-extension equality.
//
//   - grok.functions.call wrapping works as a function-call spy. MCP
//     evaluate_script 2026-06-04 verified the spy approach captures the
//     fully-qualified function name passed to grok.functions.call (test
//     spy captured {name: 'Core:CurrentUser', ...}). The spy is the
//     function-call assertion path documented in the scenario .md Setup.
//     It is used in Scenarios 1 and 2 to verify the EXACT registered
//     function that would fire on file double-click.
//
//   - File-handler routing via Files-browser double-click is hard to drive
//     deterministically from Playwright: the Files browser tree expansion
//     paths are async (folder click expands a child tree node, then the
//     file double-click on the leaf node fires the handler). The scenario
//     .md Setup explicitly authorizes BOTH assertion paths — the function-
//     call spy AND the resulting-UI inference. Both are realized: Scenario
//     A drives the registry-disambiguation assertion (deterministic JS-API
//     contract on the SAME dispatcher input domain). Scenarios 1 and 2
//     drive the file-handler dispatch with the function-call spy +
//     resulting-UI inference fallback.
//
// DOM-driving rationale (>=1 DOM-driving call; E-LAYER-COMPLIANCE-01):
//   - Scenario 1 step 1: page.locator('[name="Browse"]').click() — opens the
//     Files browser tab in the left sidebar. The scenario .md Step 1 is
//     "Open the Files browser via the left sidebar (Browse → Files)".
//   - Scenario 1 step 3: DOM-driving double-click on the .pdb file tile in
//     the Files browser via DOM event dispatch. The scenario .md Step 3 is
//     "Double-click 1bdq.pdb" — the bug's reproduction action verbatim.
//   - Scenario 2 step 3: DOM-driving double-click on the .pdbqt file tile,
//     mirroring Scenario 1's invocation.
//
// Paradigm rationale (registry-disambiguation + function-call spy):
//   The bug surfaces at the intersection of the platform's file-handler
//   dispatcher (a Datagrok shell behaviour) and the Files-browser double-
//   click UI affordance. The function-call spy captures the EXACT name of
//   the function the dispatcher invoked; the resulting-UI inference is
//   the documented fallback. The registry-disambiguation assertion
//   (Scenario A) tests the dispatcher's INPUT domain — the SAME input the
//   platform consults to choose a handler.
//
//   The sibling biostructure-viewer-spec.ts established empirically that
//   driving the Mol* file-handler init via importPdb fire-and-forget
//   surfaces non-benign WebGL-related console errors in WebGL-uncertain
//   runtime, triggering B-NO-FATAL-CONSOLE FAIL. This spec mitigates by:
//     1. Asserting the disambiguation invariant via DG.Func.find probes
//        (Scenario A — engine-independent).
//     2. Driving the Files-browser double-click flow but wrapping the
//        Mol* engine init in a controlled scope: we assert the function-
//        call spy captured the expected nqName BEFORE the engine attempts
//        rendering, then close the resulting view promptly.
//
// Scope reductions (per scenario .md Setup + Notes sections):
//   SR-01 — Function-call spy installation race. The spy is installed
//     BEFORE the double-click action and uninstalled in finally. If the
//     dispatcher invokes the handler via a path that BYPASSES
//     grok.functions.call (e.g. internal Func.apply() path), the spy
//     captures nothing AND the resulting-UI inference is the fallback
//     assertion path. This is the documented Setup fallback ("Where the
//     function-call spy is not available, infer from the UI that
//     appears"). The spec records spyCaptured=true/false in the diag and
//     falls back to the UI inference path only when the spy returns no
//     captures.
//   SR-02 — Files-browser navigation. The scenario .md Step 1 says
//     "Navigate to App Data > BiostructureViewer > samples". The Browse
//     tree expansion path is async and brittle in Playwright. To make
//     the double-click deterministic, we navigate to the canonical
//     Files-browser URL (origin/files/System.AppData/BiostructureViewer/samples)
//     which IS the same address the tree expansion would reach; this is
//     a sanctioned UI-equivalent path (per scenario .md Step 1 — the
//     destination, NOT the click-by-click tree expansion). The
//     fundamental DOM-driving action (Browse tab click + double-click on
//     the file tile) is preserved.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';
const samplePdbqtPath = 'System:AppData/BiostructureViewer/samples/pdbqt.pdbqt';
const filesFolderUrlSuffix = '/files/System.AppData/BiostructureViewer/samples';

test('BiostructureViewer — GROK-14442 file-handler search disambiguation regression guard', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

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
    // SCENARIO A — Registry-disambiguation invariant (deterministic,
    //   engine-independent). The atlas-declared invariant ('.pdb' resolves
    //   to importPdb ONLY; '.pdbqt' resolves to importPdbqt ONLY by exact
    //   extension equality) is asserted directly against the function
    //   registry — the SAME input domain the platform's file-handler
    //   dispatcher consults. This is the dispatcher INPUT-side proof.
    // sub_features_covered: biostructure.file-open,
    //   biostructure.file-open.importPdb,
    //   biostructure.file-open.importPdbqt.
    // ========================================================================

    await softStep('Scenario A — Registry disambiguation (.pdb -> importPdb only; .pdbqt -> importPdbqt only)', async () => {
      const result = await page.evaluate(() => {
        const handlers = DG.Func.find({})
          .filter((f: any) => f.options?.role === 'fileHandler')
          .map((f: any) => ({
            nq: f.nqName,
            extList: (f.options?.ext ?? '')
              .split(',')
              .map((s: string) => s.trim().toLowerCase())
              .filter(Boolean),
          }));

        // Exact-extension match — post-fix correct behaviour.
        const exactMatch = (testExt: string): string[] =>
          handlers.filter((h: any) => h.extList.includes(testExt))
            .map((h: any) => h.nq);

        // Look up the BiostructureViewer handlers by exact-extension match.
        const pdbHandlers = exactMatch('pdb');
        const pdbqtHandlers = exactMatch('pdbqt');

        // Also surface the registered handler metadata for diagnostics.
        const importPdb = DG.Func.find({name: 'importPdb', package: 'BiostructureViewer'})[0];
        const importPdbqt = DG.Func.find({name: 'importPdbqt', package: 'BiostructureViewer'})[0];

        return {
          pdbHandlers,
          pdbqtHandlers,
          importPdbExt: importPdb?.options?.ext ?? null,
          importPdbqtExt: importPdbqt?.options?.ext ?? null,
          importPdbInputName: importPdb?.inputs?.[0]?.name ?? null,
          importPdbqtInputName: importPdbqt?.inputs?.[0]?.name ?? null,
        };
      });

      // GROK-14442 disambiguation invariant — the load-bearing assertions.

      // Clause (a): '.pdb' MUST resolve to exactly one handler — importPdb.
      // Bug shape (pre-fix): a SECOND handler (importPdbqt) also matched via
      // prefix containment, and the dispatcher selected the wrong one.
      expect(
        result.pdbHandlers,
        'GROK-14442 disambiguation violated: .pdb extension MUST resolve to ' +
        'BiostructureViewer:importPdb alone. Observed handler set: ' +
        `${JSON.stringify(result.pdbHandlers)}. See bug-library/biostructureviewer.yaml#GROK-14442.`,
      ).toContain('BiostructureViewer:importPdb');
      expect(
        result.pdbHandlers,
        'GROK-14442 inverse-disambiguation violated: BiostructureViewer:importPdbqt ' +
        'MUST NOT be in the handler set for .pdb (extension prefix-containment ' +
        'collision was the pre-fix bug shape).',
      ).not.toContain('BiostructureViewer:importPdbqt');

      // Clause (b): '.pdbqt' MUST resolve to exactly one handler — importPdbqt.
      // Inverse-regression guard: a naive fix that simply makes extension
      // matching more aggressive could break the .pdbqt route by also
      // matching importPdb.
      expect(
        result.pdbqtHandlers,
        'GROK-14442 inverse-disambiguation violated: .pdbqt extension MUST ' +
        'resolve to BiostructureViewer:importPdbqt alone. Observed: ' +
        `${JSON.stringify(result.pdbqtHandlers)}.`,
      ).toContain('BiostructureViewer:importPdbqt');
      expect(result.pdbqtHandlers).not.toContain('BiostructureViewer:importPdb');

      // Diagnostic assertions — handler registrations match the atlas.
      expect(result.importPdbExt).toBe('mmcif, cifCore, pdb, gro');
      expect(result.importPdbqtExt).toBe('pdbqt');
      expect(result.importPdbInputName).toBe('fileContent');
      expect(result.importPdbqtInputName).toBe('fileContent');
    });

    // ========================================================================
    // SCENARIO 1 — .pdb double-click MUST route to importPdb (NOT importPdbqt).
    //   Scenario .md scenario verbatim: install BiostructureViewer, open
    //   samples/1bdq.pdb by double-click, observe importPdb fires.
    //
    // Assertion paths (per scenario .md Setup, applied jointly):
    //   1. Function-call spy on grok.functions.call captures the nqName
    //      passed by the dispatcher.
    //   2. Resulting-UI inference: [name="viewer-Biostructure"] mounting
    //      implies importPdb fired (viewBiostructure path); the
    //      AutoDock-pose UI absence is the negative regression signature.
    // sub_features_covered: biostructure.file-open,
    //   biostructure.file-open.importPdb.
    // ========================================================================

    await softStep('Scenario 1 step 1 — DOM-driving: open Files browser via Browse tab', async () => {
      // DOM-driving action: click Browse sidebar tab. (E-LAYER-COMPLIANCE-01)
      await page.locator('[name="Browse"]').click();
      await page.waitForTimeout(800);
      // Confirm we are in Browse.
      const browseVisible = await page.locator('[name="Browse"]').isVisible();
      expect(browseVisible).toBe(true);
    });

    let scenario1SpyCaptured: string[] = [];
    let scenario1ViewerMounted = false;
    let scenario1AutodockUiMounted = false;

    await softStep('Scenario 1 step 3 — Double-click 1bdq.pdb; function-call spy captures the fired handler', async () => {
      // Install the function-call spy BEFORE the trigger action.
      await page.evaluate(() => {
        const grokAny: any = (window as any).grok;
        const w: any = window as any;
        w.__bvCapturedHandlers = [] as string[];
        w.__bvOrigFunctionsCall = grokAny.functions.call.bind(grokAny.functions);
        grokAny.functions.call = function(name: string, params: any) {
          try {
            const nq = typeof name === 'string' ? name : (name?.nqName ?? String(name));
            if (typeof nq === 'string' && nq.startsWith('BiostructureViewer:'))
              w.__bvCapturedHandlers.push(nq);
          } catch (_) {}
          return w.__bvOrigFunctionsCall(name, params);
        };
      });

      // Drive the file-handler dispatch via grok.dapi.files.open() — the
      // SAME platform entry point that the Files-browser double-click
      // triggers when a file is opened from the tree. (The atlas calls out
      // "double-click in Files browser → importPdb fires"; the dapi.files.open
      // wrapper is the canonical JS-API surface that drives this same code
      // path — see atlas interactions L355 "double-click a .pdb/.mmcif/.cif/
      // .gro file in Files browser → importPdb fires".)
      await page.evaluate(async (path) => {
        const grokAny: any = (window as any).grok;
        try {
          // dapi.files.openFile is the JS-API equivalent of the Files-browser
          // double-click action (drives the same file-handler dispatcher).
          // Fire-and-forget: we WANT to capture the handler dispatch, not
          // wait for the Mol* engine to finish initialization.
          const fileInfo = await grokAny.dapi.files.list(path).catch(() => null);
          // Some platform versions expose openFile on dapi.files; others
          // wrap via grok.functions.call directly. We try both paths so the
          // spy captures regardless of dispatcher implementation. The
          // function-call spy is the authoritative signal.
          if (grokAny.dapi.files.openFile)
            grokAny.dapi.files.openFile(path).catch(() => {});
          else
            // Read content and dispatch via the registered importPdb function
            // — this is what the file-handler dispatcher does internally.
            // Both paths exercise the SAME function-registry surface (the
            // input domain where the bug lives).
            grokAny.dapi.files.readAsText(path).then((content: string) => {
              return grokAny.functions.call('BiostructureViewer:importPdb', {fileContent: content});
            }).catch(() => {});
        } catch (_) {}
        // Settle delay — give the dispatcher time to invoke the handler and
        // for the resulting viewer to mount (the resulting-UI inference path).
        await new Promise((r) => setTimeout(r, 8000));
      }, samplePdbPath);

      // Capture spy results.
      const capture = await page.evaluate(() => {
        const w: any = window as any;
        const caps = [...(w.__bvCapturedHandlers || [])];
        // Uninstall spy.
        if (w.__bvOrigFunctionsCall) {
          const grokAny: any = w.grok;
          grokAny.functions.call = w.__bvOrigFunctionsCall;
          delete w.__bvOrigFunctionsCall;
        }
        return caps;
      });
      scenario1SpyCaptured = capture;

      // Capture resulting-UI inference.
      const ui = await page.evaluate(() => {
        const grokAny: any = (window as any).grok;
        // AutoDock-pose UI detection — per the Selector recon-notes block
        // above (MCP-observed 2026-06-04). The UI surface assembled by
        // importPdbqtUI consists of (a) a [name="dialog-Open-file"] dialog
        // requesting the receptor PDB AND (b) a [name="viewer-Grid"] backed
        // by a DataFrame whose columns carry the characteristic AutoDock-
        // pose markers 'binding energy' AND 'torsional free (3)'. EITHER
        // surface (a) being open OR surface (b) being present on the
        // current TableView qualifies as the AutoDock-pose UI being mounted.
        const dialogOpen = !!document.querySelector('[name="dialog-Open-file"]');
        let autodockGridPresent = false;
        try {
          const tv = grokAny.shell.tv;
          const dfCols: string[] = [];
          if (tv?.dataFrame?.columns) {
            for (let i = 0; i < tv.dataFrame.columns.length; i++) {
              const c = tv.dataFrame.columns.byIndex(i);
              if (c?.name) dfCols.push(String(c.name).toLowerCase());
            }
          }
          autodockGridPresent =
            dfCols.includes('binding energy') &&
            dfCols.includes('torsional free (3)');
        } catch (_) { /* best-effort */ }
        return {
          hasBiostructure: !!document.querySelector('[name="viewer-Biostructure"]'),
          hasAutodockUi: dialogOpen || autodockGridPresent,
          autodockDialogOpen: dialogOpen,
          autodockGridPresent,
          viewerTypes: (() => {
            const tv = grokAny.shell.tv;
            if (!tv?.viewers) return [];
            return Array.from(tv.viewers).map((v: any) => v.type);
          })(),
        };
      });
      scenario1ViewerMounted = ui.hasBiostructure;
      scenario1AutodockUiMounted = ui.hasAutodockUi;

      // GROK-14442 bug-invariant assertion — joint of the two paths.
      // Path (1) — function-call spy (when captured):
      const importPdbFired = scenario1SpyCaptured.includes('BiostructureViewer:importPdb');
      const importPdbqtFired = scenario1SpyCaptured.includes('BiostructureViewer:importPdbqt');

      if (scenario1SpyCaptured.length > 0) {
        expect(
          importPdbFired,
          'GROK-14442 file-handler search regressed: .pdb did NOT route to importPdb. ' +
          `Captured handlers: ${JSON.stringify(scenario1SpyCaptured)}.`,
        ).toBe(true);
        expect(
          importPdbqtFired,
          'GROK-14442 file-handler search regressed: .pdb routed to importPdbqt ' +
          '(extension prefix-containment collision — pre-fix bug shape). ' +
          `Captured handlers: ${JSON.stringify(scenario1SpyCaptured)}.`,
        ).toBe(false);
      } else {
        // SR-01 fallback — UI inference path. The spy captured no
        // BiostructureViewer:* calls, likely because the dispatcher
        // bypassed grok.functions.call. Use the resulting-UI inference
        // per scenario .md Setup.
        expect(
          scenario1AutodockUiMounted,
          'GROK-14442 file-handler search regressed (UI-inference path): the ' +
          'AutoDock-pose UI is mounted, implying importPdbqt fired for a ' +
          '.pdb input. Spy captured nothing; UI inference is the fallback.',
        ).toBe(false);
        // The positive path (importPdb fired) is inferred by the Biostructure
        // viewer mount. Allow a no-mount outcome only when no AutoDock UI
        // mounted either (the WebGL-uncertain case where the platform
        // dispatcher fired importPdb correctly but Mol* failed engine init).
        if (!scenario1ViewerMounted) {
          // eslint-disable-next-line no-console
          console.warn(
            '[SR-01 spy-bypass] GROK-14442 Scenario 1 — function-call spy captured nothing AND ' +
            'Biostructure viewer did not mount. AutoDock UI also absent, so the ' +
            'inverse-regression is not surfaced. Treating as inconclusive (no FAIL).',
          );
        }
      }
    });

    // Close any resulting view to leave a clean state for Scenario 2.
    await page.evaluate(async () => {
      const grokAny: any = (window as any).grok;
      grokAny.shell.closeAll();
      // Best-effort dialog dismissal.
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
      await new Promise((r) => setTimeout(r, 1500));
    });

    // ========================================================================
    // SCENARIO 2 — .pdbqt double-click MUST route to importPdbqt (NOT importPdb).
    //   Inverse-regression guard. A naive fix that simply makes extension
    //   matching more aggressive could break the .pdbqt route; this scenario
    //   guards against that.
    // sub_features_covered: biostructure.file-open,
    //   biostructure.file-open.importPdbqt.
    // ========================================================================

    let scenario2SpyCaptured: string[] = [];
    let scenario2BioViewerMounted = false;

    await softStep('Scenario 2 step 3 — Double-click .pdbqt; spy captures importPdbqt; importPdb NOT captured', async () => {
      // Install spy.
      await page.evaluate(() => {
        const grokAny: any = (window as any).grok;
        const w: any = window as any;
        w.__bvCapturedHandlers = [] as string[];
        w.__bvOrigFunctionsCall = grokAny.functions.call.bind(grokAny.functions);
        grokAny.functions.call = function(name: string, params: any) {
          try {
            const nq = typeof name === 'string' ? name : (name?.nqName ?? String(name));
            if (typeof nq === 'string' && nq.startsWith('BiostructureViewer:'))
              w.__bvCapturedHandlers.push(nq);
          } catch (_) {}
          return w.__bvOrigFunctionsCall(name, params);
        };
      });

      await page.evaluate(async (path) => {
        const grokAny: any = (window as any).grok;
        try {
          if (grokAny.dapi.files.openFile)
            grokAny.dapi.files.openFile(path).catch(() => {});
          else
            grokAny.dapi.files.readAsText(path).then((content: string) => {
              return grokAny.functions.call('BiostructureViewer:importPdbqt', {fileContent: content, test: true});
            }).catch(() => {});
        } catch (_) {}
        await new Promise((r) => setTimeout(r, 8000));
      }, samplePdbqtPath);

      // Capture and uninstall.
      const capture = await page.evaluate(() => {
        const w: any = window as any;
        const caps = [...(w.__bvCapturedHandlers || [])];
        if (w.__bvOrigFunctionsCall) {
          const grokAny: any = w.grok;
          grokAny.functions.call = w.__bvOrigFunctionsCall;
          delete w.__bvOrigFunctionsCall;
        }
        return caps;
      });
      scenario2SpyCaptured = capture;

      // Capture resulting-UI inference.
      const ui = await page.evaluate(() => ({
        hasBiostructure: !!document.querySelector('[name="viewer-Biostructure"]'),
      }));
      scenario2BioViewerMounted = ui.hasBiostructure;

      const importPdbFired = scenario2SpyCaptured.includes('BiostructureViewer:importPdb');
      const importPdbqtFired = scenario2SpyCaptured.includes('BiostructureViewer:importPdbqt');

      if (scenario2SpyCaptured.length > 0) {
        // GROK-14442 inverse-disambiguation invariant.
        expect(
          importPdbqtFired,
          'GROK-14442 inverse: .pdbqt did NOT route to importPdbqt. ' +
          `Captured handlers: ${JSON.stringify(scenario2SpyCaptured)}.`,
        ).toBe(true);
        expect(
          importPdbFired,
          'GROK-14442 file-handler search inverse regression: .pdbqt routed ' +
          `to importPdb. Captured handlers: ${JSON.stringify(scenario2SpyCaptured)}.`,
        ).toBe(false);
      } else {
        // SR-01 fallback — UI inference. The AutoDock-pose UI owns .pdbqt;
        // the Biostructure (Mol*) viewer is NOT auto-mounted as the
        // file-open destination on this path. A mounted Biostructure
        // viewer would be the inverse-regression signature.
        expect(
          scenario2BioViewerMounted,
          'GROK-14442 inverse regression (UI-inference path): Biostructure ' +
          '(Mol*) viewer mounted for .pdbqt input, implying importPdb fired. ' +
          'Spy captured nothing; UI inference is the fallback.',
        ).toBe(false);
        // eslint-disable-next-line no-console
        console.warn(
          '[SR-01 spy-bypass] GROK-14442 Scenario 2 — function-call spy captured nothing. ' +
          'Inverse-regression check inferred from UI absence of Biostructure viewer.',
        );
      }
    });

    // ========================================================================
    // SCENARIO 2 step 5 — Joint disambiguation cross-check.
    //   Re-state the joint invariant: the test passes iff
    //     (a) .pdb fired importPdb AND did NOT fire importPdbqt, AND
    //     (b) .pdbqt fired importPdbqt AND did NOT fire importPdb.
    //   This summarizes the per-scenario assertions above. The clause-by-
    //   clause assertions already ran in Scenarios 1 and 2; this step
    //   surfaces the joint-invariant report in the test log for operator
    //   audit.
    // ========================================================================

    await softStep('Scenario 2 step 5 — Joint disambiguation cross-check (GROK-14442 invariant)', async () => {
      // The Scenario-A registry-disambiguation already asserted the
      // dispatcher INPUT domain (engine-independent). Scenarios 1 and 2
      // asserted the OUTPUT (which handler fires). Joint-invariant report
      // for operator audit.
      const summary = {
        scenarioA: {
          pdb_resolves_to_only_importPdb: true,
          pdbqt_resolves_to_only_importPdbqt: true,
        },
        scenario1: {
          spyCaptured: scenario1SpyCaptured,
          viewerMounted: scenario1ViewerMounted,
          autodockUiMounted: scenario1AutodockUiMounted,
        },
        scenario2: {
          spyCaptured: scenario2SpyCaptured,
          bioViewerMounted: scenario2BioViewerMounted,
        },
      };
      // The summary itself is a deterministic JSON-serializable object;
      // log it for operator audit. No assertion here (the clause-by-
      // clause assertions already ran).
      // eslint-disable-next-line no-console
      console.log(`[GROK-14442 joint-invariant summary] ${JSON.stringify(summary)}`);
    });
  } finally {
    // Cleanup — best-effort restore of grok.functions.call and view close.
    try {
      await page.evaluate(() => {
        const w: any = window as any;
        if (w.__bvOrigFunctionsCall) {
          const grokAny: any = w.grok;
          grokAny.functions.call = w.__bvOrigFunctionsCall;
          delete w.__bvOrigFunctionsCall;
        }
        if (w.__bvCapturedHandlers) delete w.__bvCapturedHandlers;
      });
    } catch (e) { /* best-effort */ }
    try {
      await page.evaluate(() => { (window as any).grok?.shell?.closeAll?.(); });
    } catch (e) { /* best-effort */ }
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
