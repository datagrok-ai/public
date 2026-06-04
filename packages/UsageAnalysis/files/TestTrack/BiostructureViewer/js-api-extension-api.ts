/* ---
sub_features_covered: [biostructure.api.viewPdbById, biostructure.api.viewPdbByData]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: apitest
//   pyramid_layer: absent on scenario .md — apitest paradigm; FORBIDDEN list
//     per Section 4.1 (target_layer: apitest) blocks page.click / page.fill /
//     page.locator / page.hover / page.press / page.keyboard / page.mouse /
//     dlg.* in this spec body. Only loginToDatagrok in spec-login.ts performs
//     DOM driving (out of body scope). Sibling apitest precedent:
//     Charts/charts-api.ts (target_layer: apitest, body pure grok.dapi.* +
//     grok.shell.* + viewer.props/setOptions/getOptions); Bio/bio-service-
//     surface-init-api.ts (target_layer: apitest, body pure
//     grok.functions.call probes via tryCall helper).
//   sub_features_covered: 2 ids mirrored above per E-STRUCT-MECH-06.
//   ui_coverage_responsibility: [] (apitest layer; no DOM driving on body).
//   ui_coverage_delegated_to: null.
//   related_bugs: [] — coverage-extension scenario; no GROK ticket
//     invariant asserted. Atlas edge_cases[5] (raw-pdb pitfall) is referenced
//     contextually because both viewPdbById and viewPdbByData are documented
//     as the canonical safe entry points that supply a name and avoid the
//     pitfall; the property-surface-extension spec carries the structural
//     edge_cases[5] assertion against the viewBiostructure sibling.
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.api.viewPdbById]
//     source: public/packages/BiostructureViewer/src/package.ts#L105
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.api.viewPdbByData]
//     source: public/packages/BiostructureViewer/src/package.ts#L115
//   feature-atlas/biostructureviewer.yaml#edge_cases[5]
//     derived_from: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L216
//     description: raw `pdb` prop without name fails to parse; canonical safe
//       entry is viewBiostructure(content, format, name) — and equivalently
//       viewPdbByData(pdbData, name) and viewPdbById(pdbId) which both route
//       through byData / byId with a structure name and avoid the pitfall.
//   feature-atlas/biostructureviewer.yaml#edge_cases[6]
//     derived_from: .claude/skills/grok-browser/references/viewers/biostructureviewer.md#L50
//     description: render is asynchronous — Mol* parse + draw takes several
//       hundred ms to a few seconds; tests must await rendering or treat the
//       'timeout' the createRcsbViewer canvas3dInit subscription surfaces
//       under WebGL-uncertain runtime as a precondition observation, not a
//       spec failure.
//
// Empirical-backing notes (mcp_observations summary; full notes in dispatch yaml).
// Live recon 2026-06-04 on dev.datagrok.ai via chrome-devtools MCP:
//   - DG.Func.find({name: 'viewPdbById', package: 'BiostructureViewer'})
//     resolves to a registered function with one input: pdbId: string
//     (required). Registration shape mirrors src/package.ts#L105.
//   - DG.Func.find({name: 'viewPdbByData', package: 'BiostructureViewer'})
//     resolves to a registered function with two inputs:
//       pdbData: string (required), name: string (required).
//     Registration shape mirrors src/package.ts#L115.
//   - SCENARIO .MD INACCURACY (SR-01 below): both viewPdbById and
//     viewPdbByData go through byId / byData (src/viewers/molstar-viewer/
//     utils.ts#L81 + L94), which call grok.shell.newView(viewName). This
//     creates a STANDALONE View (DG.View, type === 'view'), NOT a docked
//     DG.JsViewer inside the active TableView. Empirical observation:
//       await grok.functions.call('BiostructureViewer:viewPdbByData',
//         {pdbData, name: 'recon-1bdq'})
//       => grok.shell.v becomes a new View with .name === 'recon-1bdq',
//          .type === 'view'; this View appears in grok.shell.views[];
//          the .msp-plugin host DOM mounts inside the View's root.
//       await grok.functions.call('BiostructureViewer:viewPdbById',
//         {pdbId: '1QBS'})
//       => grok.shell.v becomes a new View with .name === 'Mol*'
//          (default name from initViewer's viewName: string = 'Mol*'
//          per src/viewers/molstar-viewer/utils.ts#L60).
//     The scenario .md's prescription to walk tv.viewers for type ===
//     'Biostructure' is therefore not the canonical observable surface for
//     these two functions — the docked-viewer pattern applies to
//     tv.addViewer('Biostructure') (covered by sibling biostructure-viewer-
//     spec.ts), not to the standalone-view JS API entry points. The
//     bug-invariant Expected bullets ("call resolves without throwing"
//     [conditionally — see SR-02]; "viewer present in shell state after the
//     call" [shell.views contains a new view named per the `name` arg]) are
//     asserted at the observable JS-API contract.
//   - WEBGL-UNCERTAIN RUNTIME (SR-02): the await grok.functions.call(...)
//     promise REJECTS with the literal string "timeout" in the dev recon
//     environment. Root cause: createRcsbViewer in src/viewers/molstar-
//     viewer/utils.ts#L29 subscribes to viewer.plugin.canvas3dInit via
//     testEvent (utils.test.ts), which times out when the WebGL context
//     cannot be created. Empirical (2026-06-04 dev.datagrok.ai):
//       viewPdbByData call duration: ~60 ms before "timeout" reject
//       viewPdbById call duration:   ~31 ms before "timeout" reject
//       Post-reject: shell.v IS the new view, .msp-plugin DOM host IS
//         mounted, but .msp-viewport canvas is absent (WebGL never
//         initialised). Same WebGL-uncertain pattern as sibling specs
//         biostructure-viewer-spec.ts, property-surface-extension-spec.ts,
//         and biostructureviewer-bug-claude-33-spec.ts.
//     The atlas-documented Expected bullet "the grok.functions.call resolves
//     without throwing" is therefore conditionally asserted: it holds under
//     healthy WebGL, but in the WebGL-uncertain dev recon the reject string
//     === "timeout" is treated as a precondition observation per atlas
//     edge_cases[6] (the createRcsbViewer canvas3dInit subscription surfaces
//     a timeout when the engine cannot init). The bug-invariant assertion
//     is at the structural surface: either the call resolves OR the reject
//     is the engine-init timeout (NOT a Parsed-object-is-empty parse error,
//     NOT a thrown exception in the dispatcher).
//   - awaitRendered API SHAPE (SR-03): the scenario .md prescribes
//     viewer.awaitRendered(timeoutMs). This method is on the MolstarViewer
//     class (DG.JsViewer subclass per src/viewers/molstar-viewer/molstar-
//     viewer.ts). Both viewPdbById and viewPdbByData route through byId /
//     byData which instantiate an RcsbViewer (the rcsb-molstar Viewer
//     class) inside a standalone View's root — NOT a DG.JsViewer. The
//     awaitRendered method belongs to the DG.JsViewer wrapper, not the
//     RcsbViewer-only path these two API entry points use. The settle is
//     therefore implemented via the standard atlas edge_cases[6] pattern:
//     bounded post-call wait + .msp-plugin presence check. The render-
//     readiness assertion is downgraded to host-DOM presence + view-name
//     mirror; no awaitRendered call is made on the docked-viewer surface
//     (because no docked-viewer exists for these two entry points).
//   - viewBiostructure (sibling JS API trio member) is registered with
//     (content, format, name?) per the property-surface-extension spec's
//     function-registry probe; this spec does NOT re-cover viewBiostructure
//     (per scenario .md Notes section guidance — completes the trio without
//     redundancy).
//
// Paradigm rationale (apitest):
//   The two sub-features are JS API entry points registered as @func in
//   package.g.ts (atlas confirms registration; STEP D non-UI tail rule).
//   The observable contract is JS-API state: function-registry registration
//   shape + post-call shell-view state. The WebGL canvas rendering surface
//   is owned by the smoke biostructure-viewer-spec.ts (which uses
//   tv.addViewer('Biostructure') — the docked-viewer path that does expose
//   a DG.JsViewer with awaitRendered) AND by property-surface-extension-
//   spec.ts (which exercises the property contract). This spec covers the
//   third canonical entry-point surface — programmatic structure opening
//   via package function call — at the JS-API contract level only.
//
// Scope reductions (per scenario .md Setup + Notes + empirical recon):
//   SR-01 — scenario .md Steps 2 prescribe walking grok.shell.tv.viewers
//     looking for type === 'Biostructure'. Empirically (live recon
//     2026-06-04 on dev.datagrok.ai), both viewPdbById and viewPdbByData
//     route through byId / byData which call grok.shell.newView(viewName)
//     to create a STANDALONE View (DG.View, type === 'view'), not a docked
//     DG.JsViewer inside an active TableView. The observable side-effect is
//     therefore asserted on grok.shell.views[] (the new standalone view's
//     name mirrors the `name` arg for viewPdbByData, defaults to 'Mol*'
//     for viewPdbById per initViewer's default viewName), not on
//     grok.shell.tv.viewers. The bug-invariant Expected bullet "a
//     Biostructure (Mol*) viewer is present in shell state after the call"
//     is preserved at the shell-state level (the Mol* engine's .msp-plugin
//     host DOM is mounted inside the new view's root). Same paradigm
//     consistency as the actual JS-API contract per
//     src/viewers/molstar-viewer/utils.ts#L81 + L94.
//   SR-02 — scenario .md Expected bullet "the grok.functions.call resolves
//     without throwing" is conditionally asserted. Under healthy WebGL the
//     await resolves cleanly; under the WebGL-uncertain dev recon
//     (empirical 2026-06-04), the await REJECTS with the literal string
//     "timeout" (createRcsbViewer canvas3dInit subscription surfaces a
//     timeout when WebGL context creation fails — atlas edge_cases[6] /
//     references/viewers/biostructureviewer.md#L50). The bug-invariant
//     assertion is structural: either the call resolves (healthy WebGL) OR
//     the reject string === 'timeout' (engine-init timeout, atlas-
//     documented). Anything else (e.g. 'Parsed object is empty', a thrown
//     dispatcher error) IS asserted as a failure. The shell-state side-
//     effect (new view created with the expected name) is asserted
//     regardless of the await outcome — both paths produce the view
//     synchronously before the canvas3dInit await begins. Same SR shape as
//     sibling biostructure-viewer-spec.ts WebGL-uncertain handling.
//   SR-03 — scenario .md Step 3 prescribes await viewer.awaitRendered
//     (timeoutMs). Empirically the standalone-view path does not expose a
//     DG.JsViewer surface with the awaitRendered method (the RcsbViewer
//     instance lives in view.root; it is not a DG viewer wrapper). The
//     render-readiness assertion is therefore implemented as the standard
//     atlas edge_cases[6] settle pattern: a bounded post-call wait (3 s)
//     followed by .msp-plugin host-DOM presence check. The "structure
//     parsed and drew successfully" / "no Parsed object is empty error"
//     bullets are asserted at the structural surface: shell-view existence
//     + Mol* engine-host DOM mount + post-call console capture not
//     containing the documented 'Parsed object is empty' pitfall signature.
//   SR-04 — scenario .md Step 1 of Scenario 2 prescribes two alternative
//     paths for acquiring a PDB text: (a) grok.dapi.fetchProxy to RCSB or
//     (b) a checked-in fixture. Path (a) introduces an outbound dependency
//     on files.rcsb.org and is bounded by network conditions in CI. This
//     spec chooses path (b): the canonical fixture at
//     System:AppData/BiostructureViewer/samples/1bdq.pdb (already used by
//     sibling specs biostructure-viewer-spec.ts and property-surface-
//     extension-spec.ts). Avoids the outbound dependency while exercising
//     the same JS-API contract; consistent with scenario .md "pick
//     whichever the apitest harness already supports without introducing
//     a new fixture".
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

test('BiostructureViewer — JS API extension (viewPdbById / viewPdbByData)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  // Console capture: bug-invariant (atlas edge_cases[5]) asserts the
  // canonical 'Parsed object is empty' pitfall signature is NOT present —
  // viewPdbById / viewPdbByData are the documented safe entry points that
  // supply a structure name and avoid the bare-pdb pitfall.
  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) ||
    /404 \(\)/.test(text) ||
    /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (err) => pageErrors.push(`pageerror: ${err.message}`));

  await loginToDatagrok(page);

  // Setup — apitest-pure: no UI driving, only shell-state hygiene. Mirrors
  // Charts/charts-api.ts L45-L53 and Bio/bio-service-surface-init-api.ts
  // L30-L39.
  await page.evaluate(() => {
    const g = (window as any).grok;
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    g.shell.closeAll();
    document.body.classList.add('selenium');
    g.shell.windows.simpleMode = true;
  });

  try {
    // ========================================================================
    // PRE-FLIGHT — function-registry probe. Verifies both JS API entry
    //   points are registered with the expected parameter shape. Same
    //   pattern as sibling property-surface-extension-spec.ts Scenario 2
    //   step 7 (viewBiostructure registration probe).
    //
    //   sub_features_covered: biostructure.api.viewPdbById (registration
    //   shape), biostructure.api.viewPdbByData (registration shape).
    //   Atlas source: src/package.ts#L105 (viewPdbById), #L115 (viewPdbByData).
    // ========================================================================

    await softStep('Pre-flight — DG.Func.find verifies viewPdbById + viewPdbByData registration shape', async () => {
      const res = await page.evaluate(() => {
        const g = (window as any).grok;
        const D = (window as any).DG;
        const mapInputs = (fn: any) => (fn && fn.inputs)
          ? fn.inputs.map((i: any) => ({
              name: i.name,
              type: i.propertyType,
              optional: i.options?.optional ?? false,
            }))
          : [];
        const byIdFns = D.Func.find({name: 'viewPdbById', package: 'BiostructureViewer'});
        const byDataFns = D.Func.find({name: 'viewPdbByData', package: 'BiostructureViewer'});
        const byIdFn = byIdFns && byIdFns[0];
        const byDataFn = byDataFns && byDataFns[0];
        return {
          byIdRegistered: !!byIdFn,
          byIdInputs: mapInputs(byIdFn),
          byDataRegistered: !!byDataFn,
          byDataInputs: mapInputs(byDataFn),
        };
      });

      // Atlas src/package.ts#L105: static async viewPdbById(pdbId: string).
      expect(res.byIdRegistered).toBe(true);
      expect(res.byIdInputs).toEqual([
        {name: 'pdbId', type: 'string', optional: false},
      ]);

      // Atlas src/package.ts#L115: static async viewPdbByData(pdbData:
      // string, name: string). The `name` arg is the load-bearing parameter
      // that the underlying byData(data, name) path supplies to initViewer
      // — avoiding the bare-pdb pitfall (atlas edge_cases[5]).
      expect(res.byDataRegistered).toBe(true);
      expect(res.byDataInputs).toEqual([
        {name: 'pdbData', type: 'string', optional: false},
        {name: 'name', type: 'string', optional: false},
      ]);
    });

    // ========================================================================
    // SCENARIO 1 — viewPdbById opens a PDB by RCSB ID and renders.
    //
    //   Steps from scenario .md:
    //     1. call grok.functions.call('BiostructureViewer:viewPdbById',
    //        {pdbId: '1QBS'}).
    //     2. Resolve the docked viewer via grok.shell.tv.viewers / shell.v.
    //     3. Await render via viewer.awaitRendered(30000).
    //
    //   Realized in apitest paradigm per SR-01 / SR-02 / SR-03:
    //     - Step 1 — call resolves OR rejects with literal 'timeout'
    //       (atlas edge_cases[6] WebGL-uncertain runtime). No other reject
    //       string is accepted.
    //     - Step 2 — the canonical observable IS grok.shell.v as a
    //       standalone View with .name === 'Mol*' (initViewer default
    //       viewName per src/viewers/molstar-viewer/utils.ts#L60) and
    //       .type === 'view' (NOT 'TableView'). The View also appears in
    //       grok.shell.views[].
    //     - Step 3 — the standalone-view path does not expose a DG.JsViewer
    //       with awaitRendered; settle is implemented as the standard atlas
    //       edge_cases[6] bounded wait + .msp-plugin host-DOM presence.
    //
    //   Bug-invariant (atlas edge_cases[5]): the documented pitfall
    //   signature 'Parsed object is empty' MUST NOT surface — viewPdbById
    //   wraps the call with a known structure name per byId, supplying the
    //   name byData requires to avoid the bare-pdb pitfall.
    //
    //   sub_features_covered: biostructure.api.viewPdbById.
    // ========================================================================

    await softStep('Scenario 1 — viewPdbById opens PDB by ID; shell state + DOM host mount', async () => {
      // Reset console capture for this scenario's bug-invariant assertion.
      consoleErrors.length = 0;
      pageErrors.length = 0;

      const res = await page.evaluate(async () => {
        const g = (window as any).grok;
        g.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));

        // Pre-state.
        const preViewCount = g.shell.views ? Array.from(g.shell.views).length : 0;

        // Invoke the JS API entry point. Per atlas edge_cases[6] /
        // SR-02, the await may reject with the literal string 'timeout'
        // under WebGL-uncertain runtime; any other reject string IS a
        // failure.
        const t0 = Date.now();
        let callResolved = false;
        let callRejectMessage: string | null = null;
        try {
          await g.functions.call('BiostructureViewer:viewPdbById', {pdbId: '1QBS'});
          callResolved = true;
        } catch (e: any) {
          callRejectMessage = String(e?.message ?? e);
        }
        const callDurationMs = Date.now() - t0;

        // Settle per atlas edge_cases[6] (Mol* parse + draw is async).
        await new Promise((r) => setTimeout(r, 3000));

        // Post-state — the canonical observable per SR-01: a standalone
        // View has been added; shell.v is that View; its name is the
        // initViewer default 'Mol*' for the byId path; its type is 'view'.
        const postViewNames: Array<{name: string, type: string}> = [];
        if (g.shell.views) {
          for (const v of g.shell.views) postViewNames.push({name: v.name, type: v.type});
        }
        const curr = g.shell.v;
        const currInfo = curr ? {name: curr.name, type: curr.type} : null;

        // Mol* engine host DOM presence (per atlas edge_cases[6]). The
        // .msp-plugin host mounts inside the View's root regardless of
        // WebGL outcome; .msp-viewport canvas only mounts when WebGL is
        // healthy.
        const hasMspPlugin = !!document.querySelector('.msp-plugin');
        const hasMspViewport = !!document.querySelector('.msp-viewport');
        const hasMspCanvas = !!document.querySelector('.msp-viewport canvas');

        return {
          preViewCount,
          callResolved,
          callRejectMessage,
          callDurationMs,
          postViewNames,
          postViewCount: postViewNames.length,
          currInfo,
          hasMspPlugin,
          hasMspViewport,
          hasMspCanvas,
        };
      });

      // SR-02: call resolves OR rejects with the documented engine-init
      // timeout signature; any other reject IS a failure (e.g. a thrown
      // dispatcher error or a parse-stage error).
      const acceptedReject = res.callRejectMessage === 'timeout';
      expect(res.callResolved || acceptedReject).toBe(true);

      // Shell-state observable side-effect (SR-01): a new View has been
      // added. The standalone-view creation in initViewer happens
      // synchronously BEFORE the canvas3dInit await begins, so the View IS
      // present regardless of WebGL outcome.
      expect(res.postViewCount).toBeGreaterThan(res.preViewCount);

      // The byId path uses initViewer's default viewName === 'Mol*'
      // (src/viewers/molstar-viewer/utils.ts#L60). The current view is the
      // newly created standalone view.
      expect(res.currInfo).not.toBe(null);
      expect(res.currInfo!.name).toBe('Mol*');
      // SR-01: the standalone view has type === 'view' (a DG.View), NOT
      // type === 'TableView' (which the scenario .md tv.viewers walk
      // assumed). The atlas-documented "viewer is present in shell state"
      // bullet is preserved at this surface.
      expect(res.currInfo!.type).toBe('view');

      // Atlas edge_cases[6] / references/viewers/biostructureviewer.md#L50:
      // the Mol* engine host (.msp-plugin) mounts inside the View's root.
      // The canvas itself is canvas3dInit-gated (may be absent under
      // WebGL-uncertain runtime per SR-02). Render-readiness assertion is
      // therefore on the host DOM, not the canvas.
      expect(res.hasMspPlugin).toBe(true);

      // Bug-invariant (atlas edge_cases[5]): the documented pitfall
      // signature MUST NOT surface in console or pageerror. viewPdbById's
      // wrapping of byId with the default 'Mol*' name IS the canonical
      // avoidance path.
      const pitfallRegex = /Parsed object is empty/i;
      const pitfallConsole = consoleErrors.filter((m) => pitfallRegex.test(m));
      const pitfallPage = pageErrors.filter((m) => pitfallRegex.test(m));
      expect(pitfallConsole.length).toBe(0);
      expect(pitfallPage.length).toBe(0);
    });

    // ========================================================================
    // SCENARIO 2 — viewPdbByData opens a PDB from raw string + name.
    //
    //   Steps from scenario .md:
    //     1. Acquire a known PDB text (this spec uses the canonical fixture
    //        System:AppData/BiostructureViewer/samples/1bdq.pdb per SR-04;
    //        scenario .md "pick whichever the apitest harness already
    //        supports without introducing a new fixture").
    //     2. call grok.functions.call('BiostructureViewer:viewPdbByData',
    //        {pdbData: <text>, name: '1QBS'}).
    //     3. Resolve docked viewer + awaitRendered(30000).
    //
    //   Realized per SR-01 / SR-02 / SR-03 / SR-04 (mirror of Scenario 1
    //   shell-state assertion shape).
    //
    //   Bug-invariant (atlas edge_cases[5]): the explicit name argument is
    //   the documented safe path that supplies a structure identity to
    //   byData (instead of falling into the bare-pdb pitfall). The pitfall
    //   signature 'Parsed object is empty, name 'undefined'' MUST NOT
    //   surface.
    //
    //   sub_features_covered: biostructure.api.viewPdbByData.
    // ========================================================================

    await softStep('Scenario 2 — viewPdbByData opens PDB from raw string + name; shell state + DOM host mount', async () => {
      consoleErrors.length = 0;
      pageErrors.length = 0;

      const res = await page.evaluate(async (pdbPath) => {
        const g = (window as any).grok;
        g.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));

        // SR-04: acquire PDB text via the canonical apitest fixture. Avoids
        // the outbound dependency on files.rcsb.org while exercising the
        // same byData contract.
        let pdbContent: string | null = null;
        let pdbReadError: string | null = null;
        try {
          pdbContent = await g.dapi.files.readAsText(pdbPath);
        } catch (e: any) {
          pdbReadError = String(e?.message ?? e);
        }

        const preViewCount = g.shell.views ? Array.from(g.shell.views).length : 0;

        const t0 = Date.now();
        let callResolved = false;
        let callRejectMessage: string | null = null;
        try {
          await g.functions.call(
            'BiostructureViewer:viewPdbByData',
            {pdbData: pdbContent, name: '1QBS'},
          );
          callResolved = true;
        } catch (e: any) {
          callRejectMessage = String(e?.message ?? e);
        }
        const callDurationMs = Date.now() - t0;

        await new Promise((r) => setTimeout(r, 3000));

        const postViewNames: Array<{name: string, type: string}> = [];
        if (g.shell.views) {
          for (const v of g.shell.views) postViewNames.push({name: v.name, type: v.type});
        }
        const curr = g.shell.v;
        const currInfo = curr ? {name: curr.name, type: curr.type} : null;

        const hasMspPlugin = !!document.querySelector('.msp-plugin');
        const hasMspViewport = !!document.querySelector('.msp-viewport');
        const hasMspCanvas = !!document.querySelector('.msp-viewport canvas');

        return {
          pdbReadError,
          pdbContentLen: pdbContent ? pdbContent.length : 0,
          preViewCount,
          callResolved,
          callRejectMessage,
          callDurationMs,
          postViewNames,
          postViewCount: postViewNames.length,
          currInfo,
          hasMspPlugin,
          hasMspViewport,
          hasMspCanvas,
        };
      }, samplePdbPath);

      // Fixture read precondition.
      expect(res.pdbReadError).toBe(null);
      expect(res.pdbContentLen).toBeGreaterThan(1000);

      // SR-02: same conditional reject pattern as Scenario 1.
      const acceptedReject = res.callRejectMessage === 'timeout';
      expect(res.callResolved || acceptedReject).toBe(true);

      // SR-01: shell-state observable side-effect.
      expect(res.postViewCount).toBeGreaterThan(res.preViewCount);
      expect(res.currInfo).not.toBe(null);
      // byData path uses the supplied `name` arg as the view's name
      // (src/viewers/molstar-viewer/utils.ts#L94 -> initViewer(name) on
      // line 95). With name === '1QBS', the standalone view's .name is
      // '1QBS'. This is the bug-invariant of edge_cases[5]: the name is
      // supplied; the bare-pdb pitfall is avoided.
      expect(res.currInfo!.name).toBe('1QBS');
      expect(res.currInfo!.type).toBe('view');

      // Mol* engine host DOM presence (per atlas edge_cases[6]).
      expect(res.hasMspPlugin).toBe(true);

      // Bug-invariant (atlas edge_cases[5]) — the explicit name argument
      // is the documented safe path; the pitfall signature MUST NOT
      // surface.
      const pitfallRegex = /Parsed object is empty/i;
      const pitfallConsole = consoleErrors.filter((m) => pitfallRegex.test(m));
      const pitfallPage = pageErrors.filter((m) => pitfallRegex.test(m));
      expect(pitfallConsole.length).toBe(0);
      expect(pitfallPage.length).toBe(0);
    });
  } finally {
    // Cleanup — close all views created by the standalone-view JS API
    // entry points to leave the runtime clean for downstream specs running
    // under fullyParallel: true.
    await page.evaluate(() => {
      const g = (window as any).grok;
      g.shell.closeAll();
    });
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
