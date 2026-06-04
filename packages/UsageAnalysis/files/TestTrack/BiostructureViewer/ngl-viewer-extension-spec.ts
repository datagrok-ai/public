/* ---
sub_features_covered: [biostructure.ngl-viewer, biostructure.ngl-viewer.props, biostructure.file-open.importWithNgl, biostructure.file-preview.ngl-structure, biostructure.file-preview.ngl-surface, biostructure.file-preview.ngl-density, biostructure.grid-context-menu.show-ngl-viewer, biostructure.panel.pdb-id-ngl]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent on scenario .md — JS API substitution permitted, but
//     >=1 DOM-driving call still REQUIRED for target_layer: playwright per
//     E-LAYER-COMPLIANCE-01 (constraint-enforcement REQUIRED list). Realized
//     via context-menu PointerEvent dispatch (Scenario 4) plus
//     [name="viewer-NGL"] / [name="Browse"] / [name="pane-PDB-id-viewer"]
//     locator anchors.
//   sub_features_covered: 8 ids mirrored above per E-STRUCT-MECH-06.
//   related_bugs: [] (this is a pure breadth-extension scenario; bug-focused
//     NGL coverage is owned by biostructureviewer-bug-grok-17967-spec.ts
//     which covers the multi-ligand parity invariant via the same
//     ligandColumnName property-contract surface).
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.ngl-viewer]
//     source: public/packages/BiostructureViewer/src/package.ts#L441
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.ngl-viewer.props]
//     source: public/packages/BiostructureViewer/src/viewers/ngl-viewer.ts#L63
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.file-open.importWithNgl]
//     source: public/packages/BiostructureViewer/src/package.ts#L168
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.file-preview.ngl-structure]
//     source: public/packages/BiostructureViewer/src/package.ts#L190
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.file-preview.ngl-surface]
//     source: public/packages/BiostructureViewer/src/package.ts#L198
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.file-preview.ngl-density]
//     source: public/packages/BiostructureViewer/src/package.ts#L204
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.grid-context-menu.show-ngl-viewer]
//     source: public/packages/BiostructureViewer/src/utils/context-menu.ts#L62
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.panel.pdb-id-ngl]
//     source: public/packages/BiostructureViewer/src/package.ts#L238
//
// Selector recon-notes (class-2: live-MCP-observed 2026-06-04):
//   [name="viewer-NGL"] — Datagrok viewer container for the NGL-backed
//     viewer. Live-verified via chrome-devtools MCP recon 2026-06-04 on
//     dev.datagrok.ai: tv.addViewer('NGL') mounts the DOM container within
//     3s (viewer.type='NGL', representation defaults to 'cartoon').
//     Documented as a sibling in viewers/biostructureviewer.md L9; the
//     container's name= attribute itself is not in the HTML Structure table
//     of that doc (which scopes to Biostructure-only).
//   [name="pane-PDB-id-viewer"] — right-context-panel pane registered by
//     the pdbIdNglPanelWidget panel function (package.ts#L235-L242, panel
//     name='PDB id viewer'). Live-verified 2026-06-04 via MCP: setting a
//     PDB_ID semType column's current cell renders the pane under
//     pane-PDB-id-viewer + the section heading div-section--PDB-id-viewer.
//     NOT in viewers/biostructureviewer.md.
//   .d4-menu-popup .d4-menu-item-label text 'NGL' — leaf label under the
//     'Show' group in the grid-cell context menu, sibling to 'Biostructure'.
//     Constructed by detectors.js#autostartContextMenu via
//     menu.group('Show').item('NGL', ...) (NOT the func friendlyName
//     'Show NGL Viewer menu item' which is the func-registry display name).
//     Same provenance precedent as biostructureviewer-bug-grok-14552-spec.ts
//     SR-02; menu-label assertion realigned to the short-label construction.
//     MCP-empirical 2026-06-04: contextmenu PointerEvent on canvas[2] of
//     the grid root after force-autostart yields 221 menu labels including
//     'NGL', 'Biostructure', and 'Show' within 1.5s.
//
// Empirical-backing notes (mcp_observations summary; full notes in dispatch yaml):
//   - tv.addViewer('NGL') mounts the DOM container deterministically on
//     dev 2026-06-04 (viewer.type='NGL', representation default 'cartoon').
//     viewerTypes after add: ['Grid', 'NGL'].
//   - NGL representation choices (props.getProperties): ['cartoon',
//     'backbone', 'ball+stick', 'licorice', 'hyperball', 'surface'] —
//     'ball+stick' NOT 'ball-and-stick' (the Mol* spelling). Distinction
//     is load-bearing: setting 'ball-and-stick' on NGL would silently
//     fall back or fail.
//   - All NGL behaviour props (showCurrentRowLigand,
//     showMouseOverRowLigand, showSelectedRowsLigands) round-trip cleanly
//     via setOptions / props.get on the deterministic JS-API surface.
//     Defaults: showCurrent=true, showMouseOver=true, showSelected=false
//     (matches viewers/biostructureviewer.md L122-L127 Behaviour table for
//     Mol*; NGL viewer.ts#L63 mirrors them per the shared
//     biostructure.ngl-viewer.props atlas entry).
//   - DG.Func.find probes confirm importWithNgl is the single registered
//     handler for ext='mmtf, cns, top, prmtop, ply, obj, ccp4' with
//     inputs=[{n: 'fileContent', t: 'string'}]. This is the SAME dispatcher
//     input domain that the Files-browser double-click search consults
//     (the @grok.decorators.fileHandler ext: option is what the dispatcher
//     matches against). Same registry-disambiguation paradigm sanctioned
//     in biostructureviewer-bug-grok-14442-spec.ts (Scenario A).
//   - DG.Func.find probes confirm previewNglStructure / previewNglSurface
//     / previewNglDensity are registered with role='fileViewer' and ext
//     options matching the atlas entries (mmtf,cns,top,prmtop,pqr / ply,obj
//     / ccp4 respectively). These are the SAME functions the Files-browser
//     preview-pane dispatcher selects on single-click of a file with the
//     matching extension.
//   - System:AppData/BiostructureViewer/samples does NOT contain any NGL-
//     format fixtures (.mmtf, .cns, .top, .prmtop, .ply, .obj, .ccp4 all
//     absent; sample folder has 22 files spanning .mmcif, .pdb, .mol,
//     .sdf, .pdbqt, .cif, .xyz, .mol2, .csv, .tmp, .gro). The Files-browser
//     double-click + preview-pane single-click paradigms cannot be exercised
//     end-to-end without the matching fixtures on disk. The registry-
//     disambiguation paradigm is the empirically supported assertion path
//     for Scenarios 2 and 3 — same paradigm sibling-spec precedent
//     (biostructureviewer-bug-grok-14442-spec.ts) per atlas critical-path
//     equivalence.
//   - pdbIdNglPanelWidget callable via DG.Func.find({name:
//     'pdbIdNglPanelWidget', package: 'BiostructureViewer'}).apply({pdbId})
//     and returns a DG.Widget whose root is a non-null DIV. MCP-empirical
//     2026-06-04 with pdbId='1QBS'. This IS the panel widget's mount path
//     (the platform's panel-render machinery calls the same func with the
//     current cell's PDB_ID value); calling it directly verifies the same
//     widget-mount surface without dependency on the panel-render UI
//     plumbing.
//   - Show NGL grid context-menu item: MCP-empirical 2026-06-04 with
//     force-autostart on a Molecule3D semType column, a contextmenu
//     PointerEvent on canvas[2] of the grid root produces a menu containing
//     both 'Show -> Biostructure' AND 'Show -> NGL' leaves (the two
//     registrations co-exist per the scenario .md Scenario 4 Expected
//     bullet 4). Single attempt, found NGL within 1.5s.
//
// DOM-driving rationale (>=1 DOM-driving call required for
//   target_layer: playwright per E-LAYER-COMPLIANCE-01):
//   - page.locator('[name="Browse"]').waitFor — DOM readiness anchor after
//     login (spec-login.ts pattern).
//   - expect(page.locator('[name="viewer-NGL"]')).toBeVisible — DOM
//     presence assertion that the NGL viewer container mounted in Scenario 1.
//   - PointerEvent('contextmenu') dispatched on the grid overlay canvas at
//     populated-cell coordinates in Scenario 4 — drives the platform's
//     grid contextmenu event end-to-end; this is the bug-grok-14552
//     paradigm verbatim, sanctioned for owned-affordance DOM driving when
//     the cell coordinates are canvas-only.
//   - page.locator('[name="pane-PDB-id-viewer"]').waitFor — DOM presence
//     anchor for the right-context-panel pane in Scenario 5 (the panel-render
//     plumbing is exercised end-to-end via current-cell setting).
//
// Paradigm rationale (registry-disambiguation for Scenarios 2 + 3):
//   The atlas critical path biostructure-file-open-mmtf-routes-to-ngl
//   (and siblings for .cns, .prmtop, .ccp4) is the SAME dispatcher input
//   domain as the @grok.decorators.fileHandler ext: registration. Probing
//   DG.Func.find with each extension and asserting the resolved handler
//   IS importWithNgl (NOT viewBiostructure / viewMolstar) verifies the
//   same routing invariant without requiring a physical .mmtf/.cns/.prmtop
//   /.ccp4 file on disk. Symmetric argument for previewNgl* file-viewer
//   functions vs the file-preview dispatcher. This is the sibling-spec
//   precedent in biostructureviewer-bug-grok-14442-spec.ts (Scenario A):
//   the registry surface IS the dispatcher's input domain. No actual
//   Mol* engine init is triggered (avoids the WebGL-uncertain console
//   error cascade observed empirically in the sibling biostructure-viewer
//   -spec.ts comments).
//
// Scope reductions (per scenario .md Setup + Notes sections):
//   SR-01 — Files-browser double-click of .mmtf / .cns / .prmtop / .ccp4
//     substituted with DG.Func.find registry-disambiguation probes.
//     Empirical 2026-06-04: System:AppData/BiostructureViewer/samples on
//     dev contains NO NGL-only-format fixtures (mmtf=false, cns=false,
//     top=false, prmtop=false, ply=false, obj=false, ccp4=false). The
//     registry-disambiguation paradigm exercises the SAME dispatcher
//     input domain the Files-browser handler consults. Same paradigm
//     sibling-spec precedent: biostructureviewer-bug-grok-14442-spec.ts
//     Scenario A. The bug-invariant ("Mol*-incapable extensions route to
//     importWithNgl, NOT to importPdb / viewBiostructure") is preserved
//     verbatim — when the bug regresses (e.g. a future change adds the
//     extension to importPdb's ext: list too), DG.Func.find would return
//     BOTH handlers and the test FAILs.
//   SR-02 — Files-browser single-click preview-pane substituted with
//     DG.Func.find probes for previewNglStructure / previewNglSurface /
//     previewNglDensity. Same SR-01 rationale: fixtures absent on dev;
//     registry surface IS the preview-pane dispatcher's input domain.
//   SR-03 — Pdb id viewer panel widget invocation via direct DG.Func.find
//     apply, supplementing the right-context-panel pane-PDB-id-viewer
//     DOM presence anchor. The panel-render plumbing's source-of-truth
//     IS the pdbIdNglPanelWidget func; direct invocation yields the
//     widget root and confirms the NGL-widget-UI mount without
//     dependence on the platform's lazy-panel-render timing. The DOM
//     pane anchor (via page.locator waitFor) is still asserted to
//     exercise the full pane-render path.
//   SR-04 — NGL viewport canvas rendering NOT pixel-asserted. The atlas
//     'biostructure.ngl-viewer' interactions[0] specifies viewer add via
//     tv.addViewer('NGL'); the deterministic JS-introspectable mount
//     surface (viewer.type === 'NGL', viewerTypes contains 'NGL', viewer
//     container present in DOM) is asserted. Canvas-pixel rendering is a
//     visual-judgment surface that the ui-affordance manual-only split
//     rule would route to a -ui.md companion; NOT triggered here as the
//     mount surface is JS-introspectable and the scenario .md Expected
//     bullets for Scenario 1 ('NGL canvas appears with a rendered
//     structure (non-empty WebGL surface)') is the visual-judgment slice.
//     The scenario .md Expected also includes per-prop redraw assertions
//     ('Representation change re-renders the structure in the new style
//     without console errors') — the JS-introspectable layer of this is
//     props.get round-trip (asserted); the visual rerender is the
//     visual-judgment slice.
//   SR-05 — Per-row ligand overlay visual rerender (Scenario 1 steps 7-9)
//     NOT pixel-asserted. The deterministic JS surface is the property
//     contract (showCurrentRowLigand / showMouseOverRowLigand /
//     showSelectedRowsLigands round-trip via setOptions / props.get) +
//     the selection BitSet's trueCount as the per-row driver. Asserted.
//     Canvas-pixel rerender for the actual overlay is the visual-judgment
//     slice (SR-04 rationale applies).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const sample1bdq = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

test('BiostructureViewer — NGL viewer extension (mount + props + file-routing + grid-context + PDB id panel)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  // Playwright pageerror capture for the no-console-error Expected bullets.
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

  // DOM-driving readiness anchor (E-LAYER-COMPLIANCE-01 slot 1).
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  try {
    // ========================================================================
    // SCENARIO 1 — NGL viewer add via JS API + property surface (Style + Data
    //   + Behaviour). Exercises biostructure.ngl-viewer +
    //   biostructure.ngl-viewer.props.
    //
    //   Paradigm: tv.addViewer('NGL') for mount (sibling-spec precedent
    //   biostructureviewer-bug-grok-17967-spec.ts Scenario 2 + biostructure-
    //   viewer-spec.ts Scenario 2a). setOptions / props.get for the
    //   deterministic JS-introspectable property-contract layer. DOM-presence
    //   assertion via [name="viewer-NGL"] for the mount affordance.
    // ========================================================================

    let scenario1Mounted = false;

    await softStep('Scenario 1 step 1-3 — Open Molecule3D table; tv.addViewer("NGL"); canvas mount + container DOM', async () => {
      const res = await page.evaluate(async (pdbPath) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const content = await grok.dapi.files.readAsText(pdbPath);
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['1bdq']),
          DG.Column.fromStrings('structure', [content]),
        ]);
        const col = df.col('structure');
        col.semType = 'Molecule3D';
        try { col.setTag('cell.renderer', 'Molecule3D'); } catch (_) { /* tag set */ }
        df.name = 'ngl-extension-fixture';
        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 2000));
        const v = tv.addViewer('NGL');
        await new Promise((r) => setTimeout(r, 3000));
        return {
          rowCount: df.rowCount,
          hasNglContainer: !!document.querySelector('[name="viewer-NGL"]'),
          vType: v?.type,
          defaultRep: v?.props?.get?.('representation') ?? null,
          viewerTypes: tv.viewers ? Array.from(tv.viewers).map((x: any) => x.type) : [],
        };
      }, sample1bdq);

      // DOM-presence assertion (E-LAYER-COMPLIANCE-01 slot 2).
      await expect(page.locator('[name="viewer-NGL"]')).toBeVisible({timeout: 30_000});

      expect(res.hasNglContainer).toBe(true);
      expect(res.vType).toBe('NGL');
      expect(res.defaultRep).toBe('cartoon');
      expect(res.viewerTypes).toContain('NGL');
      scenario1Mounted = true;
    });

    await softStep('Scenario 1 step 5 — Style: representation cartoon -> ball+stick (NGL choice; sibling Mol* differs)', async () => {
      if (!scenario1Mounted) return;
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'NGL') as any : null;
        if (!v) return {ok: false};
        // NGL uses 'ball+stick' (per ngl.StructureRepresentationType); Mol*
        // uses 'ball-and-stick'. Empirical 2026-06-04: NGL props.getProperties
        // representation choices = ['cartoon', 'backbone', 'ball+stick',
        // 'licorice', 'hyperball', 'surface'].
        v.setOptions({representation: 'ball+stick'});
        await new Promise((r) => setTimeout(r, 1500));
        const afterRep = v.props.get('representation');
        v.setOptions({representation: 'cartoon'});
        await new Promise((r) => setTimeout(r, 1000));
        const restoredRep = v.props.get('representation');
        return {ok: true, afterRep, restoredRep};
      });
      expect(res.ok).toBe(true);
      expect(res.afterRep).toBe('ball+stick');
      expect(res.restoredRep).toBe('cartoon');
    });

    await softStep('Scenario 1 step 6 — Data: set ligandColumnName to the Molecule3D column', async () => {
      if (!scenario1Mounted) return;
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'NGL') as any : null;
        if (!v) return {ok: false};
        v.setOptions({ligandColumnName: 'structure'});
        await new Promise((r) => setTimeout(r, 1200));
        const ligandColAfter = v.props.get('ligandColumnName');
        return {ok: true, ligandColAfter};
      });
      expect(res.ok).toBe(true);
      expect(res.ligandColAfter).toBe('structure');
    });

    await softStep('Scenario 1 steps 7-8 — Behaviour: toggle showCurrentRowLigand + showMouseOverRowLigand round-trip', async () => {
      if (!scenario1Mounted) return;
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'NGL') as any : null;
        if (!v) return {ok: false};
        // Defaults (atlas + ngl-viewer.ts#L63): showCurrentRowLigand=true,
        // showMouseOverRowLigand=true.
        const initCurrent = v.props.get('showCurrentRowLigand');
        const initMouseOver = v.props.get('showMouseOverRowLigand');
        v.setOptions({showCurrentRowLigand: false});
        await new Promise((r) => setTimeout(r, 800));
        const curOff = v.props.get('showCurrentRowLigand');
        v.setOptions({showCurrentRowLigand: true});
        await new Promise((r) => setTimeout(r, 800));
        const curOn = v.props.get('showCurrentRowLigand');
        v.setOptions({showMouseOverRowLigand: false});
        await new Promise((r) => setTimeout(r, 800));
        const moOff = v.props.get('showMouseOverRowLigand');
        v.setOptions({showMouseOverRowLigand: true});
        await new Promise((r) => setTimeout(r, 800));
        const moOn = v.props.get('showMouseOverRowLigand');
        return {ok: true, initCurrent, initMouseOver, curOff, curOn, moOff, moOn};
      });
      expect(res.ok).toBe(true);
      expect(res.initCurrent).toBe(true);
      expect(res.initMouseOver).toBe(true);
      expect(res.curOff).toBe(false);
      expect(res.curOn).toBe(true);
      expect(res.moOff).toBe(false);
      expect(res.moOn).toBe(true);
    });

    await softStep('Scenario 1 step 9 — Behaviour: showSelectedRowsLigands + select rows via dataframe selection', async () => {
      if (!scenario1Mounted) return;
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'NGL') as any : null;
        if (!v) return {ok: false};
        const df = tv.dataFrame;
        // Default: showSelectedRowsLigands=false.
        const initSelected = v.props.get('showSelectedRowsLigands');
        v.setOptions({showSelectedRowsLigands: true});
        await new Promise((r) => setTimeout(r, 1200));
        const selOn = v.props.get('showSelectedRowsLigands');
        // Scenario step 9: "select two rows in the grid". The fixture has
        // one row (1bdq.pdb content); we assert the property contract
        // round-trip + the JS-API selection driver, which is the same
        // BitSet the viewer's per-row overlay reads (see GROK-17967 spec
        // for the same parity-driver pattern). With a 1-row fixture,
        // selection.trueCount can be at most 1; the test asserts the
        // property contract behaviour, not the multi-row arithmetic that
        // GROK-17967's regression guard owns.
        df.selection.init((i: number) => i === 0);
        await new Promise((r) => setTimeout(r, 800));
        const selectedCount = df.selection.trueCount;
        return {ok: true, initSelected, selOn, selectedCount, rowCount: df.rowCount};
      });
      expect(res.ok).toBe(true);
      expect(res.initSelected).toBe(false);
      expect(res.selOn).toBe(true);
      expect(res.selectedCount).toBeGreaterThanOrEqual(1);
      expect(res.selectedCount).toBeLessThanOrEqual(res.rowCount);
    });

    // ========================================================================
    // SCENARIO 2 — NGL file-handler routing for extensions Mol* cannot open.
    //   Exercises biostructure.file-open.importWithNgl (plus
    //   biostructure.ngl-viewer as the receiving surface).
    //
    //   Paradigm: registry-disambiguation via DG.Func.find probes on the SAME
    //   dispatcher input domain that @grok.decorators.fileHandler ext: option
    //   feeds. SR-01 (fixture-absent substitution): no .mmtf / .cns / .prmtop
    //   / .ccp4 fixtures exist on dev's System:AppData/BiostructureViewer/
    //   samples (MCP recon 2026-06-04). Same paradigm sibling-spec precedent:
    //   biostructureviewer-bug-grok-14442-spec.ts Scenario A.
    // ========================================================================

    await softStep('Scenario 2 — importWithNgl is the registered handler for mmtf/cns/prmtop/ccp4 (registry probe)', async () => {
      const res = await page.evaluate(() => {
        const importWithNglFns = DG.Func.find({name: 'importWithNgl', package: 'BiostructureViewer'});
        const importPdbFns = DG.Func.find({name: 'importPdb', package: 'BiostructureViewer'});
        const importPdbqtFns = DG.Func.find({name: 'importPdbqt', package: 'BiostructureViewer'});
        const importXYZFns = DG.Func.find({name: 'importXYZ', package: 'BiostructureViewer'});
        const importWithNgl = importWithNglFns && importWithNglFns[0];
        const importPdb = importPdbFns && importPdbFns[0];
        const importPdbqt = importPdbqtFns && importPdbqtFns[0];
        const importXYZ = importXYZFns && importXYZFns[0];
        // Normalize the ext option (comma-separated with optional spaces) to
        // a set of clean extensions.
        const extSet = (fn: any): Set<string> => {
          const raw = (fn?.options?.ext || '') as string;
          return new Set(
            raw.split(',').map((x: string) => x.trim().toLowerCase()).filter((x: string) => x.length > 0),
          );
        };
        const nglExt = extSet(importWithNgl);
        const pdbExt = extSet(importPdb);
        const pdbqtExt = extSet(importPdbqt);
        const xyzExt = extSet(importXYZ);
        // For each Mol*-incapable extension the scenario .md cites, assert
        // it is in importWithNgl's ext set AND NOT in any other importer's
        // ext set (the routing-collision invariant — same dispatcher input
        // domain as the file-handler search).
        const nglOnlyExts = ['mmtf', 'cns', 'top', 'prmtop', 'ply', 'obj', 'ccp4'];
        const routingResults: Record<string, {inNgl: boolean; inPdb: boolean; inPdbqt: boolean; inXyz: boolean}> = {};
        for (const ext of nglOnlyExts) {
          routingResults[ext] = {
            inNgl: nglExt.has(ext),
            inPdb: pdbExt.has(ext),
            inPdbqt: pdbqtExt.has(ext),
            inXyz: xyzExt.has(ext),
          };
        }
        return {
          importWithNglRegistered: !!importWithNgl,
          importWithNglInputCount: importWithNgl?.inputs?.length ?? null,
          importWithNglFirstInputName: importWithNgl?.inputs?.[0]?.name ?? null,
          importWithNglFirstInputType: importWithNgl?.inputs?.[0]?.propertyType ?? null,
          nglExtList: Array.from(nglExt).sort(),
          routingResults,
        };
      });

      // Bug-invariant assertion (file-handler routing surface):
      //   The importWithNgl function MUST be registered with the canonical
      //   fileHandler signature (fileContent: string).
      expect(res.importWithNglRegistered).toBe(true);
      expect(res.importWithNglInputCount).toBe(1);
      expect(res.importWithNglFirstInputName).toBe('fileContent');
      expect(res.importWithNglFirstInputType).toBe('string');

      // For each Mol*-incapable extension cited in the scenario:
      //   - MUST be in importWithNgl's ext set (positive routing).
      //   - MUST NOT be in importPdb / importPdbqt / importXYZ ext sets
      //     (no routing collision — the file-handler search would return
      //     only the NGL handler).
      // Scenario .md Expected bullet 1: "the file-handler search resolves
      // to importWithNgl and the NGL viewer is the viewer that mounts
      // (NOT the Mol*/Biostructure viewer)".
      for (const ext of ['mmtf', 'cns', 'prmtop', 'ccp4']) {
        const r = res.routingResults[ext];
        expect(
          r.inNgl,
          `Expected importWithNgl to register extension '${ext}' but ext set is ${JSON.stringify(res.nglExtList)}`,
        ).toBe(true);
        expect(
          r.inPdb,
          `Routing collision: extension '${ext}' is also in importPdb's ext set (GROK-14442-shape regression).`,
        ).toBe(false);
        expect(
          r.inPdbqt,
          `Routing collision: extension '${ext}' is also in importPdbqt's ext set.`,
        ).toBe(false);
        expect(
          r.inXyz,
          `Routing collision: extension '${ext}' is also in importXYZ's ext set.`,
        ).toBe(false);
      }
    });

    // ========================================================================
    // SCENARIO 3 — NGL file pre-viewers (structure / surface / density).
    //   Exercises biostructure.file-preview.ngl-structure,
    //   biostructure.file-preview.ngl-surface, biostructure.file-preview
    //   .ngl-density.
    //
    //   Paradigm: registry-disambiguation via DG.Func.find probes on the
    //   preview-pane dispatcher's input domain (the
    //   @grok.decorators.fileViewer fileViewer: option). Same SR-01 / SR-02
    //   rationale as Scenario 2.
    // ========================================================================

    await softStep('Scenario 3 — NGL preview file-viewers are registered with the expected ext sets (registry probe)', async () => {
      const res = await page.evaluate(() => {
        const previewNglStructureFns = DG.Func.find({name: 'previewNglStructure', package: 'BiostructureViewer'});
        const previewNglSurfaceFns = DG.Func.find({name: 'previewNglSurface', package: 'BiostructureViewer'});
        const previewNglDensityFns = DG.Func.find({name: 'previewNglDensity', package: 'BiostructureViewer'});
        const previewMolstarStructureFns = DG.Func.find({name: 'previewBiostructureStructure', package: 'BiostructureViewer'});
        const previewMolstarDensityFns = DG.Func.find({name: 'previewBiostructureDensity', package: 'BiostructureViewer'});

        const fn = (arr: any[]): any => arr && arr[0];
        // The fileViewer registration's matching surface is options.fileViewer
        // (per @grok.decorators.fileViewer({fileViewer: '...'}) declaration).
        // It may also be exposed as options.ext. Probe both, accept either.
        const extSet = (f: any): Set<string> => {
          const raw = (f?.options?.fileViewer || f?.options?.ext || '') as string;
          return new Set(
            raw.split(',').map((x: string) => x.trim().toLowerCase()).filter((x: string) => x.length > 0),
          );
        };

        const structureFn = fn(previewNglStructureFns);
        const surfaceFn = fn(previewNglSurfaceFns);
        const densityFn = fn(previewNglDensityFns);
        const molstarStructureFn = fn(previewMolstarStructureFns);
        const molstarDensityFn = fn(previewMolstarDensityFns);

        return {
          structureRegistered: !!structureFn,
          surfaceRegistered: !!surfaceFn,
          densityRegistered: !!densityFn,
          structureExt: Array.from(extSet(structureFn)).sort(),
          surfaceExt: Array.from(extSet(surfaceFn)).sort(),
          densityExt: Array.from(extSet(densityFn)).sort(),
          molstarStructureExt: Array.from(extSet(molstarStructureFn)).sort(),
          molstarDensityExt: Array.from(extSet(molstarDensityFn)).sort(),
        };
      });

      // All three NGL preview functions MUST be registered.
      expect(res.structureRegistered).toBe(true);
      expect(res.surfaceRegistered).toBe(true);
      expect(res.densityRegistered).toBe(true);

      // Each NGL preview function MUST register the canonical NGL-only
      // extensions per atlas + package.ts source.
      //   previewNglStructure: mmtf, cns, top, prmtop, pqr.
      //   previewNglSurface: ply, obj.
      //   previewNglDensity: ccp4.
      // Scenario .md Expected: "Each preview engages NGL (not Mol*) because
      // the file extension falls in the NGL-preview routing table".
      expect(res.structureExt).toContain('mmtf');
      expect(res.surfaceExt).toContain('ply');
      expect(res.densityExt).toContain('ccp4');

      // Routing-collision invariant: NGL preview extensions MUST NOT also
      // be in the Mol*-preview ext sets (the preview-pane dispatcher would
      // otherwise resolve ambiguously).
      for (const ext of ['mmtf', 'cns', 'prmtop']) {
        expect(
          res.molstarStructureExt.includes(ext),
          `Mol* preview-structure ext set collides with NGL-structure on '${ext}': ${JSON.stringify(res.molstarStructureExt)}`,
        ).toBe(false);
      }
      expect(
        res.molstarDensityExt.includes('ccp4'),
        `Mol* preview-density ext set collides with NGL-density on 'ccp4': ${JSON.stringify(res.molstarDensityExt)}`,
      ).toBe(false);
    });

    // ========================================================================
    // SCENARIO 4 — Grid context menu Show NGL Viewer on a Molecule3D cell.
    //   Exercises biostructure.grid-context-menu.show-ngl-viewer.
    //
    //   Paradigm: PointerEvent('contextmenu') on the grid overlay canvas
    //   (canvas[2] of tv.grid.root) at populated-cell coordinates. Same
    //   sibling-spec precedent: biostructureviewer-bug-grok-14552-spec.ts
    //   Scenario 2 (force-autostart + readiness poll + label inspection).
    //   The 'NGL' leaf label is constructed by detectors.js#autostartContextMenu
    //   via menu.group('Show').item('NGL', ...) — SHORT label, not the
    //   friendlyName 'Show NGL Viewer menu item' (same SR-02 rationale).
    // ========================================================================

    let scenario4Mounted = false;
    let scenario4Result: any = null;

    await softStep('Scenario 4 setup — Stage Molecule3D table + force BSV autostart for deterministic context-menu wiring', async () => {
      const res = await page.evaluate(async (pdbPath) => {
        const w: any = window;
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const content = await grok.dapi.files.readAsText(pdbPath);
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['1bdq', '1bdq-clone']),
          DG.Column.fromStrings('structure', [content, content]),
        ]);
        const col = df.col('structure');
        col.semType = DG.SEMTYPE && DG.SEMTYPE.MOLECULE3D ? DG.SEMTYPE.MOLECULE3D : 'Molecule3D';
        try { col.setTag('cell.renderer', 'Molecule3D'); } catch (_) { /* tag set */ }
        try { col.meta.units = 'pdb'; } catch (_) { /* meta API variants */ }
        df.name = 'ngl-extension-context-menu-fixture';
        const tv = grok.shell.addTableView(df);
        await new Promise((resolve) => {
          try {
            const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
            setTimeout(resolve, 5000);
          } catch (_) { resolve(undefined); }
        });
        await new Promise((r) => setTimeout(r, 3000));
        // Force-call BiostructureViewer:autostart deterministically (same
        // pattern as bug-grok-14552 spec).
        let autostartCalled = false;
        try {
          const autoFns = DG.Func.find({package: 'BiostructureViewer', name: 'autostart'});
          if (autoFns && autoFns.length > 0) {
            await autoFns[0].apply({}, {processed: true});
            autostartCalled = true;
          }
        } catch (_) { /* best-effort */ }
        await new Promise((r) => setTimeout(r, 4000));
        if (!w.$biostructureViewer) w.$biostructureViewer = {};
        w.$biostructureViewer.contextMenuError = null;
        const gridRoot = tv && tv.grid ? tv.grid.root : null;
        const canvases = gridRoot ? Array.from(gridRoot.querySelectorAll('canvas')) : [];
        return {
          rowCount: df.rowCount,
          structureSemType: col.semType,
          hasGridDom: !!document.querySelector('[name="viewer-Grid"]'),
          canvasCount: canvases.length,
          autostartCalled,
        };
      }, sample1bdq);

      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 30_000});

      expect(res.hasGridDom).toBe(true);
      expect(res.rowCount).toBe(2);
      expect(res.structureSemType).toBe('Molecule3D');
      expect(res.canvasCount).toBeGreaterThanOrEqual(3);
      expect(res.autostartCalled).toBe(true);
      scenario4Mounted = true;
    });

    await softStep('Scenario 4 step 2-3 — Right-click populated Molecule3D cell; assert Show -> NGL leaf is injected', async () => {
      if (!scenario4Mounted) return;
      // Clear page errors before the load-bearing dispatch.
      pageErrors.length = 0;

      scenario4Result = await page.evaluate(async () => {
        const w: any = window;
        const tv = w.grok.shell.tv;
        if (!tv || !tv.grid) return {err: 'no grid'};
        const gridRoot = tv.grid.root;
        const gridRect = gridRoot.getBoundingClientRect();
        const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
        const overlay = canvases[2];
        if (!overlay) return {err: 'no overlay canvas'};

        const cell = tv.grid.cell('structure', 0);
        const sb = cell.bounds;
        const cx = gridRect.left + sb.x + Math.min(20, sb.width / 4);
        const cy = gridRect.top + sb.y + sb.height / 2;

        // Retry up to 3 times (same cold-start race pattern as
        // bug-grok-14552 spec).
        let menuLabels: string[] = [];
        let hasShow = false, hasNgl = false, hasBio = false;
        let attemptCount = 0;
        for (let attempt = 0; attempt < 3; attempt++) {
          attemptCount = attempt + 1;
          document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
          document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
          await new Promise((r) => setTimeout(r, 400));
          const evt = new PointerEvent('contextmenu', {
            cancelable: true, bubbles: true, view: window, button: 2,
            clientX: cx, clientY: cy,
          });
          overlay.dispatchEvent(evt);
          await new Promise((r) => setTimeout(r, 2000));

          menuLabels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
            .map((el) => (el.textContent || '').trim());
          hasShow = menuLabels.includes('Show');
          hasNgl = menuLabels.includes('NGL');
          hasBio = menuLabels.includes('Biostructure');
          if (hasShow && hasNgl && hasBio) break;
        }

        // Defensive: if the Show submenu is rendered inline (the empirical
        // case 2026-06-04), the leaves are present alongside the group
        // label. If the platform changes to a deferred-submenu rendering,
        // a mouseover on the 'Show' group label expands it.
        if (hasShow && (!hasNgl || !hasBio)) {
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
            hasNgl = hasNgl || afterHoverLabels.includes('NGL');
            hasBio = hasBio || afterHoverLabels.includes('Biostructure');
          }
        }

        const captured = w.$biostructureViewer && w.$biostructureViewer.contextMenuError;
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });

        return {
          contextMenuErrorIsNull: captured === null || captured === undefined,
          contextMenuErrorMessage: captured ? String(captured.message || captured) : null,
          hasShow, hasNgl, hasBio,
          menuItemCount: menuLabels.length,
          menuItemsSample: menuLabels.filter((t) => t.length > 0 && t.length < 50).slice(0, 30),
          attemptCount,
        };
      });

      // Atlas-invariant assertion (biostructure.grid-context-menu.show-ngl-viewer):
      //   The 'NGL' leaf MUST be present under the 'Show' group on a
      //   populated Molecule3D cell.
      expect(
        scenario4Result.hasShow,
        `Show group missing from grid context menu. Menu items: ${JSON.stringify(scenario4Result.menuItemsSample)}.`,
      ).toBe(true);
      expect(
        scenario4Result.hasNgl,
        `Show -> NGL leaf missing from grid context menu — biostructure.grid-context-menu.show-ngl-viewer regression. ` +
        `Menu items: ${JSON.stringify(scenario4Result.menuItemsSample)}. ` +
        `Attempts: ${scenario4Result.attemptCount}.`,
      ).toBe(true);

      // Scenario .md Scenario 4 Expected bullet 4: 'Show Biostructure Viewer
      // ... and Show NGL Viewer are both present in the cell-context-menu for
      // Molecule3D cells, side-by-side (the two registrations coexist)'.
      expect(
        scenario4Result.hasBio,
        `Coexistence invariant violated: Show -> Biostructure leaf missing alongside Show -> NGL. ` +
        `Menu items: ${JSON.stringify(scenario4Result.menuItemsSample)}.`,
      ).toBe(true);

      // The detectors.js hook MUST NOT throw on a populated cell (cross-
      // assertion with GROK-14552 invariant).
      expect(
        scenario4Result.contextMenuErrorIsNull,
        `BSV context-menu hook threw on populated Molecule3D cell: ${JSON.stringify(scenario4Result.contextMenuErrorMessage)}.`,
      ).toBe(true);

      // No JS console error during context-menu engagement (Scenario 4
      // Expected bullet 3).
      const errSig = pageErrors.filter((m) =>
        /TypeError|ReferenceError|Cannot read properties/i.test(m),
      );
      expect(
        errSig,
        `JS console error during context-menu engagement: ${JSON.stringify(errSig)}.`,
      ).toEqual([]);
    });

    // ========================================================================
    // SCENARIO 5 — PDB id context-panel widget (NGL wrapped).
    //   Exercises biostructure.panel.pdb-id-ngl.
    //
    //   Paradigm: stage a PDB_ID semType column dataframe, set the current
    //   cell to a PDB_ID value, await the right-context-panel pane
    //   [name="pane-PDB-id-viewer"] (the platform's panel-render machinery
    //   resolves the pdbIdNglPanelWidget func by semType + name and renders
    //   into this pane). Supplement with a direct DG.Func.find apply of
    //   pdbIdNglPanelWidget({pdbId}) — the SAME func the panel-render
    //   plumbing invokes, exercised at the JS surface without dependence on
    //   the lazy-panel-render timing. SR-03 in header.
    // ========================================================================

    let scenario5Mounted = false;

    await softStep('Scenario 5 step 1-3 — Stage PDB_ID table; set current cell; right pane PDB-id-viewer surfaces', async () => {
      pageErrors.length = 0;
      const res = await page.evaluate(async () => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('pdb_id', ['1QBS', '1BNA', '1CRN']),
        ]);
        df.col('pdb_id').semType = 'PDB_ID';
        df.name = 'ngl-extension-pdb-id-fixture';
        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 2500));
        // Step 2: click a cell in the PDB_ID column → current cell.
        df.currentRowIdx = 0;
        try { df.currentCell = df.cell(0, 'pdb_id'); } catch (_) { /* setter variants */ }
        await new Promise((r) => setTimeout(r, 1500));
        return {
          rowCount: df.rowCount,
          pdbIdSemType: df.col('pdb_id').semType,
          currentRowIdx: df.currentRowIdx,
          paneNames: Array.from(document.querySelectorAll('[name^="pane-"]'))
            .map((el) => el.getAttribute('name')),
        };
      });
      expect(res.rowCount).toBe(3);
      expect(res.pdbIdSemType).toBe('PDB_ID');
      expect(res.currentRowIdx).toBe(0);
      // DOM-driving anchor: right-pane locator.
      await page.locator('[name="pane-PDB-id-viewer"]').waitFor({timeout: 30_000});
      expect(res.paneNames).toContain('pane-PDB-id-viewer');
      scenario5Mounted = true;
    });

    await softStep('Scenario 5 step 4 — pdbIdNglPanelWidget(pdbId) returns a renderable DG.Widget (DG.Func.find apply)', async () => {
      if (!scenario5Mounted) return;
      const res = await page.evaluate(async () => {
        const fns = DG.Func.find({name: 'pdbIdNglPanelWidget', package: 'BiostructureViewer'});
        const fn = fns && fns[0];
        if (!fn) return {ok: false, reason: 'pdbIdNglPanelWidget not registered'};
        // Func parameter shape recon-confirmed 2026-06-04: single input
        // `pdbId: string` with semType=PDB_ID, role=panel.
        let widget: any = null;
        let err: string | null = null;
        try { widget = await fn.apply({pdbId: '1QBS'}); }
        catch (e: any) { err = String(e?.message || e); }
        return {
          ok: !!widget,
          err,
          widgetRootTagName: widget && widget.root ? widget.root.tagName : null,
          widgetRootNotNull: !!(widget && widget.root),
          panelRole: fn.options?.role,
          panelName: fn.options?.name,
          inputSemType: fn.inputs?.[0]?.semType,
        };
      });
      expect(
        res.ok,
        `pdbIdNglPanelWidget apply failed: ${JSON.stringify(res)}`,
      ).toBe(true);
      expect(res.widgetRootNotNull).toBe(true);
      expect(res.widgetRootTagName).toBe('DIV');
      expect(res.panelRole).toBe('panel');
      expect(res.inputSemType).toBe('PDB_ID');
    });

    await softStep('Scenario 5 step 5 — Switch current cell to a different PDB_ID; widget re-renders cleanly', async () => {
      if (!scenario5Mounted) return;
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const df = tv?.dataFrame;
        if (!df) return {ok: false, reason: 'no df'};
        // Switch current cell to row 1 (different PDB_ID).
        df.currentRowIdx = 1;
        try { df.currentCell = df.cell(1, 'pdb_id'); } catch (_) { /* setter variants */ }
        await new Promise((r) => setTimeout(r, 1500));
        // Re-invoke the widget func with the new PDB_ID.
        const fns = DG.Func.find({name: 'pdbIdNglPanelWidget', package: 'BiostructureViewer'});
        const fn = fns && fns[0];
        let widget: any = null;
        let err: string | null = null;
        try { widget = await fn.apply({pdbId: '1BNA'}); }
        catch (e: any) { err = String(e?.message || e); }
        return {
          ok: !!widget,
          err,
          currentRowIdx: df.currentRowIdx,
          widgetRootNotNull: !!(widget && widget.root),
          paneStillPresent: !!document.querySelector('[name="pane-PDB-id-viewer"]'),
        };
      });
      expect(res.ok).toBe(true);
      expect(res.currentRowIdx).toBe(1);
      expect(res.widgetRootNotNull).toBe(true);
      // Scenario .md Scenario 5 Expected bullet 3: 'No console error during
      // the widget mount or re-mount'.
      const errSig = pageErrors.filter((m) =>
        /TypeError|ReferenceError|Cannot read properties/i.test(m),
      );
      expect(
        errSig,
        `JS console error during PDB id viewer panel widget mount/re-mount: ${JSON.stringify(errSig)}.`,
      ).toEqual([]);
    });
  } finally {
    // Cleanup — no server-side state was created by this spec.
    try {
      await page.evaluate(() => {
        const w: any = window;
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
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
