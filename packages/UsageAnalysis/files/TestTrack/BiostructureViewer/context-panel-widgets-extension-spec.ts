/* ---
sub_features_covered: [biostructure.panel.structure-3d, biostructure.panel.pdb-file-info, biostructure.panel.pdb-info, biostructure.panel.prolif, biostructure.panel.link-molecule-column]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent on scenario .md — JS API substitution permitted, but
//     >=1 DOM-driving call still REQUIRED for target_layer: playwright per
//     E-LAYER-COMPLIANCE-01 (constraint-enforcement REQUIRED list). Realized
//     via [name="Browse"] waitFor + accordion-pane header clicks for the
//     panel surfaces (Scenarios 1, 2, 3, 4).
//   sub_features_covered: 5 ids mirrored above per E-STRUCT-MECH-06.
//   related_bugs: [] — pure breadth coverage of the context-panel widget
//     family; F-BUG-COVERAGE-01 already closed for this cycle.
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.panel.structure-3d]
//     source: public/packages/BiostructureViewer/src/package.ts#L854
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.panel.pdb-file-info]
//     source: public/packages/BiostructureViewer/src/package.ts#L878
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.panel.pdb-info]
//     source: public/packages/BiostructureViewer/src/package.ts#L247
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.panel.prolif]
//     source: public/packages/BiostructureViewer/src/package.ts#L262
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.panel.link-molecule-column]
//     source: public/packages/BiostructureViewer/src/package.ts#L893
//
// Selector recon-notes (class-2: live-MCP-observed 2026-06-04 on dev.datagrok.ai
//   via chrome-devtools take_snapshot + evaluate_script):
//   [name="div-section--3D-Structure"] — accordion-pane HEADER (role=button,
//     aria-expanded) for the `3D Structure` info panel. The pane container is
//     [name="pane-3D-Structure"]; the content area is
//     [name="pane-3D-Structure"] > .d4-accordion-pane-content. The accordion is
//     part of the right context panel (.grok-prop-panel), gated by
//     grok.shell.o being a DG.SemanticValue with semType=Molecule3D.
//     Click on the header flips aria-expanded false->true and mounts the
//     widget body (text marker "Biostructure Viewer:3D Structure" in
//     data-source).
//   [name="div-section--PDB-Information"] — accordion-pane HEADER for the
//     `PDB Information` panel. Same accordion-pane structure. Same registration
//     under Molecule3D (pdbFileInfoPanel: pdbFileInfoWidget(molecule.value)) and
//     PDB_ID (pdbInfoPanel: async pdbInfoWidget(pdbId)) — the two
//     registrations differ by param semType, NOT by display name, per
//     scenario .md Scenario 2 IMPORTANT note.
//   [name="div-section--Protein-Ligand-Interactions"] — accordion-pane HEADER
//     for the ProLIF panel family. Three registrations: Molecule3D +
//     hasNonWaterHetatm (pdbInteractionsWidget), Molecule3D + isAutoDockPose
//     (dockingInteractionsWidget), PDB_ID (pdbIdInteractionsWidget). Live
//     2026-06-04: 1bdq.pdb has hasNonWaterHetatm=true, isAutoDockPose=false.
//   [name="pane-3D-Structure"] / [name="pane-PDB-Information"] /
//     [name="pane-Protein-Ligand-Interactions"] — pane containers
//     (.d4-accordion-pane). Document-scope querySelector resolves; nested
//     accordion-pane-content has the widget body.
//   [name="input-host-SMILES-column"] — Link With Molecule Column widget's
//     SMILES-column ChoiceInput. Live 2026-06-04: widget root carries
//     data-source="Biostructure Viewer:Link With Molecule Column"; the
//     ChoiceInput's select.ui-input-editor has options ['', 'ligand']
//     (the only Molecule semType column in the fixture).
//
// Empirical-backing notes (mcp_observations summary; full notes in dispatch yaml):
//   - grok.shell.o = DG.SemanticValue.fromTableCell(cell) IS the canonical
//     way to surface a Molecule3D / PDB_ID value to the right context-panel
//     accordion (not df.currentCell directly — the cell click maps to the
//     COLUMN object, semType-gated panels then fire on the cell VALUE). Live
//     2026-06-04: setting shell.o to a SemanticValue with semType=Molecule3D
//     surfaces [pane-3D-Structure, pane-PDB-Information, pane-Protein-Ligand-
//     Interactions]; setting shell.o to a SemanticValue with semType=PDB_ID
//     surfaces [pane-PDB-id-viewer, pane-PDB-Information, pane-Protein-
//     Ligand-Interactions].
//   - 1bdq.pdb empirically has hasNonWaterHetatm=true (21 ligands found:
//     IM1 R 400, IM1 C 401, ...) and isAutoDockPose=false. The Molecule3D
//     ProLIF widget renders verbatim "21 ligands found: IM1 ..." with
//     data-source="Biostructure Viewer:Protein-Ligand Interactions"
//     (pdbInteractionsWidget). The Docking-flavour registration
//     (dockingInteractionsWidget) is deferred per the scenario .md cited
//     deferral (no AutoDock pose fixture available without a co-installed
//     Docking package + receptor AppData).
//   - PDB_ID = '1QBS' triggers pdbInfoPanel (async pdbInfoWidget): verbatim
//     rendered text begins "GeneralDescriptionHIV-1 PROTEASE INHIBITORS
//     WIIH LOW NANOMOLAR POTENCYClassificationASPARTYL PROTEASEMethodX-ray
//     Resolution1.8 ÅPDB URLrcsb.org/structure/1QBS...".
//   - PDB_ID = '1QBS' triggers pdbIdInteractionsWidget (RCSB fetchProxy):
//     verbatim rendered text begins "3 ligands found: DMP A 323, CSO A 67,
//     CSO B 67..." with data-source="Biostructure Viewer:Protein-Ligand
//     Interactions". This is the third ProLIF registration's resolution path.
//   - Molecule3D `pdbFileInfoPanel` returns DG.Widget with verbatim rendered
//     text beginning "GeneralAtom Count2240ChainsA, B, C, D, E, F, ...".
//     data-source="Biostructure Viewer:PDB Information" (the panel name is
//     shared across the two Molecule3D / PDB_ID registrations; the
//     differentiator is the body text content — file-info vs RCSB metadata).
//   - mol3dAtomPickerLinkWidget(mol3DCol) returns DG.Widget whose root
//     carries class="d4-flex-col ui-div", data-source="Biostructure Viewer:
//     Link With Molecule Column", and contains a ChoiceInput
//     [name="input-host-SMILES-column"] whose select.ui-input-editor options
//     are ['', 'ligand'] (filtered to Molecule semType columns in the host
//     DataFrame). Function registration shape per DG.Func.find: inputs
//     [{name:'mol3DCol', type:'column', semType:'Molecule3D'}], outputs
//     [{name:'result', type:'widget'}], friendlyName='Link With Molecule
//     Column'.
//   - WebGL-uncertain dev environment caveat: the Mol* engine init logs
//     "Could not create a WebGL rendering context" (same as sibling
//     biostructure-viewer-spec.ts / property-surface-extension-spec.ts /
//     ngl-viewer-extension-spec.ts). 3D Structure panel's WebGL canvas
//     pixel content is therefore NOT pixel-asserted (SR-01 below); the
//     widget MOUNT + .bsv-container-info-panel container presence + data-
//     source="Biostructure Viewer:3D Structure" data attribute ARE
//     deterministically asserted.
//   - 3D Structure widget container live 2026-06-04: header click on
//     [name="div-section--3D-Structure"] flips aria-expanded false->true;
//     the pane content area mounts a `.bsv-container-info-panel` element
//     containing `.d4-molstar-viewer` + `.msp-plugin` DOM. The
//     bsv-container-info-panel class IS the assert anchor (atlas Scenario 1
//     "the panel root carries the bsv-container-info-panel class added by
//     structure3D after render").
//
// DOM-driving rationale (>=1 DOM-driving call REQUIRED for target_layer:
//   playwright per E-LAYER-COMPLIANCE-01):
//   - page.locator('[name="Browse"]').waitFor — DOM readiness anchor after
//     login (spec-login.ts pattern; same as sibling specs).
//   - DOM header clicks on [name="div-section--3D-Structure"] /
//     [name="div-section--PDB-Information"] / [name="div-section--Protein-
//     Ligand-Interactions"] — Scenarios 1, 2, 3, 4 expand the relevant
//     accordion pane. A real user opens these panels exactly this way
//     (click the section header).
//   - DOM presence verification via page.locator(...).waitFor on
//     [name="pane-3D-Structure"], [name="pane-PDB-Information"],
//     [name="pane-Protein-Ligand-Interactions"], and [name="pane-PDB-id-
//     viewer"] (PDB_ID flow).
//
// Paradigm rationale (mostly DOM-driven panel expansion with deterministic
//   text-content assertions plus a JS-API function-call surface for
//   Scenario 5):
//   - Scenarios 1-4 exercise the right-hand context panel — DOM accordion
//     pane mount + text content assertion. The widget bodies' rendered text
//     is deterministically introspectable via .textContent on the
//     .d4-accordion-pane-content node (the panel functions are NOT WebGL-
//     dependent for their text output; only Scenario 1's 3D Structure mounts
//     a Mol* canvas, and the widget MOUNT is asserted at the container-
//     class level rather than the canvas pixel level).
//   - Scenario 5 exercises a function-widget — JS-API call surface is the
//     contract per the scenario .md IMPORTANT note ("the scenario exercises
//     only the direct JS-API call surface (the contract)"). The widget root
//     IS attached to a transient host inside the page.evaluate context for
//     DOM-shape verification (data-source attribute, SMILES column picker
//     [name="input-host-SMILES-column"]).
//
// Scope reductions (per scenario .md Setup + Deferral cite + Notes):
//   SR-01 — Scenario 1 WebGL canvas pixel content NOT asserted. The Mol*
//     engine fails to create a WebGL context in the dev recon environment
//     (same pattern as biostructure-viewer-spec.ts SR / property-surface-
//     extension-spec.ts SR-01 / ngl-viewer-extension-spec.ts SR-04). The
//     widget MOUNT (bsv-container-info-panel container present + Mol*
//     plugin DOM root present + data-source attribute exact) IS
//     deterministically asserted; the "non-blank canvas" pixel-content
//     expectation in the scenario .md Expected bullet is the visual-
//     judgment slice that the ui-affordance manual-only split rule routes
//     to a -ui.md companion (precedent: sibling SR-01 / SR-04 in the
//     section). The widget's correct registration + mount + container-
//     class invariant IS the structural assertion the bug-invariant of the
//     panel hinges on.
//   SR-02 — Scenario 4 Docking-flavour (isAutoDockPose) sub-step
//     deferred per the scenario .md explicit Deferral note. The fixture
//     corpus does not include an AutoDock pose Molecule3D row (a pose PDB
//     with a `REMARK ... binding energy` line and a corresponding receptor
//     under `System:AppData/Docking/targets/`). The deferral cites the
//     real technical dependency (a co-installed Docking package + receptor
//     AppData files). Live 2026-06-04: isAutoDockPose returns false on
//     1bdq.pdb (the section's only Molecule3D fixture). The first ProLIF
//     registration (pdbInteractionsWidget for Molecule3D + hasNonWater-
//     Hetatm) AND the third (pdbIdInteractionsWidget for PDB_ID) ARE
//     fully exercised. The bug-invariant of "three condition-gated
//     registrations under the same panel name" is structurally captured
//     via DG.Func.find probes on the panel function-registry surface
//     (all three panel registrations exist with the documented
//     condition: predicate strings).
//   SR-03 — Scenario 1 re-mount sub-step "Switching the current cell to a
//     different Molecule3D value re-mounts the panel against the new cell
//     value" — asserted at the JS-API panel-function re-invocation level
//     (a fresh structure3D widget produced for a different cell value
//     yields a fresh widget root) rather than at the live-DOM accordion
//     re-render level. The widget body IS WebGL-uncertain on re-mount in
//     this recon environment; the function-contract is asserted in lieu
//     of the canvas pixel re-render.
//   SR-04 — Scenario 5 widget picker interaction ("Selecting a SMILES
//     column in the picker updates the widget's internal state") — the
//     internal-state contract IS in-memory per the scenario .md Expected
//     bullet; we assert the picker structure (input-host-SMILES-column
//     present + select options filtered to Molecule semType columns) but
//     do NOT drive a UI selection event. The bug-invariant (the widget
//     surfaces a SMILES-column picker tied to the Molecule3D column
//     reference) is structurally captured at the DOM shape level. The
//     scenario .md IMPORTANT note explicitly scopes Scenario 5 to "the
//     direct JS-API call surface (the contract)".
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

test('BiostructureViewer — context-panel widgets extension (3D Structure / PDB Information x2 / ProLIF x3 / Link With Molecule Column)', async ({page}) => {
  test.setTimeout(900_000);
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

  // DOM-driving readiness anchor (E-LAYER-COMPLIANCE-01 slot 1).
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  // Prepare the shared fixture DataFrame (Setup step 3): three columns,
  // semType-tagged for Molecule3D / PDB_ID / Molecule respectively. Reused
  // across Scenarios 1-5.
  await page.evaluate(async (pdbPath) => {
    grok.shell.closeAll();
    await new Promise((r) => setTimeout(r, 1500));
    const pdbContent = await grok.dapi.files.readAsText(pdbPath);
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['row-1', 'row-2', 'row-3']),
      DG.Column.fromStrings('structure', [pdbContent, pdbContent, pdbContent]),
      DG.Column.fromStrings('pdb_id', ['1QBS', '1BNA', '2J1X']),
      DG.Column.fromStrings('ligand', ['CCO', 'CCC', 'CCN']),
    ]);
    try { df.col('structure').semType = 'Molecule3D'; } catch (_e) { /* best-effort */ }
    try { df.col('pdb_id').semType = 'PDB_ID'; } catch (_e) { /* best-effort */ }
    try { df.col('ligand').semType = 'Molecule'; } catch (_e) { /* best-effort */ }
    df.name = 'context-panel-fixture';
    grok.shell.addTableView(df);
    await new Promise((r) => setTimeout(r, 2500));
  }, samplePdbPath);

  try {
    // ========================================================================
    // SCENARIO 1 — `3D Structure` panel widget for a Molecule3D cell.
    //
    //   DOM-driving paradigm: surface the Molecule3D SemanticValue to the
    //   right-context-panel accordion via grok.shell.o; click the
    //   [name="div-section--3D-Structure"] HEADER to expand the pane (real
    //   user opens the section this way); assert the widget MOUNT via the
    //   .bsv-container-info-panel container class (atlas-prescribed marker
    //   added by structure3D() after render). Canvas pixel content NOT
    //   pixel-asserted per SR-01.
    //
    //   sub_features_covered: biostructure.panel.structure-3d.
    // ========================================================================

    await softStep('Scenario 1 — 3D Structure pane mounts for Molecule3D cell; container class .bsv-container-info-panel present', async () => {
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const df = tv.dataFrame;
        // Step 2: surface Molecule3D cell value to the right context panel.
        const cell = df.cell(0, 'structure');
        grok.shell.o = DG.SemanticValue.fromTableCell(cell);
        await new Promise((r) => setTimeout(r, 3500));
        return {
          pane3DPresent: !!document.querySelector('[name="pane-3D-Structure"]'),
          header3DPresent: !!document.querySelector('[name="div-section--3D-Structure"]'),
          header3DAriaBefore: document.querySelector('[name="div-section--3D-Structure"]')?.getAttribute('aria-expanded') ?? null,
        };
      });

      // Step 4: pane container present in the DOM.
      await page.locator('[name="pane-3D-Structure"]').waitFor({timeout: 30_000});

      expect(res.pane3DPresent).toBe(true);
      expect(res.header3DPresent).toBe(true);
      // Pane collapsed by default — header aria-expanded='false'.
      expect(res.header3DAriaBefore).toBe('false');

      // DOM-driving slot: click the accordion section header (real-user
      // expansion action). E-LAYER-COMPLIANCE-01 slot 2.
      await page.locator('[name="div-section--3D-Structure"]').click({timeout: 10_000});
      await page.waitForTimeout(5000); // allow widget mount + Mol* engine attempt.

      const after = await page.evaluate(() => {
        const pane3D = document.querySelector('[name="pane-3D-Structure"]');
        const content = pane3D?.querySelector('.d4-accordion-pane-content');
        return {
          header3DAriaAfter: document.querySelector('[name="div-section--3D-Structure"]')?.getAttribute('aria-expanded') ?? null,
          contentChildren: content?.children.length ?? 0,
          hasBsvContainer: !!content?.querySelector('.bsv-container-info-panel'),
          dataSource: content?.querySelector('[data-source]')?.getAttribute('data-source') ?? null,
        };
      });

      // Step 5: widget MOUNT verified at the structural level.
      expect(after.header3DAriaAfter).toBe('true');
      expect(after.contentChildren).toBeGreaterThan(0);
      // Atlas-prescribed marker (Scenario 1 Expected): the panel root carries
      // the `bsv-container-info-panel` class added by structure3D() after
      // render. The class IS added inside the renderer.createViewer().then()
      // callback per package.ts#L867 — its presence is the structural
      // mount-success invariant.
      expect(after.hasBsvContainer).toBe(true);
      // The widget root carries data-source="Biostructure Viewer:3D Structure"
      // — the registration tag that identifies the panel function.
      expect(after.dataSource).toBe('Biostructure Viewer:3D Structure');
    });

    await softStep('Scenario 1 — re-mount: structure3D(SemanticValue) yields a fresh DG.Widget for a different cell value (SR-03)', async () => {
      // SR-03: live-DOM accordion re-render is WebGL-uncertain on re-mount in
      // this recon environment; we assert the function-contract re-invocation
      // surface instead (the same registration the cell switch dispatches to).
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const df = tv.dataFrame;
        const fns = DG.Func.find({name: 'structure3D', package: 'BiostructureViewer'});
        const fn = fns && fns[0];
        const cell0 = df.cell(0, 'structure');
        const cell1 = df.cell(1, 'structure');
        const w0 = await fn?.apply({molecule: DG.SemanticValue.fromTableCell(cell0)});
        const w1 = await fn?.apply({molecule: DG.SemanticValue.fromTableCell(cell1)});
        return {
          fnFound: !!fn,
          w0Present: !!w0,
          w1Present: !!w1,
          // Distinct widget roots (DOM nodes are not reference-equal across
          // invocations even when input value happens to be identical).
          rootsDistinct: !!w0 && !!w1 && w0.root !== w1.root,
        };
      });

      expect(res.fnFound).toBe(true);
      expect(res.w0Present).toBe(true);
      expect(res.w1Present).toBe(true);
      expect(res.rootsDistinct).toBe(true);
    });

    // ========================================================================
    // SCENARIO 2 — `PDB Information` panel on a Molecule3D cell (header
    //   / file-info from PDB string via pdbFileInfoWidget).
    //
    //   DOM-driving paradigm: header click on
    //   [name="div-section--PDB-Information"] expands the pane; the
    //   widget body text begins with verbatim "GeneralAtom Count2240Chains
    //   A, B, C, D, E, F, ..." (live-observed 2026-06-04). The same display
    //   name is reused by Scenario 3's PDB_ID-flavoured registration — the
    //   two registrations are differentiated by the param semType the
    //   panel decorator gates on, NOT by the display name (scenario .md
    //   IMPORTANT note).
    //
    //   sub_features_covered: biostructure.panel.pdb-file-info.
    // ========================================================================

    await softStep('Scenario 2 — PDB Information pane (Molecule3D) renders pdbFileInfoWidget content', async () => {
      const surface = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const df = tv.dataFrame;
        const cell = df.cell(0, 'structure');
        grok.shell.o = DG.SemanticValue.fromTableCell(cell);
        await new Promise((r) => setTimeout(r, 3000));
        return {
          panePdbPresent: !!document.querySelector('[name="pane-PDB-Information"]'),
        };
      });

      await page.locator('[name="pane-PDB-Information"]').waitFor({timeout: 30_000});
      expect(surface.panePdbPresent).toBe(true);

      // DOM-driving slot: expand PDB Information pane via header click.
      await page.locator('[name="div-section--PDB-Information"]').click({timeout: 10_000});
      await page.waitForTimeout(3500);

      const after = await page.evaluate(() => {
        const pane = document.querySelector('[name="pane-PDB-Information"]');
        const content = pane?.querySelector('.d4-accordion-pane-content');
        return {
          ariaAfter: document.querySelector('[name="div-section--PDB-Information"]')?.getAttribute('aria-expanded') ?? null,
          textLen: content?.textContent?.length ?? 0,
          textSnippet: (content?.textContent || '').slice(0, 200),
          // Live 2026-06-04: pdbFileInfoWidget output begins
          // "GeneralAtom Count2240ChainsA, B, C, D, E, F, ...".
          hasGeneralAtomCount: (content?.textContent || '').includes('Atom Count'),
          hasChains: (content?.textContent || '').includes('Chains'),
        };
      });

      expect(after.ariaAfter).toBe('true');
      expect(after.textLen).toBeGreaterThan(50);
      // Bug-invariant of the pdb-file-info widget: file-info derived from the
      // PDB string (HEADER / ATOM count / chain listing). The "Atom Count"
      // + "Chains" markers are the deterministic surface.
      expect(after.hasGeneralAtomCount).toBe(true);
      expect(after.hasChains).toBe(true);
    });

    // ========================================================================
    // SCENARIO 3 — `PDB Information` panel on a PDB_ID cell (async
    //   pdbInfoWidget assembling RCSB metadata).
    //
    //   DOM-driving paradigm: surface a PDB_ID SemanticValue (1QBS) to the
    //   right context panel; expand the PDB Information pane via header
    //   click. The pdbInfoWidget IS async (returns Promise<DG.Widget>);
    //   the body assembles after the widget mount with verbatim text
    //   beginning "GeneralDescriptionHIV-1 PROTEASE INHIBITORS WIIH LOW
    //   NANOMOLAR POTENCYClassificationASPARTYL PROTEASEMethodX-rayResolution
    //   1.8 Å..." (live-observed 2026-06-04). The pdbInfoWidget surfaces the
    //   PDB URL anchor "rcsb.org/structure/1QBS" — the load-bearing PDB-ID-
    //   addressed metadata.
    //
    //   sub_features_covered: biostructure.panel.pdb-info.
    // ========================================================================

    await softStep('Scenario 3 — PDB Information pane (PDB_ID) renders async pdbInfoWidget metadata', async () => {
      const surface = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const df = tv.dataFrame;
        const cell = df.cell(0, 'pdb_id');
        grok.shell.o = DG.SemanticValue.fromTableCell(cell);
        await new Promise((r) => setTimeout(r, 3500));
        return {
          panePdbInfoIdPresent: !!document.querySelector('[name="pane-PDB-Information"]'),
          panePdbIdViewerPresent: !!document.querySelector('[name="pane-PDB-id-viewer"]'),
        };
      });

      // Both the PDB Information (PDB_ID flavour) AND the PDB id viewer
      // panes register on a PDB_ID SemanticValue.
      expect(surface.panePdbInfoIdPresent).toBe(true);
      expect(surface.panePdbIdViewerPresent).toBe(true);

      // DOM-driving slot: expand PDB Information (PDB_ID-flavour now).
      await page.locator('[name="div-section--PDB-Information"]').click({timeout: 10_000});
      // Bounded wait for the async pdbInfoWidget body (RCSB metadata
      // assembly — empirically settled in 5-8s on the recon env).
      await page.waitForTimeout(8000);

      const after = await page.evaluate(() => {
        const pane = document.querySelector('[name="pane-PDB-Information"]');
        const content = pane?.querySelector('.d4-accordion-pane-content');
        const text = content?.textContent || '';
        return {
          ariaAfter: document.querySelector('[name="div-section--PDB-Information"]')?.getAttribute('aria-expanded') ?? null,
          textLen: text.length,
          textSnippet: text.slice(0, 200),
          // Empirically observed 2026-06-04 on 1QBS:
          hasDescription: text.includes('Description'),
          hasRcsbUrl: text.includes('rcsb.org/structure/1QBS'),
          // Per scenario Expected: "No `Could not fetch PDB <id>` error string
          // in the widget body" — the failure-mode marker we negative-assert.
          noFetchError: !text.includes('Could not fetch PDB'),
        };
      });

      expect(after.ariaAfter).toBe('true');
      expect(after.textLen).toBeGreaterThan(50);
      expect(after.hasDescription).toBe(true);
      // Load-bearing PDB-ID-addressed metadata marker.
      expect(after.hasRcsbUrl).toBe(true);
      expect(after.noFetchError).toBe(true);
    });

    // ========================================================================
    // SCENARIO 4 — `Protein-Ligand Interactions` (ProLIF) — three
    //   condition-gated registrations under the same panel name resolve per
    //   cell semType + condition.
    //
    //   Path A (Molecule3D + hasNonWaterHetatm): pdbInteractionsWidget →
    //     makeProlifWidget({protein}); 1bdq.pdb has hasNonWaterHetatm=true
    //     (21 ligands found: IM1 R 400, ...). Live-DOM mount + verbatim
    //     "ligands found" marker.
    //   Path B (Molecule3D + isAutoDockPose): dockingInteractionsWidget;
    //     deferred per scenario .md Deferral (no AutoDock pose fixture
    //     without co-installed Docking package + receptor AppData). The
    //     registration's existence IS asserted via DG.Func.find on the
    //     panel function-registry surface (SR-02).
    //   Path C (PDB_ID): pdbIdInteractionsWidget — fetches RCSB PDB via
    //     dapi.fetchProxy and renders verbatim "3 ligands found: DMP A
    //     323, CSO A 67, CSO B 67" for 1QBS.
    //
    //   sub_features_covered: biostructure.panel.prolif.
    // ========================================================================

    await softStep('Scenario 4 Path A — ProLIF panel (Molecule3D + hasNonWaterHetatm) renders pdbInteractionsWidget LigNetwork', async () => {
      const surface = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const df = tv.dataFrame;
        const cell = df.cell(0, 'structure');
        grok.shell.o = DG.SemanticValue.fromTableCell(cell);
        await new Promise((r) => setTimeout(r, 3000));
        // Live-confirm the gating predicate on the fixture.
        const pdbContent = await grok.dapi.files.readAsText('System:AppData/BiostructureViewer/samples/1bdq.pdb');
        const hasHet = await grok.functions.call('BiostructureViewer:hasNonWaterHetatm', {molecule: pdbContent});
        const isAuto = await grok.functions.call('BiostructureViewer:isAutoDockPose', {molecule: pdbContent});
        return {
          paneProlifPresent: !!document.querySelector('[name="pane-Protein-Ligand-Interactions"]'),
          hasNonWaterHetatm: hasHet,
          isAutoDockPose: isAuto,
        };
      });

      // Gate-precondition empirically established (1bdq.pdb has
      // hasNonWaterHetatm=true; isAutoDockPose=false). Path A is the
      // resolution path for this fixture.
      expect(surface.paneProlifPresent).toBe(true);
      expect(surface.hasNonWaterHetatm).toBe(true);
      // SR-02: isAutoDockPose=false on the fixture; Path B deferred per
      // scenario .md.
      expect(surface.isAutoDockPose).toBe(false);

      // DOM-driving slot: expand the ProLIF pane via header click.
      await page.locator('[name="div-section--Protein-Ligand-Interactions"]').click({timeout: 10_000});
      // pdbInteractionsWidget is async (await makeProlifWidget); allow a
      // bounded wait for LigNetwork assembly (empirically 4-8s on the recon
      // env).
      await page.waitForTimeout(8000);

      const after = await page.evaluate(() => {
        const pane = document.querySelector('[name="pane-Protein-Ligand-Interactions"]');
        const content = pane?.querySelector('.d4-accordion-pane-content');
        const text = content?.textContent || '';
        return {
          ariaAfter: document.querySelector('[name="div-section--Protein-Ligand-Interactions"]')?.getAttribute('aria-expanded') ?? null,
          textLen: text.length,
          textSnippet: text.slice(0, 200),
          dataSource: content?.querySelector('[data-source]')?.getAttribute('data-source') ?? null,
          // Live 2026-06-04 on 1bdq: "21 ligands found: IM1 R 400, ...".
          hasLigandsFound: /\d+\s+ligands?\s+found/.test(text),
          hasIM1: text.includes('IM1'),
        };
      });

      expect(after.ariaAfter).toBe('true');
      expect(after.textLen).toBeGreaterThan(50);
      // Atlas-prescribed widget body marker (LigNetwork interactions diagram).
      expect(after.hasLigandsFound).toBe(true);
      // 1bdq's ligand IM1 — fixture-specific verbatim marker.
      expect(after.hasIM1).toBe(true);
      // data-source identifies the panel function (shared across the three
      // registrations — the gating predicate is what differentiates them).
      expect(after.dataSource).toBe('Biostructure Viewer:Protein-Ligand Interactions');
    });

    await softStep('Scenario 4 Path B — DG.Func registration probe for dockingInteractionsWidget (SR-02: AutoDock pose fixture not available)', async () => {
      // SR-02: the Docking-flavour registration is structurally exercised at
      // the DG.Func.find / panel decorator registration level. The fixture
      // corpus does not include an AutoDock pose Molecule3D row (real
      // technical dependency on a co-installed Docking package + receptor
      // AppData files). The registration's existence — with the correct
      // condition predicate and input semType — IS the structural surface
      // the bug-invariant of "three condition-gated registrations" hinges on.
      //
      // Probe shape note: panel decorators register with role at
      // fn.options.role === 'panel' (NOT in fn.tags). fn.friendlyName carries
      // the display name 'Protein-Ligand Interactions' (the @panel({name:...})
      // value) while fn.name carries the export function name
      // 'dockingInteractionsWidget'. fn.options.condition carries the gating
      // predicate string 'BiostructureViewer:isAutoDockPose(molecule)' —
      // the load-bearing differentiator across the three same-display-name
      // registrations. Live-verified 2026-06-04 on dev.datagrok.ai via
      // chrome-devtools MCP evaluate_script.
      const res = await page.evaluate(() => {
        const fns = DG.Func.find({name: 'dockingInteractionsWidget', package: 'BiostructureViewer'});
        const fn = fns && fns[0];
        const inputs = (fn && fn.inputs) ? fn.inputs.map((i: any) => ({
          name: i.name,
          type: i.propertyType,
          semType: i.options?.semType ?? null,
        })) : [];
        let tagsArr: string[] = [];
        try { tagsArr = (fn && fn.tags) ? Array.from(fn.tags as Iterable<string>) : []; } catch (_e) { tagsArr = []; }
        const optionsRole = fn?.options?.role ?? null;
        const friendlyName = fn?.friendlyName ?? null;
        const condition = fn?.options?.condition ?? null;
        return {
          registered: !!fn,
          inputCount: inputs.length,
          inputSemTypes: inputs.map((i: any) => i.semType),
          // Panel role exposed via fn.options.role (canonical) OR fn.tags
          // (legacy fallback). The two-source disjunction makes the probe
          // robust to API surface evolution.
          isPanel: optionsRole === 'panel' || tagsArr.includes('panel'),
          optionsRole,
          friendlyName,
          condition,
        };
      });

      // The Docking-flavour ProLIF registration EXISTS at the function-
      // registry level (the panel decorator is processed at package load).
      expect(res.registered).toBe(true);
      expect(res.inputCount).toBeGreaterThan(0);
      // The molecule input is gated to semType=Molecule3D — same as the
      // hasNonWaterHetatm registration; the differentiator is the
      // BiostructureViewer:isAutoDockPose condition predicate, separately
      // verified above (returned false on the fixture).
      expect(res.inputSemTypes).toContain('Molecule3D');
      // Panel role assertion — accepts either canonical surface
      // (fn.options.role === 'panel'; observed 2026-06-04) or legacy
      // tag-list surface.
      expect(res.isPanel).toBe(true);
      // friendlyName carries the @panel({name: 'Protein-Ligand Interactions'})
      // display value — confirms the panel decorator was processed and bound
      // its display name to the export function dockingInteractionsWidget.
      expect(res.friendlyName).toBe('Protein-Ligand Interactions');
      // The condition predicate is the load-bearing differentiator across
      // the three same-display-name 'Protein-Ligand Interactions'
      // registrations (hasNonWaterHetatm / isAutoDockPose / no-condition
      // PDB_ID). The Docking-flavour registration MUST carry the
      // isAutoDockPose predicate per package.ts#L297.
      expect(res.condition).toBe('BiostructureViewer:isAutoDockPose(molecule)');
    });

    await softStep('Scenario 4 Path C — ProLIF panel (PDB_ID) renders pdbIdInteractionsWidget via RCSB fetchProxy', async () => {
      const surface = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const df = tv.dataFrame;
        const cell = df.cell(0, 'pdb_id');
        grok.shell.o = DG.SemanticValue.fromTableCell(cell);
        await new Promise((r) => setTimeout(r, 3500));
        return {
          paneProlifPresent: !!document.querySelector('[name="pane-Protein-Ligand-Interactions"]'),
        };
      });

      expect(surface.paneProlifPresent).toBe(true);

      // DOM-driving slot: expand the ProLIF pane (PDB_ID flavour) via header
      // click. fetchProxy to files.rcsb.org takes longer than Path A — allow
      // a bounded 15s window (empirically 10-12s on the recon env).
      await page.locator('[name="div-section--Protein-Ligand-Interactions"]').click({timeout: 10_000});
      await page.waitForTimeout(15000);

      const after = await page.evaluate(() => {
        const pane = document.querySelector('[name="pane-Protein-Ligand-Interactions"]');
        const content = pane?.querySelector('.d4-accordion-pane-content');
        const text = content?.textContent || '';
        return {
          ariaAfter: document.querySelector('[name="div-section--Protein-Ligand-Interactions"]')?.getAttribute('aria-expanded') ?? null,
          textLen: text.length,
          textSnippet: text.slice(0, 200),
          dataSource: content?.querySelector('[data-source]')?.getAttribute('data-source') ?? null,
          // Live 2026-06-04 on 1QBS: "3 ligands found: DMP A 323, CSO A 67,
          // CSO B 67".
          hasLigandsFound: /\d+\s+ligands?\s+found/.test(text),
          hasDmp: text.includes('DMP'),
          // Per scenario Expected: "If the fetch is non-ok, the widget body
          // reads `Could not fetch PDB <id>`; treat that as flake-retry once".
          // Here we assert the absence of the failure marker (success path).
          noFetchError: !text.includes('Could not fetch PDB'),
        };
      });

      expect(after.ariaAfter).toBe('true');
      expect(after.textLen).toBeGreaterThan(20);
      expect(after.hasLigandsFound).toBe(true);
      // 1QBS's known ligand DMP — fixture-specific verbatim marker proving
      // the RCSB fetch succeeded and ProLIF rendered against the fetched
      // structure.
      expect(after.hasDmp).toBe(true);
      // Same data-source registration tag as Path A — confirms the three
      // ProLIF registrations share the panel display name.
      expect(after.dataSource).toBe('Biostructure Viewer:Protein-Ligand Interactions');
      expect(after.noFetchError).toBe(true);
    });

    // ========================================================================
    // SCENARIO 5 — `Link With Molecule Column` — cross-package widget
    //   injection (function-widget).
    //
    //   Paradigm: direct JS-API call surface per the scenario .md IMPORTANT
    //   note. grok.functions.call('BiostructureViewer:mol3dAtomPickerLinkWidget',
    //   {mol3DCol}) returns a DG.Widget whose root carries data-source=
    //   "Biostructure Viewer:Link With Molecule Column" and contains a
    //   ChoiceInput [name="input-host-SMILES-column"] tied to the host
    //   DataFrame's Molecule semType columns (here: `ligand`). The function
    //   is registered with friendlyName='Link With Molecule Column' and
    //   outputs [{name:'result', type:'widget'}] per the
    //   @grok.decorators.func({name: 'Link With Molecule Column', outputs:
    //   [{name:'result', type:'widget'}]}) declaration at package.ts#L889.
    //
    //   sub_features_covered: biostructure.panel.link-molecule-column.
    // ========================================================================

    await softStep('Scenario 5 — mol3dAtomPickerLinkWidget(mol3DCol) returns DG.Widget with SMILES-column ChoiceInput', async () => {
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const df = tv.dataFrame;

        // Step 2: invoke directly via JS API (scenario .md IMPORTANT note —
        // this is the contract surface).
        let callErr: string | null = null;
        let widgetOk = false;
        let widgetRootClass: string | null = null;
        let widgetDataSource: string | null = null;
        let smilesPickerPresent = false;
        let smilesPickerOptions: string[] | null = null;
        try {
          const w = await grok.functions.call(
            'BiostructureViewer:mol3dAtomPickerLinkWidget',
            {mol3DCol: df.col('structure')},
          );
          widgetOk = !!w;
          widgetRootClass = w?.root?.className ?? null;
          widgetDataSource = w?.root?.querySelector('[data-source]')?.getAttribute('data-source')
            ?? w?.root?.getAttribute?.('data-source')
            ?? null;
          // Step 3-4: SMILES column picker UI inside the widget.
          smilesPickerPresent = !!w?.root?.querySelector('[name="input-host-SMILES-column"]');
          const select = w?.root?.querySelector('[name="input-host-SMILES-column"] select.ui-input-editor');
          if (select) {
            smilesPickerOptions = Array.from(select.options).map((o: any) => o.value);
          }
        } catch (e: any) { callErr = String(e?.message ?? e); }

        // Function-registry contract check.
        const fns = DG.Func.find({name: 'mol3dAtomPickerLinkWidget', package: 'BiostructureViewer'});
        const fn = fns && fns[0];
        const inputs = (fn && fn.inputs) ? fn.inputs.map((i: any) => ({
          name: i.name,
          type: i.propertyType,
          semType: i.options?.semType ?? null,
        })) : [];
        const outputs = (fn && fn.outputs) ? fn.outputs.map((o: any) => ({
          name: o.name,
          type: o.propertyType,
        })) : [];

        return {
          callErr,
          widgetOk,
          widgetRootClass,
          widgetDataSource,
          smilesPickerPresent,
          smilesPickerOptions,
          fnRegistered: !!fn,
          fnFriendlyName: fn?.friendlyName ?? null,
          inputCount: inputs.length,
          inputNames: inputs.map((i: any) => i.name),
          inputSemTypes: inputs.map((i: any) => i.semType),
          outputs,
        };
      });

      // Step 2 Expected: the call resolves to a DG.Widget; no exception.
      expect(res.callErr).toBe(null);
      expect(res.widgetOk).toBe(true);

      // The widget root data-source identifies the registration. The
      // attribute is set by Datagrok's function-widget host wrapper around
      // the returned DG.Widget; live-verified 2026-06-04.
      expect(res.widgetDataSource).toBe('Biostructure Viewer:Link With Molecule Column');

      // Step 3-4 Expected: SMILES-column picker scoped to Molecule semType
      // columns in the host DataFrame.
      expect(res.smilesPickerPresent).toBe(true);
      // SR-04: the in-memory selection contract is asserted at the picker
      // structure level (options include the Molecule-semType column
      // 'ligand'). We do NOT drive a UI selection event — the IMPORTANT
      // note scopes the test to the contract surface.
      expect(res.smilesPickerOptions).not.toBe(null);
      expect((res.smilesPickerOptions || []).includes('ligand')).toBe(true);

      // Function-registry contract: register name='Link With Molecule
      // Column', outputs [{name:'result', type:'widget'}], inputs [{name:
      // 'mol3DCol', type:'column', semType:'Molecule3D'}] per
      // package.ts#L889-L893.
      expect(res.fnRegistered).toBe(true);
      expect(res.fnFriendlyName).toBe('Link With Molecule Column');
      expect(res.inputCount).toBe(1);
      expect(res.inputNames).toEqual(['mol3DCol']);
      expect(res.inputSemTypes).toEqual(['Molecule3D']);
      expect(res.outputs).toEqual([{name: 'result', type: 'widget'}]);
    });
  } finally {
    // Cleanup.
    await page.evaluate(() => { grok.shell.closeAll(); });
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
