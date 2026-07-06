/* ---
sub_features_covered: [biostructure.panel.link-molecule-column, biostructure.panel.pdb-file-info, biostructure.panel.pdb-info, biostructure.panel.prolif, biostructure.panel.structure-3d]
--- */
// Context-panel widgets: 3D Structure / PDB Information (Molecule3D + PDB_ID) / ProLIF (3 gated
// registrations) / Link With Molecule Column. Panes are surfaced via grok.shell.o =
// SemanticValue.fromTableCell(cell); accordion headers ([name="div-section--..."]) expand on click.
// WebGL canvas pixels are not asserted (CI WebGL is unreliable) — widget mount + data-source + text
// content are. 1bdq.pdb: hasNonWaterHetatm=true, isAutoDockPose=false (Docking ProLIF path deferred).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';
// A second, genuinely-different structure so the re-mount step (SR-03) exercises a real value change.
const samplePdbPath2 = 'System:AppData/BiostructureViewer/samples/protease.pdb';

test('BiostructureViewer / context-panel widgets extension (3D Structure / PDB Information x2 / ProLIF x3 / Link With Molecule Column)', async ({page}) => {
  test.setTimeout(180_000);
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

  // Shared fixture DataFrame: columns semType-tagged Molecule3D / PDB_ID / Molecule. Reused by all.
  // Row 0 is 1bdq (all downstream scenarios read row 0); rows 1-2 hold a different structure so the
  // re-mount step operates on genuinely distinct cell values.
  await page.evaluate(async (paths) => {
    grok.shell.closeAll();
    const pdbContent = await grok.dapi.files.readAsText(paths.pdb1);
    const pdbContent2 = await grok.dapi.files.readAsText(paths.pdb2);
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['row-1', 'row-2', 'row-3']),
      DG.Column.fromStrings('structure', [pdbContent, pdbContent2, pdbContent2]),
      DG.Column.fromStrings('pdb_id', ['1QBS', '1BNA', '2J1X']),
      DG.Column.fromStrings('ligand', ['CCO', 'CCC', 'CCN']),
    ]);
    try { df.col('structure').semType = 'Molecule3D'; } catch (_e) { /* best-effort */ }
    try { df.col('pdb_id').semType = 'PDB_ID'; } catch (_e) { /* best-effort */ }
    try { df.col('ligand').semType = 'Molecule'; } catch (_e) { /* best-effort */ }
    df.name = 'context-panel-fixture';
    grok.shell.addTableView(df);
  }, {pdb1: samplePdbPath, pdb2: samplePdbPath2});

  await page.waitForFunction(() =>
    grok.shell.tv?.dataFrame?.name === 'context-panel-fixture' &&
    grok.shell.tv.dataFrame.columns.length === 4, null, {timeout: 15_000});

  // Prerequisite guard: the 3D Structure / PDB Information panels live in BiostructureViewer.
  // If they aren't registered on this deployment, skip fast with a clear message instead of
  // hanging the full 180s test budget on accordion panes that can never appear.
  const bsvPanelsReady = await page.evaluate(() => {
    const has = (name: string) => (DG.Func.find({package: 'BiostructureViewer', name}) || []).length > 0;
    return has('structure3D') && has('pdbFileInfoPanel');
  });
  test.skip(!bsvPanelsReady, 'BiostructureViewer context-panel widgets not registered on this deployment');

  try {
    // SCENARIO 1 — 3D Structure pane for a Molecule3D cell; assert mount via .bsv-container-info-panel.
    await softStep('Scenario 1 — 3D Structure pane mounts for Molecule3D cell; container class .bsv-container-info-panel present', async () => {
      await page.evaluate(() => {
        // Force the context panel visible: the baseline sets simpleMode=true which hides it,
        // and the accordion panes ([name="pane-*"]) only render into a shown context panel.
        grok.shell.windows.showProperties = true;
        const cell = grok.shell.tv.dataFrame.cell(0, 'structure');
        grok.shell.o = DG.SemanticValue.fromTableCell(cell);
      });

      await page.locator('[name="pane-3D-Structure"]').waitFor({timeout: 30_000});

      const res = await page.evaluate(() => ({
        pane3DPresent: !!document.querySelector('[name="pane-3D-Structure"]'),
        header3DPresent: !!document.querySelector('[name="div-section--3D-Structure"]'),
        header3DAriaBefore: document.querySelector('[name="div-section--3D-Structure"]')?.getAttribute('aria-expanded') ?? null,
      }));

      expect(res.pane3DPresent).toBe(true);
      expect(res.header3DPresent).toBe(true);
      // Collapsed by default.
      expect(res.header3DAriaBefore).toBe('false');

      // Expand the accordion section header.
      await page.locator('[name="div-section--3D-Structure"]').click({timeout: 10_000});
      // structure3D() appends .bsv-container-info-panel once the Mol* viewer mounts — poll for it.
      await page.locator('[name="pane-3D-Structure"] .d4-accordion-pane-content .bsv-container-info-panel')
        .waitFor({timeout: 15_000});

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

      expect(after.header3DAriaAfter).toBe('true');
      expect(after.contentChildren).toBeGreaterThan(0);
      // structure3D() adds .bsv-container-info-panel after render — the mount-success marker.
      expect(after.hasBsvContainer).toBe(true);
      expect(after.dataSource).toBe('Biostructure Viewer:3D Structure');
    });

    await softStep('Scenario 1 — re-mount: structure3D(SemanticValue) yields a fresh DG.Widget that renders each distinct cell value (SR-03)', async () => {
      await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        const fns = DG.Func.find({name: 'structure3D', package: 'BiostructureViewer'});
        const fn = fns && fns[0];
        const cell0 = df.cell(0, 'structure'); // 1bdq
        const cell1 = df.cell(1, 'structure'); // protease — genuinely different structure
        const w0 = await fn?.apply({molecule: DG.SemanticValue.fromTableCell(cell0)});
        const w1 = await fn?.apply({molecule: DG.SemanticValue.fromTableCell(cell1)});
        (window as any).__remount = {fn, w0, w1, val0: cell0.value, val1: cell1.value};
      });

      // Each re-mounted widget must actually render its own Mol* viewer (.bsv-container-info-panel),
      // not merely return a distinct object — the same mount marker Scenario 1 already asserts.
      await page.waitForFunction(() => {
        const r = (window as any).__remount;
        return !!r?.w0?.root?.querySelector('.bsv-container-info-panel') &&
          !!r?.w1?.root?.querySelector('.bsv-container-info-panel');
      }, null, {timeout: 30_000});

      const res = await page.evaluate(() => {
        const r = (window as any).__remount;
        return {
          fnFound: !!r.fn,
          w0Present: !!r.w0,
          w1Present: !!r.w1,
          valuesDistinct: r.val0 !== r.val1,
          rootsDistinct: r.w0.root !== r.w1.root,
          w0Mounted: !!r.w0.root.querySelector('.bsv-container-info-panel'),
          w1Mounted: !!r.w1.root.querySelector('.bsv-container-info-panel'),
        };
      });

      expect(res.fnFound).toBe(true);
      // The two cells carry genuinely different structures — the re-mount is a real value change.
      expect(res.valuesDistinct).toBe(true);
      expect(res.w0Present).toBe(true);
      expect(res.w1Present).toBe(true);
      expect(res.rootsDistinct).toBe(true);
      expect(res.w0Mounted).toBe(true);
      expect(res.w1Mounted).toBe(true);
    });

    // SCENARIO 2 — PDB Information pane on a Molecule3D cell (file-info from the PDB string).
    await softStep('Scenario 2 — PDB Information pane (Molecule3D) renders pdbFileInfoWidget content', async () => {
      await page.evaluate(() => {
        // Force the context panel visible: the baseline sets simpleMode=true which hides it,
        // and the accordion panes ([name="pane-*"]) only render into a shown context panel.
        grok.shell.windows.showProperties = true;
        const cell = grok.shell.tv.dataFrame.cell(0, 'structure');
        grok.shell.o = DG.SemanticValue.fromTableCell(cell);
      });

      await page.locator('[name="pane-PDB-Information"]').waitFor({timeout: 30_000});

      await page.locator('[name="div-section--PDB-Information"]').click({timeout: 10_000});
      // pdb-file-info derives synchronously from the PDB string; poll for its deterministic markers.
      await expect.poll(async () => page.evaluate(() => {
        const c = document.querySelector('[name="pane-PDB-Information"] .d4-accordion-pane-content');
        const t = c?.textContent || '';
        return t.includes('Atom Count') && t.includes('Chains');
      }), {timeout: 8_000, intervals: [250, 500, 1000]}).toBe(true);

      const after = await page.evaluate(() => {
        const pane = document.querySelector('[name="pane-PDB-Information"]');
        const content = pane?.querySelector('.d4-accordion-pane-content');
        return {
          ariaAfter: document.querySelector('[name="div-section--PDB-Information"]')?.getAttribute('aria-expanded') ?? null,
          textLen: content?.textContent?.length ?? 0,
          textSnippet: (content?.textContent || '').slice(0, 200),
          hasGeneralAtomCount: (content?.textContent || '').includes('Atom Count'),
          hasChains: (content?.textContent || '').includes('Chains'),
        };
      });

      expect(after.ariaAfter).toBe('true');
      expect(after.textLen).toBeGreaterThan(50);
      // pdb-file-info derives from the PDB string: "Atom Count" + "Chains" are the deterministic markers.
      expect(after.hasGeneralAtomCount).toBe(true);
      expect(after.hasChains).toBe(true);
    });

    // SCENARIO 3 — PDB Information pane on a PDB_ID cell (async pdbInfoWidget assembling RCSB metadata).
    await softStep('Scenario 3 — PDB Information pane (PDB_ID) renders async pdbInfoWidget metadata', async () => {
      await page.evaluate(() => {
        // Context panel must be shown for the accordion panes to render (see scenario 1).
        grok.shell.windows.showProperties = true;
        const cell = grok.shell.tv.dataFrame.cell(0, 'pdb_id');
        grok.shell.o = DG.SemanticValue.fromTableCell(cell);
      });

      // Both panes register on a PDB_ID SemanticValue.
      await page.waitForFunction(() =>
        !!document.querySelector('[name="pane-PDB-Information"]') &&
        !!document.querySelector('[name="pane-PDB-id-viewer"]'), null, {timeout: 15_000});

      await page.locator('[name="div-section--PDB-Information"]').click({timeout: 10_000});
      // Poll for the async pdbInfoWidget body (RCSB metadata assembly) instead of a blind wait.
      await expect.poll(async () => page.evaluate(() => {
        const c = document.querySelector('[name="pane-PDB-Information"] .d4-accordion-pane-content');
        const t = c?.textContent || '';
        return t.includes('Description') && t.includes('rcsb.org/structure/1QBS') && !t.includes('Could not fetch PDB');
      }), {timeout: 12_000, intervals: [500, 1000, 2000]}).toBe(true);

      const after = await page.evaluate(() => {
        const pane = document.querySelector('[name="pane-PDB-Information"]');
        const content = pane?.querySelector('.d4-accordion-pane-content');
        const text = content?.textContent || '';
        return {
          ariaAfter: document.querySelector('[name="div-section--PDB-Information"]')?.getAttribute('aria-expanded') ?? null,
          textLen: text.length,
          textSnippet: text.slice(0, 200),
          hasDescription: text.includes('Description'),
          hasRcsbUrl: text.includes('rcsb.org/structure/1QBS'),
          noFetchError: !text.includes('Could not fetch PDB'),
        };
      });

      expect(after.ariaAfter).toBe('true');
      expect(after.textLen).toBeGreaterThan(50);
      expect(after.hasDescription).toBe(true);
      expect(after.hasRcsbUrl).toBe(true);
      expect(after.noFetchError).toBe(true);
    });

    // SCENARIO 4 — Protein-Ligand Interactions (ProLIF): three condition-gated registrations.
    //   Path A: Molecule3D + hasNonWaterHetatm; Path B: isAutoDockPose (deferred, registry probe);
    //   Path C: PDB_ID via RCSB fetchProxy.
    await softStep('Scenario 4 Path A — ProLIF panel (Molecule3D + hasNonWaterHetatm) renders pdbInteractionsWidget LigNetwork', async () => {
      const surface = await page.evaluate(async () => {
        // Force the context panel visible: the baseline sets simpleMode=true which hides it,
        // and the accordion panes ([name="pane-*"]) only render into a shown context panel.
        grok.shell.windows.showProperties = true;
        const cell = grok.shell.tv.dataFrame.cell(0, 'structure');
        grok.shell.o = DG.SemanticValue.fromTableCell(cell);
        // Live-confirm the gating predicate on the fixture.
        const pdbContent = await grok.dapi.files.readAsText('System:AppData/BiostructureViewer/samples/1bdq.pdb');
        const hasHet = await grok.functions.call('BiostructureViewer:hasNonWaterHetatm', {molecule: pdbContent});
        const isAuto = await grok.functions.call('BiostructureViewer:isAutoDockPose', {molecule: pdbContent});
        return {hasNonWaterHetatm: hasHet, isAutoDockPose: isAuto};
      });

      await page.waitForFunction(() =>
        !!document.querySelector('[name="pane-Protein-Ligand-Interactions"]'), null, {timeout: 15_000});

      // Path A is the resolution path for this fixture (hasNonWaterHetatm=true; isAutoDockPose=false).
      expect(surface.hasNonWaterHetatm).toBe(true);
      expect(surface.isAutoDockPose).toBe(false);

      await page.locator('[name="div-section--Protein-Ligand-Interactions"]').click({timeout: 10_000});
      // Poll for the async LigNetwork assembly (fixture ligand IM1) instead of a blind wait.
      await expect.poll(async () => page.evaluate(() => {
        const c = document.querySelector('[name="pane-Protein-Ligand-Interactions"] .d4-accordion-pane-content');
        const t = c?.textContent || '';
        return /\d+\s+ligands?\s+found/.test(t) && t.includes('IM1');
      }), {timeout: 12_000, intervals: [500, 1000, 2000]}).toBe(true);

      const after = await page.evaluate(() => {
        const pane = document.querySelector('[name="pane-Protein-Ligand-Interactions"]');
        const content = pane?.querySelector('.d4-accordion-pane-content');
        const text = content?.textContent || '';
        return {
          ariaAfter: document.querySelector('[name="div-section--Protein-Ligand-Interactions"]')?.getAttribute('aria-expanded') ?? null,
          textLen: text.length,
          textSnippet: text.slice(0, 200),
          dataSource: content?.querySelector('[data-source]')?.getAttribute('data-source') ?? null,
          hasLigandsFound: /\d+\s+ligands?\s+found/.test(text),
          hasIM1: text.includes('IM1'),
        };
      });

      expect(after.ariaAfter).toBe('true');
      expect(after.textLen).toBeGreaterThan(50);
      expect(after.hasLigandsFound).toBe(true);
      // 1bdq's ligand IM1 — fixture-specific marker.
      expect(after.hasIM1).toBe(true);
      expect(after.dataSource).toBe('Biostructure Viewer:Protein-Ligand Interactions');
    });

    await softStep('Scenario 4 Path B — DG.Func registration probe for dockingInteractionsWidget (SR-02: AutoDock pose fixture not available)', async () => {
      // Docking-flavour ProLIF asserted at the registry level: role=panel, friendlyName, and the
      // isAutoDockPose condition predicate (no AutoDock pose fixture available).
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

      expect(res.registered).toBe(true);
      expect(res.inputCount).toBeGreaterThan(0);
      expect(res.inputSemTypes).toContain('Molecule3D');
      expect(res.isPanel).toBe(true);
      expect(res.friendlyName).toBe('Protein-Ligand Interactions');
      // The condition predicate differentiates the three same-named registrations.
      expect(res.condition).toBe('BiostructureViewer:isAutoDockPose(molecule)');
    });

    await softStep('Scenario 4 Path C — ProLIF panel (PDB_ID) renders pdbIdInteractionsWidget via RCSB fetchProxy', async () => {
      await page.evaluate(() => {
        // Context panel must be shown for the accordion panes to render (see scenario 1).
        grok.shell.windows.showProperties = true;
        const cell = grok.shell.tv.dataFrame.cell(0, 'pdb_id');
        grok.shell.o = DG.SemanticValue.fromTableCell(cell);
      });

      await page.waitForFunction(() =>
        !!document.querySelector('[name="pane-Protein-Ligand-Interactions"]'), null, {timeout: 15_000});

      // Expand the PDB_ID-flavour ProLIF pane; poll for the RCSB fetchProxy + LigNetwork result.
      await page.locator('[name="div-section--Protein-Ligand-Interactions"]').click({timeout: 10_000});
      await expect.poll(async () => page.evaluate(() => {
        const c = document.querySelector('[name="pane-Protein-Ligand-Interactions"] .d4-accordion-pane-content');
        const t = c?.textContent || '';
        return /\d+\s+ligands?\s+found/.test(t) && t.includes('DMP') && !t.includes('Could not fetch PDB');
      }), {timeout: 20_000, intervals: [500, 1000, 2000]}).toBe(true);

      const after = await page.evaluate(() => {
        const pane = document.querySelector('[name="pane-Protein-Ligand-Interactions"]');
        const content = pane?.querySelector('.d4-accordion-pane-content');
        const text = content?.textContent || '';
        return {
          ariaAfter: document.querySelector('[name="div-section--Protein-Ligand-Interactions"]')?.getAttribute('aria-expanded') ?? null,
          textLen: text.length,
          textSnippet: text.slice(0, 200),
          dataSource: content?.querySelector('[data-source]')?.getAttribute('data-source') ?? null,
          hasLigandsFound: /\d+\s+ligands?\s+found/.test(text),
          hasDmp: text.includes('DMP'),
          noFetchError: !text.includes('Could not fetch PDB'),
        };
      });

      expect(after.ariaAfter).toBe('true');
      expect(after.hasLigandsFound).toBe(true);
      // 1QBS's known ligand DMP — proves the RCSB fetch succeeded.
      expect(after.hasDmp).toBe(true);
      expect(after.dataSource).toBe('Biostructure Viewer:Protein-Ligand Interactions');
      expect(after.noFetchError).toBe(true);
    });

    // SCENARIO 5 — Link With Molecule Column: mol3dAtomPickerLinkWidget returns a DG.Widget with the
    //   SMILES-column ChoiceInput (asserted via the direct JS-API call surface).
    await softStep('Scenario 5 — mol3dAtomPickerLinkWidget(mol3DCol) returns DG.Widget with SMILES-column ChoiceInput', async () => {
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const df = tv.dataFrame;

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

      expect(res.callErr).toBe(null);
      expect(res.widgetOk).toBe(true);

      expect(res.widgetDataSource).toBe('Biostructure Viewer:Link With Molecule Column');

      // SMILES-column picker scoped to Molecule semType columns (options include 'ligand').
      expect(res.smilesPickerPresent).toBe(true);
      expect(res.smilesPickerOptions).not.toBe(null);
      expect((res.smilesPickerOptions || []).includes('ligand')).toBe(true);

      // Function-registry contract.
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
