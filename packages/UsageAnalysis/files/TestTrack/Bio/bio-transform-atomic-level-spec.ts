/* ---
sub_features_covered:
  - bio.transform.to-atomic-level
  - bio.transform.to-atomic-level.action
  - bio.transform.to-atomic-level.single
  - bio.transform.to-atomic-level.api
  - bio.transform.helm-to-mol
  - bio.transform.molecules-to-helm
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: regression — source-matrix-like
//     runner with broad JS-API permission + ≥1 DOM-driving call required;
//     same body shape as sibling convert-spec.ts on coverage_type: regression)
//   sub_features_covered: [bio.transform.to-atomic-level, .action, .single,
//     .api, bio.transform.helm-to-mol, bio.transform.molecules-to-helm]
//   ui_coverage_responsibility: absent (delegated_to: null) — no specialty
//     ui-smoke ownership; Bio top-menu + dialog OK is the DOM-driving slice;
//     API wrappers + GROK-15176 isotope-flag regression guard run as JS-API
//     assertions inside the same Playwright runner per scenario Notes.
//   related_bugs: [GROK-15176] — bug-focused regression guard inline
//     (Scenario 5); molfile must not carry MASS=1 isotope flags on heavy
//     atoms; ties to cross-feature interaction bio.x.bio-to-chem-via-atomic-level.
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/bio.yaml#sub_features[bio.transform.to-atomic-level]
//     interactions[0] = "Bio | Transform | To Atomic Level..."
//     source = public/packages/Bio/src/package.ts#L863
//   feature-atlas/bio.yaml#sub_features[bio.transform.to-atomic-level.action]
//     interactions[0] = "column action: to atomic level"
//     source = public/packages/Bio/src/package.ts#L886 (toAtomicLevelAction)
//   feature-atlas/bio.yaml#sub_features[bio.transform.to-atomic-level.single]
//     interactions[0] = grok.functions.call('Bio:toAtomicLevelSingleSeq')
//     source = public/packages/Bio/src/package.ts#L911
//   feature-atlas/bio.yaml#sub_features[bio.transform.to-atomic-level.api]
//     interactions[0] = grok.functions.call('Bio:seq2atomic')
//     source = public/packages/Bio/src/package.ts#L1605
//   feature-atlas/bio.yaml#sub_features[bio.transform.helm-to-mol]
//     interactions[0] = HelmToMolfileConverter pipeline
//     source = public/packages/Bio/src/package.ts#L1684 (helper-converter wrapper)
//   feature-atlas/bio.yaml#sub_features[bio.transform.molecules-to-helm]
//     interactions[0] = "Bio | Transform | Molecules to HELM..."
//     source = public/packages/Bio/src/package.ts#L824 (moleculesToHelmTopMenu)
//   feature-atlas/bio.yaml#critical_paths[bio.cp.to-atomic-level] (p1)
//   feature-atlas/bio.yaml#interactions[bio.x.bio-to-chem-via-atomic-level]
//     coverage_type: regression, related_bugs: [GROK-15176]
//
// Bug-library cross-reference:
//   bug-library/bio.yaml :: GROK-15176 — to-atomic-level molfile must not
//   contain heavy-atom oxygens with isotope=1 (MASS=1 flag); downstream
//   PubChem standardization (Chem panel) rejects illegal isotopes,
//   breaking the bio.x.bio-to-chem-via-atomic-level renderer contract.
//
// All selectors used here are class-1 — present in
// grok-browser/references/bio.md (lines 379-415 for Transform top-menu
// dialogs; lines 499-521 for Macromolecule column Context Pane actions).
// No class-2 selector recon-notes required.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Fixtures match scenario Setup section — the two canonical Bio test
// fixtures cover the HELM and FASTA branches of toAtomicLevel.
const datasets = [
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv', units: 'helm'},
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv', units: 'fasta'},
];

// Heavy-atom isotope-flag scanner (GROK-15176 invariant — see scenario
// Scenario 5). V3000 atom lines look like:
//   "M  V30 <idx> <symbol> x y z 0 ..."
// Isotope is encoded via the "MASS=<n>" key on the atom line. The bug
// surfaced as MASS=1 on heavy-atom oxygens (illegal isotope of O), which
// PubChem rejects. We exclude H / D explicitly: a "1" mass on hydrogen
// is the natural isotope and not the regression signal.
function heavyAtomIsotopeFlags(mol: string): string[] {
  if (!mol) return [];
  return mol.split(/\r?\n/)
    .filter((l) => /^M  V30 \d+ [A-Za-z]+ /.test(l))
    .filter((l) => / MASS=1\b/.test(l))
    .filter((l) => !/ [HD] /.test(l));
}

for (const ds of datasets) {
  test(`Bio Transform To Atomic Level + round-trip on ${ds.name}`, async ({page}) => {
    test.setTimeout(600_000);
    stepErrors.length = 0;

    await loginToDatagrok(page);

    // Setup phase: open dataset, wait for Macromolecule semType detection
    // + Bio package init. Mirrors convert-spec.ts cold-start tolerance —
    // the Bio top-menu doesn't appear until Bio package functions are
    // registered against the active Macromolecule TableView.
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
    }, ds.path);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

    // Bio top-menu readiness poll — two-layer guard (same as convert-spec.ts
    // and analyze-spec.ts cold-start stabilization).
    //
    // Layer 1 — DOM-level top-menu visibility.
    // Layer 2 — init-order serialization probe via Bio:getSeqHelper /
    // getMonomerLibHelper / getBioLib (atlas bio.cp.bio-service-surface-init).
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await page.evaluate(async () => {
      const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
      }
      await new Promise((r) => setTimeout(r, 3000));
    });

    // Sanity: Macromolecule column is tagged with the expected units
    // (atlas bio.detector). This anchors the scenario's Setup expectation
    // before any Transform action runs.
    await softStep(`${ds.name}: Macromolecule column has units=${ds.units}`, async () => {
      const info: {hasMacro: boolean, units: string | null} = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        return {hasMacro: !!macro, units: macro?.meta?.units ?? null};
      });
      expect(info.hasMacro).toBe(true);
      expect(info.units).toBe(ds.units);
    });

    // Scenario 1 (HELM dataset) / Scenario 3 (FASTA dataset) — Top-menu
    // To Atomic Level... (atlas bio.transform.to-atomic-level,
    // package.ts#L863). On HELM input the pipeline routes through
    // HelmToMolfileConverter (atlas bio.transform.helm-to-mol,
    // package.ts#L1684); on FASTA it takes the linear _toAtomicLevel
    // branch. Column-level contract is identical: a Molecule/molblock
    // column appears whose first non-null V3K cell starts with the
    // "M  V30 BEGIN CTAB" header.
    await softStep(`${ds.name}: Bio | Transform | To Atomic Level adds a Molecule molblock column`, async () => {
      const beforeMolCount: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Molecule').length;
      });
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Transform"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector(
          '[name="div-Bio---Transform---To-Atomic-Level..."]') as HTMLElement).click();
      });
      await page.locator('[name="dialog-To-Atomic-Level"]').waitFor({timeout: 60_000});
      // Dialog prefills the active Macromolecule column. Accept defaults
      // (Non-linear=true, Highlight=false per package.ts#L863 defaults).
      await page.locator('[name="dialog-To-Atomic-Level"] [name="button-OK"]').click();
      // Atomic-level conversion on the 540-row HELM dataset clocks
      // ~20-30s per convert-run.md 2026-04-23; 120s ceiling is defensive.
      await page.waitForFunction((base) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Molecule').length > base;
      }, beforeMolCount, {timeout: 120_000});

      // Assert column contract: semType=Molecule, units=molblock,
      // first non-null cell is a V3K molfile.
      const info: {hasMol: boolean, units: string | null, name: string | null, firstCell: string | null} =
        await page.evaluate(() => {
          const df = grok.shell.tv.dataFrame;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const mol: any = cols.reverse().find((c: any) => c.semType === 'Molecule');
          let firstCell: string | null = null;
          if (mol) {
            for (let r = 0; r < Math.min(df.rowCount, 50); r++) {
              const v = mol.get(r);
              if (v != null && String(v).length > 0) { firstCell = String(v); break; }
            }
          }
          return {hasMol: !!mol, units: mol?.meta?.units ?? null, name: mol?.name ?? null, firstCell};
        });
      expect(info.hasMol).toBe(true);
      expect(info.units).toBe('molblock');
      // GROK-15176-related shape check: produced cell is a V3K molfile.
      expect(info.firstCell).not.toBeNull();
      expect(info.firstCell!).toMatch(/M\s+V30\s+BEGIN\s+CTAB/);

      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-To-Atomic-Level"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });

    // Scenario 2 — Column-header cell-action "to atomic level"
    // (atlas bio.transform.to-atomic-level.action, package.ts#L886).
    //
    // Per bio.md "Macromolecule cell context menu — hit-test caveat"
    // (lines 531-548): synthetic dispatchEvent on the canvas reaches
    // only the Grid title-bar context menu, NOT the cell or
    // column-header context menu — the hit-test code path needs a real
    // OS mouse event. The `.action` is also surfaced via the column
    // Context Panel (lines 499-521): focus the column with
    // `grok.shell.o = df.col(...)`, expand the Actions accordion pane,
    // click the `label.d4-link-action` whose text starts with
    // "To Atomic Level". This DOM-clicking path drives the same
    // toAtomicLevelAction function via the platform's action registry
    // (the action's func.prepare().edit() opens the same
    // dialog-To-Atomic-Level surface).
    await softStep(`${ds.name}: Context Panel column-action "To Atomic Level" opens the same dialog`, async () => {
      // Set context object = the Macromolecule column to populate the
      // column-scoped Context Panel.
      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        if (macro) grok.shell.o = macro;
      });
      // Give the Context Panel a moment to render the accordion panes.
      await page.waitForTimeout(800);

      // Expand the column-scoped Actions pane. There are two `Actions`
      // accordion panes in the right column (toolbox-section + column-scoped);
      // we click them all — expanding the toolbox one is a harmless no-op.
      await page.evaluate(() => {
        for (const h of Array.from(document.querySelectorAll('.d4-accordion-pane-header'))) {
          if ((h as HTMLElement).textContent === 'Actions')
            h.dispatchEvent(new MouseEvent('click', {bubbles: true}));
        }
      });
      await page.waitForTimeout(800);

      // Click the "To Atomic Level..." link action. Once expanded the
      // label is present in the DOM; before expand it's `null`. Find by
      // text-starts-with per bio.md line 517.
      const clicked: boolean = await page.evaluate(() => {
        const labels = Array.from(document.querySelectorAll('label.d4-link-action'));
        const link = labels.find((l) => (l as HTMLElement).textContent?.trim()
          .toLowerCase().startsWith('to atomic level'));
        if (link) {
          (link as HTMLElement).click();
          return true;
        }
        return false;
      });
      expect(clicked).toBe(true);

      // The action invokes toAtomicLevelAction → func.prepare(...).edit()
      // which opens the same dialog-To-Atomic-Level surface as the
      // top-menu path. Close the dialog without re-running the conversion
      // (the previous step already exercised the conversion path).
      await page.locator('[name="dialog-To-Atomic-Level"]').waitFor({timeout: 30_000});
      const cancelBtn = page.locator('[name="dialog-To-Atomic-Level"] [name="button-CANCEL"]');
      if (await cancelBtn.count() > 0)
        await cancelBtn.first().click();
      else
        await page.keyboard.press('Escape');
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-To-Atomic-Level"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });

    // Scenario 4 — Molecules to HELM round-trip (atlas
    // bio.transform.molecules-to-helm, package.ts#L824). With the
    // molfile column from the top-menu step present on the active
    // dataframe, open Bio | Transform | Molecules to HELM... — the
    // dialog appears prefilled with the molfile column; OK runs the
    // Python-script + monomer-library-matching converter that produces
    // a HELM column from molecule input. The round-trip exercises the
    // bio.x.bio-to-chem-via-atomic-level cross-package contract end-to-end.
    await softStep(`${ds.name}: Bio | Transform | Molecules to HELM round-trip adds a Macromolecule column`, async () => {
      const beforeMacroCount: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Macromolecule').length;
      });
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Transform"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector(
          '[name="div-Bio---Transform---Molecules-to-HELM..."]') as HTMLElement).click();
      });
      await page.locator('[name="dialog-Molecules-to-HELM"]').waitFor({timeout: 60_000});
      await page.locator('[name="dialog-Molecules-to-HELM"] [name="button-OK"]').click();
      // Python-script + monomer-library-matching can take 30-60s on a
      // 540-row dataset; 180s ceiling is defensive (matches the
      // 2x molfile conversion-time-equivalent observed in convert-run.md).
      await page.waitForFunction((base) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Macromolecule').length > base;
      }, beforeMacroCount, {timeout: 180_000});

      // Assert: a new Macromolecule column appears. The moleculesToHelm
      // implementation sets meta.units='helm' + cell.renderer='helm' on
      // the regenerated-sequence column per package.ts#L840-844.
      const info: {newMacroCount: number, lastUnits: string | null, lastName: string | null} =
        await page.evaluate(() => {
          const df = grok.shell.tv.dataFrame;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const macros = cols.filter((c: any) => c.semType === 'Macromolecule');
          const last: any = macros[macros.length - 1];
          return {
            newMacroCount: macros.length,
            lastUnits: last?.meta?.units ?? null,
            lastName: last?.name ?? null,
          };
        });
      expect(info.newMacroCount).toBeGreaterThan(beforeMacroCount);
      expect(info.lastUnits).toBe('helm');

      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-Molecules-to-HELM"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });

    // Scenario 5 — Public API wrappers + GROK-15176 isotope-flag
    // regression guard. Calls the two public wrappers directly via
    // grok.functions.call (no dataframe round-trip needed) and asserts
    // the molfile-validity invariant GROK-15176 documents: no heavy
    // atom may carry MASS=1.
    //
    // Bug invariant (GROK-15176): Bio's to-atomic-level produces V3000
    // molfiles whose heavy-atom oxygens carry illegal isotope=1 flags
    // that downstream PubChem standardization (Chem identifier lookup
    // panel) rejects, breaking bio.x.bio-to-chem-via-atomic-level.
    await softStep(`${ds.name}: API wrappers seq2atomic + toAtomicLevelSingleSeq produce V3K molfiles (GROK-15176 isotope guard)`, async () => {
      const out: {molLinear: string | null, molNonLinear: string | null, linearErr: string | null, nonLinearErr: string | null} =
        await page.evaluate(async () => {
          const result: any = {molLinear: null, molNonLinear: null, linearErr: null, nonLinearErr: null};
          try {
            const r = await (grok as any).functions.call('Bio:toAtomicLevelSingleSeq', {sequence: 'ACDEFGHIK'});
            result.molLinear = typeof r === 'string' ? r : (r?.molfile ?? r?.result ?? null);
          } catch (e: any) {
            result.linearErr = String(e?.message ?? e);
          }
          try {
            const r = await (grok as any).functions.call('Bio:seq2atomic', {
              seq: 'PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0', nonlinear: true});
            result.molNonLinear = typeof r === 'string' ? r : (r?.molfile ?? r?.result ?? null);
          } catch (e: any) {
            result.nonLinearErr = String(e?.message ?? e);
          }
          return result;
        });

      // API wrappers MUST be callable and return non-empty string molfiles.
      expect(out.linearErr).toBeNull();
      expect(out.nonLinearErr).toBeNull();
      expect(out.molLinear).not.toBeNull();
      expect(out.molNonLinear).not.toBeNull();
      expect(out.molLinear!.length).toBeGreaterThan(0);
      expect(out.molNonLinear!.length).toBeGreaterThan(0);

      // V3K shape: the molfile contains V3000 markers (header line
      // includes "V3000"; CTAB block opens with "M  V30 BEGIN CTAB").
      expect(out.molLinear!).toMatch(/V3000/);
      expect(out.molNonLinear!).toMatch(/V3000/);
      expect(out.molLinear!).toMatch(/M\s+V30\s+BEGIN\s+CTAB/);
      expect(out.molNonLinear!).toMatch(/M\s+V30\s+BEGIN\s+CTAB/);

      // GROK-15176 invariant: no heavy atom carries MASS=1.
      const offendersLinear = heavyAtomIsotopeFlags(out.molLinear!);
      const offendersNonLinear = heavyAtomIsotopeFlags(out.molNonLinear!);

      // The bug-library entry GROK-15176 names the expected behaviour:
      // a clean V3K molfile that PubChem standardization accepts.
      // Any heavy-atom MASS=1 flag means the bug has regressed.
      expect(offendersLinear, `GROK-15176 regression on linear wrapper: ${offendersLinear.length} heavy-atom MASS=1 flag(s). Offending atom lines:\n${offendersLinear.join('\n')}`).toEqual([]);
      expect(offendersNonLinear, `GROK-15176 regression on non-linear wrapper: ${offendersNonLinear.length} heavy-atom MASS=1 flag(s). Offending atom lines:\n${offendersNonLinear.join('\n')}`).toEqual([]);
    });

    if (stepErrors.length > 0) {
      const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
    }
  });
}
