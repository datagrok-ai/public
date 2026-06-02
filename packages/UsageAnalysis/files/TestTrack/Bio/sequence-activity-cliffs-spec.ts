/* ---
sub_features_covered:
  - bio.analyze.activity-cliffs
  - bio.analyze.activity-cliffs.top-menu
  - bio.analyze.activity-cliffs.editor
  - bio.analyze.activity-cliffs.transform
  - bio.analyze.activity-cliffs.init
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: regression — deeper Activity-Cliffs-specific scenario)
//   sub_features_covered: [bio.analyze.activity-cliffs, .top-menu, .editor,
//     .transform, .init]
//   ui_coverage_responsibility: absent (delegated_to: null)
//   related_bugs: [GROK-19150, GROK-19928, GROK-16111] (cross-cutting invariants
//     delegated to chain bug_focused_candidates: bio-grok-19150-spec.ts /
//     bio-grok-19928-spec.ts; GROK-16111 empty-input is atlas-cross-cutting
//     bio.x.empty-input-on-row-viewers covered separately).
//   produced_from: migrated
//   coverage_type: regression
//
// Atlas provenance: scenario realises atlas critical path bio.cp.activity-cliffs
// (p0, derived_from package.ts#L537) across three canonical Macromolecule
// notations (FASTA, HELM, MSA), exercising the multi-subsystem path: editor
// dialog (bio.analyze.activity-cliffs.editor, SeqActivityCliffsEditor —
// package.ts#L268) → embedding compute (bio.analyze.activity-cliffs.transform,
// seqActivityCliffsTransform — package.ts#L658) → ScatterPlot with cliff
// overlay (bio.analyze.activity-cliffs.init, seqActivityCliffsInitFunction —
// package.ts#L625).
//
// SCOPE_REDUCTIONS honoured from scenario frontmatter:
//   SR-01 (A-CONT-01) — "arbitrary Similarity/Method edit set" not defined in
//     atlas. The edit-then-run flow is preserved (dialog re-opens, inputs are
//     editable via the <select> widget surface, a second cliff result docks);
//     the correctness assertion on the edited-parameter SALI distribution is
//     deferred until atlas or operator supplies a concrete edit set. Concrete
//     picks below (UMAP → t-SNE, Hamming → Levenshtein) come from the prior
//     run log (sequence-activity-cliffs-run.md) — they exercise the input
//     re-binding contract but are NOT canonical for correctness assertion.
//
// Sister specs: sequence-space-spec.ts (parallel deeper scenario for the
// Sequence-Space top menu); analyze-spec.ts (umbrella runner covering Activity
// Cliffs alongside Sequence Space and Composition).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Dataset family pinned to `samples/*.csv` (NOT `tests/filter_*.csv`).
//
// Round-1 retry test-bug fix (cycle 2026-06-01-bio-migrate-02): the original
// migrated paths (`tests/filter_FASTA.csv`, `tests/filter_HELM.csv`) point at
// single-column fixtures that carry ONLY the sequence column with no numeric
// "Activity" sibling. The Activity Cliffs dialog defaults its Activities
// input to the first numeric column; with no numeric column present the
// default resolves to null, and on OK the seqActivityCliffsTransform call
// throws `Cannot read properties of null (reading 'name')` at
// package.ts L679 (`activities.name`). MSA passed only because
// `tests/filter_MSA.csv` happens to carry a real "Activity" column; FASTA
// and HELM `tests/filter_*.csv` lack one. Evidence captured from
// test-playwright-output/.../trace.zip under attempts 1–3, all three
// attempts reproducing the same NullError.
//
// `samples/FASTA.csv` carries (Entry, Length, UniProtKB, Sequence, Activity,
// Cluster) — auto-default picks Activity. `samples/HELM.csv` carries
// (HELM, Activity); `samples/MSA.csv` carries (MSA, Activity). All three
// share the (sequence, Activity) shape — the canonical fixture family for
// Activity Cliffs and the family the sibling sequence-space-spec.ts uses
// for the parallel analyze-search axis. The prior run log
// (sequence-activity-cliffs-run.md) also used the samples/ family
// (FASTA.csv 64 rows) and PASSed end-to-end.
//
// Same-paradigm tactical fix per the cheap-checks-via-trace investigation;
// NOT a paradigm pivot.
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/samples/FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/samples/HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/samples/MSA.csv'},
];

for (const ds of datasets) {
  test(`Bio Sequence Activity Cliffs on ${ds.name}`, async ({page}) => {
    test.setTimeout(600_000);
    stepErrors.length = 0;

    await loginToDatagrok(page);

    // Setup phase: open dataset, wait for Macromolecule semType detection +
    // Bio package init (cell renderer + filter registration).
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

    // Bio top-menu + init-completion readiness (cycle-2 retry pattern from
    // analyze-spec.ts — closes the analogous MSA cold-boot flake where the
    // first softStep dialog never materialized within 15s on a truly-cold
    // Bio package init).
    //
    // Layer 1 — DOM-level top-menu visibility: the Bio top-menu entry appears
    // only once the Bio package functions are registered against the active
    // Macromolecule TableView.
    //
    // Layer 2 — init-order serialization probe: per bio.md "Init-order
    // invariant" (line 566), initBio registers _package.completeInit after
    // building the SeqHelper / MonomerLibManager singletons; the runtime
    // serializes grok.functions.call('Bio:<...>') after init. Awaiting
    // Bio:getSeqHelper is a deterministic "Bio package ready" probe — it
    // blocks until init completes, not error.
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await page.evaluate(async () => {
      const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
      }
      await new Promise((r) => setTimeout(r, 3000));
    });

    // Per-leaf function-registration probe (cycle-2 retry refinement). The
    // init probe above guarantees init COMPLETION; it does NOT guarantee that
    // the Bio:activityCliffsTopMenu leaf is findable in the function registry
    // yet. Bounded poll until findable, with a short defensive settle if not.
    await page.evaluate(async () => {
      const candidates = ['Bio:activityCliffsTopMenu', 'Bio:activityCliffs', 'Bio:macromoleculeActivityCliffs'];
      const findAny = (names: string[]): boolean => {
        for (const n of names) {
          try {
            if ((grok as any).functions.find && (grok as any).functions.find(n)) return true;
          } catch { /* try next */ }
        }
        return false;
      };
      const deadline = Date.now() + 15_000;
      while (Date.now() < deadline) {
        if (findAny(candidates)) return;
        await new Promise((r) => setTimeout(r, 300));
      }
      // Defensive settle if no candidate name resolves (function-rename across
      // Bio versions) — the 60s dialog tolerance below is the defensive ceiling.
      await new Promise((r) => setTimeout(r, 1500));
    });
    await page.waitForTimeout(2000);

    // Scenario 1, Step 2 — Open Bio > Analyze > Activity Cliffs... (defaults run).
    // Atlas: bio.analyze.activity-cliffs.top-menu (package.ts#L537),
    // bio.analyze.activity-cliffs.editor (SeqActivityCliffsEditor — package.ts#L268).
    //
    // Click pattern per bio.md ("Click pattern (MCP-validated, mirrors chem.md)"):
    // click [name="div-Bio"] → 400ms → mouseover [name="div-Bio---Analyze"] →
    // 300ms → click leaf. Hover-not-click on the group is required to surface
    // the leaves. Top-menu submenu opening is hover-driven; CDP synthetic
    // clicks alone do not trigger the Dart onMouseEnter listener that calls
    // _showSubMenu (per the prior run-log retrospective).
    await softStep(`${ds.name}: Open Bio > Analyze > Activity Cliffs (defaults)`, async () => {
      // Source text says "Bio > Search > Sequence Activity Cliffs" but atlas
      // and live platform have it at "Bio > Analyze > Activity Cliffs..."
      // (the function is internally named "Sequence Activity Cliffs"). The
      // migrated scenario uses the atlas-canonical path.
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Analyze"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Activity-Cliffs..."]') as HTMLElement).click();
      });
      // Cold-start tolerance: 60s tolerates the empirical cold ceiling
      // (analyze-spec.ts cycle-2 evidence; per bio.md the Bio package
      // function dispatch can take >15s on a first-ever cold MSA boot).
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
      const title = await page.locator('.d4-dialog .d4-dialog-title').textContent();
      expect(title?.trim()).toBe('Activity Cliffs');
    });

    // Scenario 1, Step 3 — Click OK to run with default parameters.
    // Atlas: bio.analyze.activity-cliffs.transform (seqActivityCliffsTransform
    // — package.ts#L658). Stamps seqActivityCliffsParams tag on the table;
    // embedding X/Y columns + SALI scoring column appended (per bio.md
    // post-OK invariants line 157).
    //
    // Scenario 1, Step 4 — Verify ScatterPlot viewer titled "Activity cliffs"
    // docks (atlas bio.analyze.activity-cliffs.init,
    // seqActivityCliffsInitFunction — package.ts#L625) and embedding columns
    // are appended.
    await softStep(`${ds.name}: Run with default parameters — ScatterPlot + embeddings appended`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      // Compound invariant: column count grows (embedding X/Y + SALI) AND a
      // Scatter plot viewer mounts on the active TableView. The 240s budget
      // tolerates the embedding compute on a cold first-run boot.
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 240_000});
      const result = await page.evaluate(() => ({
        hasScatter: Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        scatterCount: Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length,
        cols: grok.shell.tv.dataFrame.columns.length,
        hasParamsTag: !!grok.shell.tv.dataFrame.getTag('seqActivityCliffsParams'),
      }));
      expect(result.hasScatter).toBe(true);
      expect(result.scatterCount).toBe(1);
      // Tag is a transform-side invariant per bio.md (line 157, "atlas
      // bio.analyze.activity-cliffs.transform writes the tag").
      expect(result.hasParamsTag).toBe(true);
    });

    // Close the scatter plot from the first run so the second run's
    // structural invariant (distinct viewer mount) is clearly observable.
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });

    // Scenario 2, Step 5 — Re-open Bio > Analyze > Activity Cliffs.
    // Editor input re-binding contract (bio.analyze.activity-cliffs.editor).
    await softStep(`${ds.name}: Re-open Bio > Analyze > Activity Cliffs`, async () => {
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Analyze"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Activity-Cliffs..."]') as HTMLElement).click();
      });
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
    });

    // Scenario 2, Step 6 — Change Similarity + Method (edit-then-run flow).
    // Per SR-01: source text says "arbitrarily" — atlas does not pin a
    // canonical edit set. The concrete picks below (UMAP → t-SNE, Hamming →
    // Levenshtein) come from the prior run log; they exercise the editor
    // input re-binding contract (inputs must accept new selections via the
    // <select> widget surface) but are NOT canonical for correctness.
    //
    // Per bio.md "Bool toggle UX" + Activity Cliffs editor table (line 152-153):
    // Method options [UMAP, t-SNE, MCL] and Similarity options
    // [Hamming, Levenshtein, Monomer chemical distance, Needlemann-Wunsch]
    // are addressed by [name="input-Method"] / [name="input-Similarity"]
    // SELECTs directly (NOT input-host wrappers).
    await softStep(`${ds.name}: Change Method to t-SNE and Similarity to Levenshtein`, async () => {
      await page.evaluate(async () => {
        const dlg = document.querySelector('.d4-dialog')!;
        // Try direct [name="input-*"] SELECT first (per bio.md table), fall back
        // to [name="input-host-*"] wrapper if a build wraps the SELECT differently.
        const methodSel = (dlg.querySelector('[name="input-Method"]') as HTMLSelectElement)
          ?? (dlg.querySelector('[name="input-host-Method"] select') as HTMLSelectElement);
        const simSel = (dlg.querySelector('[name="input-Similarity"]') as HTMLSelectElement)
          ?? (dlg.querySelector('[name="input-host-Similarity"] select') as HTMLSelectElement);
        methodSel.value = 't-SNE';
        methodSel.dispatchEvent(new Event('input', {bubbles: true}));
        methodSel.dispatchEvent(new Event('change', {bubbles: true}));
        simSel.value = 'Levenshtein';
        simSel.dispatchEvent(new Event('input', {bubbles: true}));
        simSel.dispatchEvent(new Event('change', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
      });
      const verified = await page.evaluate(() => {
        const dlg = document.querySelector('.d4-dialog')!;
        const methodSel = (dlg.querySelector('[name="input-Method"]') as HTMLSelectElement)
          ?? (dlg.querySelector('[name="input-host-Method"] select') as HTMLSelectElement);
        const simSel = (dlg.querySelector('[name="input-Similarity"]') as HTMLSelectElement)
          ?? (dlg.querySelector('[name="input-host-Similarity"] select') as HTMLSelectElement);
        return {method: methodSel.value, sim: simSel.value};
      });
      expect(verified.method).toBe('t-SNE');
      expect(verified.sim).toBe('Levenshtein');
    });

    // Scenario 2, Step 7 — Click OK to run with edited parameters.
    // Structural invariant survives SR-01 deferral: a second Activity Cliffs
    // Scatter plot must dock (distinct from the first run's viewer — the
    // first was closed above). Total Scatter-plot viewer count = 1 at this
    // point (the closed first run); after the edited-parameter run it must
    // be 1 again (the new viewer) with the dataframe columns growing by a
    // fresh embedding+SALI triplet. The correctness assertion on the edited
    // SALI distribution is deferred per SR-01.
    await softStep(`${ds.name}: Run with edited parameters — second cliff result docks`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 240_000});
      const result = await page.evaluate(() => ({
        hasScatter: Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        cols: grok.shell.tv.dataFrame.columns.length,
      }));
      expect(result.hasScatter).toBe(true);
      // Structural-only assertion: a Scatter plot is present and the dataframe
      // grew by an embedding+SALI triplet. The "edited result differs
      // meaningfully from defaults" assertion is deferred per SR-01.
      expect(result.cols).toBeGreaterThan(baseCols);
    });

    // Final cleanup: close non-Grid viewers so subsequent dataset iterations
    // observe a clean TableView.
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });

    if (stepErrors.length > 0) {
      const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
    }
  });
}
