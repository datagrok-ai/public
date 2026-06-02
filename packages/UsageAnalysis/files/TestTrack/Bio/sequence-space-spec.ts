/* ---
sub_features_covered:
  - bio.analyze.sequence-space
  - bio.analyze.sequence-space.top-menu
  - bio.analyze.sequence-space.editor
  - bio.analyze.sequence-space.transform
--- */
//   related_bugs: [GROK-18616, GROK-19928] (cross-cutting invariants delegated
//     to chain bug_focused_candidates: bio-grok-18616-spec.ts /
//     bio-grok-19928-spec.ts per scenario Notes)
// SCOPE_REDUCTIONS honoured from scenario frontmatter:
//   SR-01 (A-CONT-01) — "arbitrary Similarity/Method edit set" not defined in
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
test.use(specTestOptions);
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/tests/filter_MSA.csv'},
];
for (const ds of datasets) {
  test(`Bio Sequence Space on ${ds.name}`, async ({page}) => {
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
    // Bio top-menu + init-completion readiness (cold-start stabilization,
    // mirrors the analyze-spec.ts / sequence-activity-cliffs-spec.ts cycle-2
    // retry pattern). Two-layer guard:
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
    // Per-leaf function-registration probe (cycle-2 retry refinement from
    // analyze-spec.ts). The init probe above guarantees init COMPLETION; it
    // does NOT guarantee the Bio:sequenceSpaceTopMenu leaf is findable in the
    // function registry yet. Bounded poll until findable, short defensive
    // settle if the candidate names diverged across a Bio version.
    await page.evaluate(async () => {
      const candidates = ['Bio:sequenceSpaceTopMenu', 'Bio:sequenceSpace'];
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
      await new Promise((r) => setTimeout(r, 1500));
    });
    await page.waitForTimeout(2000);
    // Scenario 1, Step 2 — Open Bio > Analyze > Sequence Space (defaults run).
    // Atlas: bio.analyze.sequence-space.top-menu (package.ts#L740),
    // bio.analyze.sequence-space.editor (SequenceSpaceEditor — package.ts#L237).
    //
    // Click pattern per bio.md ("Click pattern (MCP-validated, mirrors chem.md)"):
    // click [name="div-Bio"] → 400ms → mouseover [name="div-Bio---Analyze"] →
    // 300ms → click leaf. Hover-not-click on the group is required to surface
    // the leaves. Top-menu submenu opening is hover-driven; CDP synthetic
    // clicks alone do not trigger the Dart onMouseEnter listener that calls
    // _showSubMenu.
    await softStep(`${ds.name}: Open Bio > Analyze > Sequence Space (defaults)`, async () => {
      // Source text says "Bio > Search > Sequence Space" but atlas and live
      // platform have it at "Bio > Analyze > Sequence Space..."; the migrated
      // scenario uses the atlas-canonical path
      // (bio.md L94 source-text adjudication).
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Analyze"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Sequence-Space..."]') as HTMLElement).click();
      });
      // Cold-start tolerance: 60s tolerates the empirical cold ceiling per
      // analyze-spec.ts cycle-2 retry evidence (Bio package function dispatch
      // can exceed 15s on a first-ever cold MSA boot).
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
      const title = await page.locator('.d4-dialog .d4-dialog-title').textContent();
      expect(title?.trim()).toBe('Sequence Space');
    });
    // Scenario 1, Step 3 — Click OK to run with default parameters.
    // Atlas: bio.analyze.sequence-space.transform (sequenceSpaceTransform —
    // package.ts#L789, via bio.engines.preprocess-encode).
    //
    // Scenario 1, Step 4 — Verify ScatterPlot viewer titled "Embeddings"
    // docks and embedding columns (Embed_X_1, Embed_Y_1) are appended.
    // Cluster (DBSCAN) column is conditionally appended when the dialog's
    // `Cluster embeddings` checkbox was left on by default.
    await softStep(`${ds.name}: Run with default parameters — ScatterPlot + embeddings appended`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      // Compound invariant: column count grows (Embed_X_1 + Embed_Y_1 [+
      // Cluster (DBSCAN)]) AND a Scatter plot viewer mounts on the active
      // TableView. The 240s budget tolerates the embedding compute on a cold
      // first-run boot.
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 240_000});
      const result = await page.evaluate(() => ({
        hasScatter: Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        scatterCount: Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length,
        cols: grok.shell.tv.dataFrame.columns.length,
        colNames: Array.from({length: grok.shell.tv.dataFrame.columns.length},
          (_, i) => grok.shell.tv.dataFrame.columns.byIndex(i).name),
      }));
      expect(result.hasScatter).toBe(true);
      expect(result.scatterCount).toBe(1);
      expect(result.cols).toBeGreaterThan(baseCols);
      // Embedding column-name invariant (bio.md L133 post-OK invariants:
      // default-named Embed_X_1, Embed_Y_1). Structural assertion only; the
      // numeric distribution of the embedding is deferred per SR-01-style
      // boundary for defaults too (atlas does not pin a canonical embedding
      // distribution for these fixtures).
      const hasEmbedX = result.colNames.some((n: string) => /^Embed_X_\d+$/.test(n));
      const hasEmbedY = result.colNames.some((n: string) => /^Embed_Y_\d+$/.test(n));
      expect(hasEmbedX).toBe(true);
      expect(hasEmbedY).toBe(true);
    });
    // Close the scatter plot from the first run so the second run's
    // structural invariant (distinct viewer mount) is clearly observable.
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });
    // Scenario 2, Step 5 — Re-open Bio > Analyze > Sequence Space.
    // Editor input re-binding contract (bio.analyze.sequence-space.editor).
    await softStep(`${ds.name}: Re-open Bio > Analyze > Sequence Space`, async () => {
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Analyze"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Sequence-Space..."]') as HTMLElement).click();
      });
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
    });
    // Scenario 2, Step 6 — Change Similarity + Method (edit-then-run flow).
    // Per SR-01: source text says "arbitrarily" — atlas does not pin a
    // canonical edit set. The concrete picks below (UMAP → t-SNE, Hamming →
    // Levenshtein) come from the prior run log (sequence-space-run.md); they
    // exercise the editor input re-binding contract (inputs must accept new
    // selections via the <select> widget surface) but are NOT canonical for
    // correctness.
    //
    // Per bio.md L127-L128: Method options [UMAP, t-SNE], Similarity options
    // [Hamming, Levenshtein, Monomer chemical distance, Needlemann-Wunsch].
    // Selectors addressed by [name="input-Method"] / [name="input-Similarity"]
    // SELECTs directly (mirror sequence-activity-cliffs-spec.ts pattern); fall
    // back to [name="input-host-*"] wrapper if a build wraps the SELECT.
    await softStep(`${ds.name}: Change Method to t-SNE and Similarity to Levenshtein`, async () => {
      await page.evaluate(async () => {
        const dlg = document.querySelector('.d4-dialog')!;
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
    // Structural invariant survives SR-01 deferral: a second "Embeddings"
    // Scatter plot must dock (distinct from the first run's viewer — the
    // first was closed above). Total Scatter-plot viewer count = 1 after the
    // edited-parameter run; dataframe grows by a fresh embedding pair
    // (Embed_X_2, Embed_Y_2 [+ Cluster (DBSCAN) #2 if clustering on]). The
    // correctness assertion on the edited-parameter embedding distribution
    // is deferred per SR-01.
    await softStep(`${ds.name}: Run with edited parameters — second embedding result docks`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 240_000});
      const result = await page.evaluate(() => ({
        hasScatter: Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        cols: grok.shell.tv.dataFrame.columns.length,
        colNames: Array.from({length: grok.shell.tv.dataFrame.columns.length},
          (_, i) => grok.shell.tv.dataFrame.columns.byIndex(i).name),
      }));
      expect(result.hasScatter).toBe(true);
      // Structural-only assertion: a Scatter plot is present and the dataframe
      // grew by an embedding pair. The "edited result differs meaningfully
      // from defaults" assertion is deferred per SR-01.
      expect(result.cols).toBeGreaterThan(baseCols);
      // The suffix-incremented embedding pair invariant (bio.md L133): the
      // second run produces Embed_X_2 / Embed_Y_2 alongside the first run's
      // Embed_X_1 / Embed_Y_1 (the first run's columns remain because the
      // first run's viewer was closed but the column append is dataframe-
      // level and persists).
      const embedXCount = result.colNames.filter((n: string) => /^Embed_X_\d+$/.test(n)).length;
      const embedYCount = result.colNames.filter((n: string) => /^Embed_Y_\d+$/.test(n)).length;
      expect(embedXCount).toBeGreaterThanOrEqual(2);
      expect(embedYCount).toBeGreaterThanOrEqual(2);
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
