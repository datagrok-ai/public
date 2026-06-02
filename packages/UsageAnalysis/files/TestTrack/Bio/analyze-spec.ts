/* ---
sub_features_covered:
  - bio.analyze.sequence-space
  - bio.analyze.sequence-space.top-menu
  - bio.analyze.sequence-space.editor
  - bio.analyze.activity-cliffs
  - bio.analyze.activity-cliffs.top-menu
  - bio.analyze.activity-cliffs.editor
  - bio.analyze.composition
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: smoke — umbrella integration scenario)
//   sub_features_covered: [bio.analyze.sequence-space, .sequence-space.top-menu,
//     .sequence-space.editor, bio.analyze.activity-cliffs, .activity-cliffs.top-menu,
//     .activity-cliffs.editor, bio.analyze.composition]
//   ui_coverage_responsibility: absent (delegated_to: null)
//   related_bugs: [GROK-18616, GROK-19928, GROK-19150] (umbrella surfaces only —
//     per-bug repro slices are delegated to chain bug_focused_candidates:
//     bio-grok-18616-spec.ts / bio-grok-19928-spec.ts / bio-grok-19150-spec.ts).
//   produced_from: migrated
//   coverage_type: smoke
//
// Atlas provenance: scenario realises the umbrella runner for the three
// Bio | Analyze top-menu functions across three Macromolecule notations
// (FASTA, HELM, MSA). Deep per-function coverage lives in sibling specs:
// sequence-space-spec.ts, sequence-activity-cliffs-spec.ts, composition-analysis-spec.ts.
//
// SCOPE_REDUCTIONS honoured from scenario frontmatter:
//   SR-01 (A-CONT-01) — "arbitrary changed parameters" run deferred (no
//     deterministic edit set in atlas). Only run-with-defaults asserted here.
//   SR-02 (A-CONT-01) — Composition Context-Pane property checklist deferred.
//     Step 5 verifies Context Panel property-editor surface presence only.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/tests/filter_MSA.csv'},
];

for (const ds of datasets) {
  test(`Bio Analyze umbrella on ${ds.name}`, async ({page}) => {
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

    // Bio top-menu readiness poll (cold-start stabilization per Gate-B FLAKY
    // attempt-1 evidence: MSA dataset's first softStep dialog never materialized
    // within 15s on a truly-cold Bio package init). Two-layer guard:
    //
    // Layer 1 — DOM-level top-menu visibility: the Bio top-menu entry appears
    // only once the Bio package functions are registered against the active
    // Macromolecule TableView. Polling for it mirrors the retried-readiness-
    // poll-on-stable-selectors pattern that resolved the analogous
    // mutation-cliffs cold-start flake.
    //
    // Layer 2 — init-order serialization probe: per bio.md "Init-order
    // invariant" (line 566), `initBio` registers `_package.completeInit` after
    // building the SeqHelper / MonomerLibManager singletons; the runtime
    // serializes `grok.functions.call('Bio:<...>')` after init. Awaiting
    // `Bio:getSeqHelper` (bio.md line 560) is a deterministic "Bio package
    // ready" probe — it blocks until init completes, not error. The dialog
    // 60_000 ms timeout below is the defensive ceiling; this probe addresses
    // the same race at its root cause rather than only widening the timeout.
    // This is the addition the cycle-2 retry layers on top of the cycle-1
    // tactical fix (timeout bump 15s → 60s + Bio top-menu poll alone), which
    // Gate-B's MSA-cold-start attempt-1 still raced past.
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await page.evaluate(async () => {
      // Tolerate a missing function name across Bio releases — fall back through
      // a small list of safe init-completion probes before giving up.
      const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
      }
      // Fallback: short settle so spec doesn't immediately race the dialog open
      // if every named probe is unavailable on this Bio build.
      await new Promise((r) => setTimeout(r, 3000));
    });

    // Per-leaf function-registration probe (cycle-2 retry refinement).
    //
    // The Bio package init probe above guarantees init COMPLETION (via
    // grok.functions.call('Bio:getSeqHelper') — atlas bio.cp.bio-service-
    // surface-init). It does NOT guarantee that the top-menu LEAF functions
    // dispatched by `[name="div-Bio---Analyze---Sequence-Space..."]` /
    // `Activity-Cliffs...` / `Composition` are findable in the function
    // registry yet — leaf registration is a distinct code path that can lag
    // init completion by a short window on truly-cold MSA boots (the
    // Gate-B FLAKY attempt-1 signature: 137.9s wall-clock with the editor
    // dialog never materialising within the 15s wait). Bounded poll until
    // each target Analyze leaf is findable as a Datagrok function, with a
    // short defensive settle if any are missing on this Bio build.
    //
    // Same-paradigm tactical reinforcement: only tightens an existing wait
    // step on the existing trigger mechanism. NOT a paradigm pivot.
    await page.evaluate(async () => {
      const candidates = [
        ['Bio:sequenceSpaceTopMenu', 'Bio:sequenceSpace'],
        ['Bio:activityCliffsTopMenu', 'Bio:activityCliffs', 'Bio:macromoleculeActivityCliffs'],
        ['Bio:compositionAnalysisWidget', 'Bio:composition', 'Bio:compositionAnalysis'],
      ];
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
        if (candidates.every(findAny)) return;
        await new Promise((r) => setTimeout(r, 300));
      }
      // Even if some candidate names are not findable (function rename across
      // Bio versions), do not error — the per-step 60s dialog tolerance below
      // is the defensive ceiling. Short settle so the menu dispatch doesn't
      // race a still-completing function-table install.
      await new Promise((r) => setTimeout(r, 1500));
    });
    await page.waitForTimeout(2000);

    // Scenario 1 — Open dataset (done above) + dispatch each Analyze function.
    // The "submenu exposes Sequence Space.../Activity Cliffs.../Composition"
    // assertion is verified implicitly by clicking the leaf in each
    // subsequent softStep — a missing leaf would surface as a click failure.

    // Scenario 2 — Run each function with default parameters.

    await softStep(`${ds.name}: Bio > Analyze > Sequence Space — run with defaults`, async () => {
      // Top-menu navigation per bio.md ("Click pattern (MCP-validated, mirrors chem.md)"):
      // click [name="div-Bio"] → 400ms → mouseover [name="div-Bio---Analyze"] → 300ms →
      // click leaf. Hover-not-click on the group is required to surface the leaves.
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Analyze"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Sequence-Space..."]') as HTMLElement).click();
      });
      // Dialog opens (atlas bio.analyze.sequence-space.editor). Cold-start tolerance:
      // the Bio package function dispatch can take >15s on a first-ever cold MSA
      // boot (Gate-B FLAKY attempt-1 evidence: dialog never materialized within
      // 15s on cold MSA; attempts 2-3 warm: <12s). 60s tolerates the cold ceiling
      // observed empirically while staying well below the 240s viewer-compute
      // budget used for the post-OK assertion.
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
      const title = await page.locator('.d4-dialog .d4-dialog-title').textContent();
      expect(title?.trim()).toBe('Sequence Space');

      // Click OK — run with default parameters.
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();

      // Verify ScatterPlot opens with embedding columns appended
      // (atlas bio.analyze.sequence-space.transform).
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 240_000});
      const hasScatter = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'));
      expect(hasScatter).toBe(true);
    });

    // Close non-Grid viewers so the next softStep observes its own Scatter plot.
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });

    await softStep(`${ds.name}: Bio > Analyze > Activity Cliffs — run with defaults`, async () => {
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Analyze"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Activity-Cliffs..."]') as HTMLElement).click();
      });
      // Cold-start tolerance: same rationale as Sequence Space step above —
      // 60s tolerates the empirical cold ceiling. Activity Cliffs follows
      // Sequence Space in the same test() body, so by the time we reach this
      // step the Bio package is already warm; the 60s tolerance is defensive.
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
      const title = await page.locator('.d4-dialog .d4-dialog-title').textContent();
      expect(title?.trim()).toBe('Activity Cliffs');

      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();

      // Verify ScatterPlot with cliff overlay (atlas bio.analyze.activity-cliffs.init).
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 240_000});
      const hasScatter = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'));
      expect(hasScatter).toBe(true);
    });

    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });

    await softStep(`${ds.name}: Bio > Analyze > Composition — WebLogo docks (no dialog)`, async () => {
      // Composition has NO "..." suffix — opens directly with NO dialog
      // per bio.md ("opens directly with NO dialog, mirroring the Chem
      // Scaffold Tree pattern"). Docks a WebLogo viewer.
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Analyze"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Composition"]') as HTMLElement).click();
      });

      // Verify WebLogo viewer docks (atlas bio.analyze.composition →
      // bio.viewers.web-logo). Pixel-level WebLogo paint is atlas
      // manual_only — only presence asserted here.
      await page.waitForFunction(
        () => Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'WebLogo'),
        null, {timeout: 60_000});
      const hasWebLogo = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'WebLogo'));
      expect(hasWebLogo).toBe(true);
    });

    // Scenario 3 — Composition Gear → Context Panel wiring.
    // Per scenario: only run on filter_FASTA.csv after Composition above.
    if (ds.name === 'FASTA') {
      await softStep(`${ds.name}: Composition Gear → Context Panel binds WebLogo property surface`, async () => {
        // Per bio.md "WebLogo Gear icon" recon-note (line 622): WebLogo title
        // bar does NOT surface a Gear icon under body.selenium on Bio 2.26.5.
        // The deterministic equivalent is `grok.shell.o = wlInstance` — opens
        // the same Context Panel property surface that the (absent) Gear
        // would open. Scenario Step 5 asserts only "Context Panel opens with
        // the Composition viewer's property surface bound (presence of
        // property editor surface only)" — see scenario SR-02 deferral of
        // the concrete property checklist.
        const result: {hasPropPanel: boolean, hasPropertyGrid: boolean} = await page.evaluate(async () => {
          const wl = (grok.shell.tv as any).viewers.find((v: any) => v.type === 'WebLogo');
          (grok as any).shell.o = wl;
          // Allow Context Panel to materialize.
          await new Promise((r) => setTimeout(r, 800));
          const propPanel = document.querySelector('.grok-prop-panel');
          const propertyGrid = document.querySelector('.grok-prop-panel .property-grid');
          return {hasPropPanel: !!propPanel, hasPropertyGrid: !!propertyGrid};
        });
        expect(result.hasPropPanel).toBe(true);
        expect(result.hasPropertyGrid).toBe(true);
      });
    }

    // Final cleanup: close non-Grid viewers.
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
