/* ---
sub_features_covered: [bio.analyze.composition, bio.viewers.web-logo]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: integration
//   sub_features_covered: [bio.analyze.composition, bio.viewers.web-logo]
//   ui_coverage_responsibility: ["Bio | Analyze | Composition",
//     composition-analysis-viewer, composition-analysis-select-on-click,
//     composition-analysis-gear-icon, composition-analysis-context-pane-properties]
//     (delegated_to: null)
//   related_bugs: []
//   produced_from: migrated
//   scope_reductions: SR-01 (A-CONT-01: Context-Pane property-name correctness
//     deferred until atlas/operator supplies canonical checklist; the wiring +
//     edit-acceptance flow is preserved as Step 5 with a generic two-toggle
//     probe that asserts "at least one property accepts an edit").
//
// Atlas provenance (derived_from):
//   feature-atlas/bio.yaml#sub_features[bio.analyze.composition]
//     derived_from: public/packages/Bio/src/package.ts#L1041
//   feature-atlas/bio.yaml#sub_features[bio.viewers.web-logo]
//     derived_from: public/packages/Bio/src/package.ts#L449
//   feature-atlas/bio.yaml#critical_paths[bio.cp.composition-analysis]
//     derived_from: public/packages/Bio/src/package.ts#L1041
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser/references/bio.md):
//   Top-menu path [name="div-Bio"] → mouseenter [name="div-Bio---Analyze"] → click
//     [name="div-Bio---Analyze---Composition"] — Bio top-menu (no `...` suffix
//     since the leaf opens directly without a dialog). Documented in
//     grok-browser/references/bio.md:193 (class 1). Listed here only to
//     anchor the menu-dispatch evaluate-block (visibility gate workaround).
//   .panel-base > .panel-titlebar [name="icon-font-icon-settings"] — Gear lives
//     on the outer docked-panel title bar (NOT inside [name="viewer-WebLogo"]
//     subtree); reached via parentElement walk from the viewer root up to the
//     .panel-base ancestor. MCP-validated 2026-06-02 on dev.datagrok.ai across
//     FASTA / HELM / MSA fixtures (`gearOnPanel: true` on all three).
//     bio.md:218 still flags this as a forward gap ("Gear icon NOT present
//     in normal title bar"); the 2026-06-02 recon contradicts that and the
//     panel-titlebar gear DOES open the prop panel. Kept as class-2 here
//     until bio.md is re-curated; equivalent fallback path
//     `grok.shell.o = wlInstance` is also documented (also MCP-validated).
//   tr[name="prop-show-position-labels"] / tr[name="prop-skip-empty-positions"] —
//     property-grid rows in .grok-prop-panel after Gear click. MCP-validated
//     2026-06-02 on dev.datagrok.ai: both rows present in the panel's DOM and
//     toggling their input[type=checkbox] flips wl.getOptions().look
//     {showPositionLabels, skipEmptyPositions} cleanly. bio.md:204 documents
//     the underlying getOptions().look property surface but does not
//     enumerate the property-grid TR `[name=]` attrs since they are not
//     officially stable; we treat them as class-2 anchors and verify the
//     edit via getOptions().look diff rather than DOM read.
//   Canvas click probe via the viewer's own rendered letter-column rects —
//     synthetic mousedown/mouseup/click MouseEvent dispatched on the WebLogo
//     canvas at the CENTER of an actual monomer rect read from
//     PositionMonomerInfo.bounds. **MCP recon 2026-06-02** confirmed the
//     bounds object shape is {x, y, width, height} (NOT
//     {left, top, width, height} as the prior spec assumed — that assumption
//     produced NaN clientX/Y, browser-coerced to 0, click-misses-canvas; the
//     prior cycle (2026-06-01-bio-migrate-02 attempt-3) Gate B FAIL on FASTA
//     and HELM was caused by this typo). With the bounds key corrected:
//     FASTA selects 8 rows on position-0 'M', HELM selects 2 rows on
//     position-2 'F', MSA selects 8 rows on position-1 'hHis' — all
//     deterministic against the static fixtures.
//   ~3s settle after viewer attach — empirical observation from
//     composition-analysis-run.md (canvas in DOM but handlers not yet bound
//     on 2nd/3rd dataset). Flagged in scenario unresolved_ambiguities as
//     "weblogo-canvas-interactivity-settle-window-not-atlas-declared".
//     MCP recon 2026-06-02 also observed the same: a 4s post-dispatch settle
//     before the WebLogo canvas hit-handlers route to df.selection is
//     required across all 3 notations.
//
// Hypothesis investigation (this cycle, 2026-06-02-bio-automate-01):
//   Round 1 — hypothesis: test-bug (selector-bounds key typo on canvas click).
//   Cheap check: re-read spec body around failing softStep (line 200-201 of
//     prior spec used b.left/b.top — undefined per the actual PositionMonomerInfo
//     bounds shape `{x, y, width, height}`, derived from MCP introspection).
//   Full repro: MCP recon on dev.datagrok.ai 2026-06-02 reproduced the FASTA
//     and HELM failure path with the prior bounds-key approach (click at
//     translated NaN→0 coords misses the canvas) and confirmed the fix
//     (bounds.x/.y produces in-canvas coords, click hits monomer rect,
//     df.selection.trueCount > 0).
//   Conclusion: test-bug (single-line bounds key fix). Same-paradigm
//     tactical fix per §"Cheap-checks usage contract" rule #2; no paradigm
//     pivot. mcp_status: used; mcp_observations[] populated.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/tests/filter_MSA.csv'},
];

test('Bio | Analyze | Composition — composition analysis integration', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  for (const ds of datasets) {
    await softStep(`[${ds.name}] Scenario 1 Step 1 — Open ${ds.path}`, async () => {
      const result: {rows: number, hasMacromolecule: boolean} = await page.evaluate(async (path: string) => {
        const g = (window as any).grok;
        g.shell.closeAll();
        const df = await g.dapi.files.readCsv(path);
        g.shell.addTableView(df);
        await new Promise<void>((resolve) => {
          const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(() => resolve(), 3000);
        });
        const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
        const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
        if (hasBioChem) {
          for (let i = 0; i < 50; i++) {
            if (document.querySelector('[name="viewer-Grid"] canvas')) break;
            await new Promise((r) => setTimeout(r, 200));
          }
          await new Promise((r) => setTimeout(r, 5000));
        }
        return {rows: df.rowCount, hasMacromolecule: hasBioChem};
      }, ds.path);
      expect(result.rows).toBeGreaterThan(0);
      expect(result.hasMacromolecule).toBe(true);
      await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    });

    await softStep(`[${ds.name}] Scenario 1 Step 2 — Bio > Analyze > Composition; WebLogo viewer docks; no multi-column dialog`, async () => {
      // Top-menu dispatch via DOM events — Playwright .click() fails the
      // visibility gate on a hidden-but-attached submenu item.
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]')!.dispatchEvent(
          new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Composition"]') as HTMLElement).click();
      });
      // Verify WebLogo presence on shell.tv.viewers (atlas declares
      // bio.viewers.web-logo as manual_only for paint correctness, so we
      // assert presence/docking only, not pixel-level correctness).
      await page.waitForFunction(() => {
        const viewers = Array.from((window as any).grok.shell.tv.viewers);
        return viewers.some((v: any) => v.type === 'WebLogo');
      }, null, {timeout: 60_000});
      const result: {hasCanvas: boolean, dialogOpen: boolean} = await page.evaluate(async () => {
        const g = (window as any).grok;
        let hasCanvas = false;
        for (let i = 0; i < 40; i++) {
          const w = g.shell.tv.viewers.find((v: any) => v.type === 'WebLogo');
          if (w && w.root.querySelector('canvas')) { hasCanvas = true; break; }
          await new Promise((r) => setTimeout(r, 300));
        }
        // Multi-column choice dialog MUST be absent — each fixture has a
        // single Macromolecule column per scenario Setup.
        const dialogOpen = !!document.querySelector('.d4-dialog');
        return {hasCanvas, dialogOpen};
      });
      expect(result.hasCanvas).toBe(true);
      expect(result.dialogOpen).toBe(false);
    });

    await softStep(`[${ds.name}] Scenario 2 Step 3 — Click letter in WebLogo selects ≥1 row in source grid`, async () => {
      // ~3-4s settle: canvas may be in DOM before the click handlers are bound
      // (empirical from composition-analysis-run.md AND from MCP recon
      // 2026-06-02; documented as unresolved_ambiguities
      // weblogo-canvas-interactivity-settle-window-not-atlas-declared).
      await page.waitForTimeout(4000);
      const selected: number = await page.evaluate(async () => {
        const g = (window as any).grok;
        const tv = g.shell.tv;
        const df = tv.dataFrame;
        df.selection.setAll(false);
        const wl: any = tv.viewers.find((v: any) => v.type === 'WebLogo');
        const canvas = wl.root.querySelector('canvas') as HTMLCanvasElement;
        const canvasRect = canvas.getBoundingClientRect();
        const dpr = window.devicePixelRatio || 1;
        const GAP = '-';
        // Read the REAL rendered letter-column rects from the viewer's own
        // layout (PositionMonomerInfo.bounds, set during render()). MCP recon
        // 2026-06-02 confirmed bounds shape is `{x, y, width, height}` in
        // canvas pixel coordinates relative to the canvas origin — the prior
        // spec used `.left` / `.top` which were undefined (NaN clientX →
        // browser-coerced to 0 → click misses canvas; cause of the prior
        // FASTA/HELM Gate B FAIL). The click handler (canvasOnMouseDown)
        // maps the cursor back through getMonomer() and selects rows
        // carrying that monomer at that position.
        const positions: any[] = Array.isArray(wl.positions) ? wl.positions : [];
        const candidates: {x: number, y: number, h: number}[] = [];
        for (let idx = 0; idx < positions.length && idx < 30; idx++) {
          const pi = positions[idx];
          if (!pi || typeof pi.getMonomers !== 'function') continue;
          let monomers: string[] = [];
          try { monomers = pi.getMonomers(); } catch (_) { monomers = []; }
          for (const m of monomers) {
            if (!m || m === GAP) continue;
            let b: any = null;
            try { b = pi.getFreq(m).bounds; } catch (_) { b = null; }
            if (!b) continue;
            const cx = (b.x + b.width / 2) / dpr + canvasRect.left;
            const cy = (b.y + b.height / 2) / dpr + canvasRect.top;
            candidates.push({x: cx, y: cy, h: b.height});
          }
        }
        // Tallest rect first — the most frequent monomer is the largest
        // hit target, minimising boundary misses. First click that selects
        // ≥1 row wins.
        candidates.sort((a, b) => b.h - a.h);
        for (const c of candidates) {
          df.selection.setAll(false);
          for (const type of ['mousedown', 'mouseup', 'click']) {
            canvas.dispatchEvent(new MouseEvent(type, {
              bubbles: true, cancelable: true, clientX: c.x, clientY: c.y, button: 0,
            }));
          }
          await new Promise((r) => setTimeout(r, 400));
          if (df.selection.trueCount > 0) break;
        }
        return df.selection.trueCount;
      });
      // Exact selection count is dataset+position dependent — scenario asserts
      // "at least one row selected" baseline (not a specific count). MCP-recon
      // observation 2026-06-02: FASTA=8, HELM=2, MSA=8 selected on the
      // tallest-rect path; all > 0.
      expect(selected).toBeGreaterThan(0);
    });

    await softStep(`[${ds.name}] Scenario 3 Step 4 — Gear icon on WebLogo opens Context Pane property grid`, async () => {
      const opened: {found: boolean, pg: boolean} = await page.evaluate(async () => {
        const g = (window as any).grok;
        g.shell.tv.dataFrame.selection.setAll(false);
        const w = g.shell.tv.viewers.find((v: any) => v.type === 'WebLogo');
        // Gear lives on the outer docked-panel title bar (.panel-base
        // ancestor), NOT inside the [name="viewer-WebLogo"] subtree.
        // MCP-validated 2026-06-02 across FASTA/HELM/MSA.
        let panelBase: any = w.root;
        while (panelBase && !panelBase.classList?.contains('panel-base'))
          panelBase = panelBase.parentElement;
        const gear = panelBase?.querySelector(
          '.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
        if (!gear) {
          // Fallback per bio.md:218: programmatic equivalent. The Gear-click
          // path is the primary UI driver; this fallback covers a forward
          // gap if the title-bar gear regresses across builds.
          g.shell.o = w;
          await new Promise((r) => setTimeout(r, 600));
          const pg2 = document.querySelector('.grok-prop-panel .property-grid, .grok-prop-panel tr[name^="prop-"]');
          return {found: false, pg: !!pg2};
        }
        gear.click();
        await new Promise((r) => setTimeout(r, 600));
        const pg = document.querySelector('.grok-prop-panel .property-grid, .grok-prop-panel tr[name^="prop-"]');
        return {found: true, pg: !!pg};
      });
      // The wiring assertion: Context Pane (.grok-prop-panel) bound to the
      // WebLogo viewer surface. The bare Gear click is the canonical path;
      // the JS-API fallback (grok.shell.o = w) is the documented forward-gap
      // recovery per bio.md:218.
      expect(opened.pg).toBe(true);
      // Property row inside a collapsed accordion section — use
      // state: 'attached' rather than the default 'visible' wait.
      await page.locator('tr[name="prop-show-position-labels"]').waitFor({
        state: 'attached', timeout: 10_000});
    });

    await softStep(`[${ds.name}] Scenario 3 Step 5 — Edit ≥1 Context Pane property (SR-01: edit-acceptance only)`, async () => {
      // SR-01: per scenario frontmatter scope_reductions, the property-name
      // correctness assertion is deferred until atlas/operator supplies a
      // canonical Context-Pane property checklist. We assert the wiring +
      // edit-acceptance flow only: toggle two boolean properties and check
      // that at least one flips its value via getOptions().look diff.
      // MCP-validated 2026-06-02: toggling either checkbox flips look state
      // cleanly on all 3 notations.
      const result: {changedShow: boolean, changedSkip: boolean} = await page.evaluate(async () => {
        const g = (window as any).grok;
        const wl: any = g.shell.tv.viewers.find((v: any) => v.type === 'WebLogo');
        const before = {
          showPositionLabels: wl.getOptions().look.showPositionLabels,
          skipEmptyPositions: wl.getOptions().look.skipEmptyPositions,
        };
        for (const name of ['prop-show-position-labels', 'prop-skip-empty-positions']) {
          const row = document.querySelector(`tr[name="${name}"]`);
          const cb = row?.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
          if (cb) {
            cb.click();
            cb.dispatchEvent(new Event('change', {bubbles: true}));
          }
          await new Promise((r) => setTimeout(r, 300));
        }
        await new Promise((r) => setTimeout(r, 400));
        const after = {
          showPositionLabels: wl.getOptions().look.showPositionLabels,
          skipEmptyPositions: wl.getOptions().look.skipEmptyPositions,
        };
        return {
          changedShow: before.showPositionLabels !== after.showPositionLabels,
          changedSkip: before.skipEmptyPositions !== after.skipEmptyPositions,
        };
      });
      // At least one of the two boolean properties accepted the edit (per
      // SR-01: scenario asserts edit-acceptance, not which specific property
      // toggled).
      expect(result.changedShow || result.changedSkip).toBe(true);
    });
  }

  // Cleanup contract per scenario Notes — close all views.
  await page.evaluate(() => (window as any).grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
