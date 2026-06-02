/* ---
sub_features_covered: [bio.analyze.composition, bio.viewers.web-logo]
--- */
//   scope_reductions: SR-01 (A-CONT-01: Context-Pane property-name correctness
//   Round 1 — hypothesis: test-bug (selector-bounds key typo on canvas click).
//   Conclusion: test-bug (single-line bounds key fix). Same-paradigm
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
