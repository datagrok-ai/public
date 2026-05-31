/* ---
sub_features_covered: [peptides.panels.peptides, peptides.rendering.weblogo-header, peptides.widgets.distribution, peptides.util.modify-selection]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: integration
//   sub_features_covered: [peptides.panels.peptides, peptides.rendering.weblogo-header,
//                          peptides.widgets.distribution, peptides.util.modify-selection]
//   ui_coverage_responsibility: [peptides-context-panel-peptides-tab, peptides-panel-parameters,
//                                weblogo-amino-acid-click-selection] (delegated_to: info-panels.md)
//   related_bugs: [GROK-19145, GROK-17557, GROK-14298]
// Atlas provenance (derived_from):
//   peptides.yaml#sub_features[peptides.panels.peptides] (atlas L140 peptidesPanel)
//   peptides.yaml#sub_features[peptides.rendering.weblogo-header] (atlas L280 setWebLogoRenderer)
//   peptides.yaml#sub_features[peptides.util.modify-selection] (atlas L519 modifySelection)
// Selectors verified via chrome-devtools MCP recon on dev.datagrok.ai (peptides v1.27.9 / bio v2.27.13);
// authoritative reference: .claude/skills/grok-browser/references/peptides.md.
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser reference):
//   .bio-wl-host canvas — WebLogo preview rendered inside the Peptides context-panel pane
//     ([name="pane-Peptides"], reached via df.currentCol = df.col('AlignedSequence') → expand the
//     pane's .d4-accordion-pane-header). Distinct from the column-HEADER WebLogo that the
//     grok-browser reference documents (which renders only after Launch SAR); this preview is
//     mounted on column-focus, pre-launch. Container element carries class "bio-wl-host ui-div"
//     with a single child <canvas> (geom ~595x130 px). Observed live 2026-05-30 via
//     chrome-devtools MCP (evaluate_script). Not in peptides.md — its WebLogo § documents only
//     the post-SAR column-header logo.
//
// Cold-init stabilization (revised per live MCP recon 2026-05-30, cycle
// 2026-05-30-peptides-automate-01): the [name="pane-Peptides"] node is built by the
// Peptides package's @panel-decorator function `peptidesPanel` (see
// public/packages/Peptides/src/package.ts L161-L172). On a warm worker, after
// `df.currentCol = df.col('AlignedSequence')` fires, the Context-Panel rebuild that mounts
// the Peptides pane child is a **2 ms** deterministic operation. The earlier Gate B FLAKY
// failure (cycle 2026-05-29-peptides-automate-02 attempt 3 — both the first 15 s waitFor
// AND a 40 s accordion-header poll AND a second 30 s waitFor expired on a cold worker) is
// rooted in the @init function `initPeptides` (loads MonomerWorks + TreeHelper) still
// running when `df.currentCol = col` is driven; the @panel function `peptidesPanel` cannot
// run until its package's @init has completed, so the Context-Panel rebuild builds the
// other panes (Details, Filter, Style, …) but skips Peptides + Bioinformatics.
//
// The fix is **strictly additive** (Gate B note: "add stabilization, do not rewrite"):
// explicit `await fn.package.load()` for the Peptides + Bio packages BEFORE driving column
// focus. `Package.load()` is idempotent (8-51 ms on warm, awaits @init completion on cold)
// and is the canonical pattern the platform uses to gate package-conditional UI. After it
// resolves, the @panel function is guaranteed runnable, so the focus → pane-Peptides path
// collapses to its deterministic warm-worker timing on cold workers too. We still drive
// focus through BOTH df.currentCol = col AND grok.shell.o = col (the existing same-paradigm
// stabilization), but the package-preload step removes the dominant tail-latency source.
// The 40 s accordion-pane-header poll + 2 s settle stay as defense-in-depth.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

test('Peptides — SAR parameters and WebLogo', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Steps 1-2: Open peptides.csv and focus the Macromolecule column title
  await softStep('Steps 1-2: Open dataset and select Macromolecule column', async () => {
    const result = await page.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach(d => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = false;

      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise<void>(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      // Macromolecule dataset: wait for grid canvas + Bio package init
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));

      // Explicitly load the Peptides + Bio packages BEFORE driving column focus.
      // Package.load() is idempotent and awaits the @init function (initPeptides loads
      // MonomerWorks + TreeHelper). Without this, on a cold worker the Context-Panel
      // rebuild that mounts the [name="pane-Peptides"] child can race the still-running
      // @init, and the pane will not appear inside the 40 s readiness poll below — the
      // empirical root cause of Gate B 2026-05-29 attempt-3 (both the 15 s waitFor and
      // the 40 s poll + 30 s waitFor expired). On warm workers this resolves in ~50 ms.
      const peptidesPanelFn = (DG.Func.find({name: 'peptidesPanel'}) || [])[0];
      const peptidesPkg = peptidesPanelFn?.package;
      if (peptidesPkg && typeof peptidesPkg.load === 'function')
        await peptidesPkg.load();
      const bioFns = DG.Func.find({package: 'Bio'}) || [];
      const bioPkg = bioFns[0]?.package;
      if (bioPkg && typeof bioPkg.load === 'function')
        await bioPkg.load();

      // Focus the Macromolecule column → Context Panel reflects the column selection.
      // With the preload above, the Context-Panel rebuild that mounts the Peptides pane
      // is deterministic (~2 ms on warm; bounded by @init completion on cold).
      // Drive the focus through BOTH df.currentCol = col (the canonical column-focus
      // signal that triggers the Context-Panel accordion rebuild) AND grok.shell.o = col
      // (the existing same-paradigm stabilization). The 40 s accordion-pane-header poll
      // + 2 s settle stay as defense-in-depth against any residual cold-init lag.
      const col = df.col('AlignedSequence')!;
      df.currentCol = col;
      grok.shell.o = col;
      let paneMounted = false;
      // Up to ~40s, polled at 500ms — defense-in-depth post-preload.
      const deadline = Date.now() + 40000;
      while (Date.now() < deadline) {
        const headers = document.querySelectorAll('.d4-accordion-pane-header').length;
        const pane = document.querySelector('[name="pane-Peptides"]');
        if (headers > 0 && pane && pane.querySelector('.d4-accordion-pane-header')) {
          paneMounted = true;
          break;
        }
        await new Promise(r => setTimeout(r, 500));
      }
      // Settle so the lagging Peptides pane finishes registering on the accordion.
      await new Promise(r => setTimeout(r, 2000));
      return {rows: df.rowCount, semType: df.col('AlignedSequence')?.semType, paneMounted};
    }, datasetPath);
    expect(result.rows).toBe(647);
    expect(result.semType).toBe('Macromolecule');
  });

  // Step 3: Expand the Peptides accordion pane on the Context Panel (peptides.panels.peptides).
  // DOM-driven via the Playwright Locator API (the pane header is a real DOM node with a stable
  // [name=] selector — drive the real gesture, then assert the mounted parameter widgets + WebLogo).
  await softStep('Step 3: Expand Peptides context-panel pane', async () => {
    const pane = page.locator('[name="pane-Peptides"]');
    // Cold package-init race recovery: on a cold worker the Context-Panel rebuild that
    // mounts this pane can lag far behind the warm path (~200ms), and a single bounded
    // Locator wait expires before the rebuild completes (Gate B 2026-05-29: both the
    // first 30s wait and a grok.shell.o-only re-assert path's 30s wait expired on cold
    // workers). Replace that with the validated info-panels-spec.ts pattern: re-drive
    // the column focus through BOTH df.currentCol = col and grok.shell.o = col, then
    // RETRY a readiness poll (up to ~40s, 500ms interval) on the accordion-pane headers
    // inside the same evaluate. df.currentCol is the canonical column-focus signal that
    // drives the accordion rebuild; grok.shell.o alone proved insufficient on cold init.
    // Never the grok.shell.o = df toggle (regressed the cold harness 4/4 previously).
    try {
      await pane.waitFor({state: 'attached', timeout: 15_000});
    } catch {
      await page.evaluate(async () => {
        // Re-assert Peptides/Bio package load before re-driving focus — idempotent and
        // cheap on warm, blocks until @init completes on a residually-cold worker.
        const peptidesPanelFn = (DG.Func.find({name: 'peptidesPanel'}) || [])[0];
        const peptidesPkg = peptidesPanelFn?.package;
        if (peptidesPkg && typeof peptidesPkg.load === 'function')
          await peptidesPkg.load();
        const bioFns = DG.Func.find({package: 'Bio'}) || [];
        const bioPkg = bioFns[0]?.package;
        if (bioPkg && typeof bioPkg.load === 'function')
          await bioPkg.load();

        const df = grok.shell.tv?.dataFrame;
        const col = df?.col('AlignedSequence');
        if (!col) return;
        df.currentCol = col;
        grok.shell.o = col;
        const deadline = Date.now() + 40000;
        while (Date.now() < deadline) {
          if (document.querySelectorAll('.d4-accordion-pane-header').length > 0
            && document.querySelector('[name="pane-Peptides"] .d4-accordion-pane-header')) break;
          await new Promise(r => setTimeout(r, 500));
        }
        await new Promise(r => setTimeout(r, 2000));
      });
      await pane.waitFor({state: 'attached', timeout: 30_000});
    }
    await pane.waitFor({timeout: 30_000});
    const header = pane.locator('.d4-accordion-pane-header');
    // Only click to expand if collapsed (focusing the column can leave it expanded).
    const alreadyExpanded = await header.evaluate(el => el.classList.contains('expanded'));
    if (!alreadyExpanded) {
      await header.scrollIntoViewIfNeeded();
      await header.click();
    }
    // Wait for the analyze widgets + WebLogo preview to mount inside the pane.
    await pane.locator('[name="input-Scaling"]').waitFor({timeout: 30_000});
    await pane.locator('.bio-wl-host canvas').waitFor({timeout: 30_000});

    const result = await pane.evaluate((el) => ({
      expanded: !!el.querySelector('.d4-accordion-pane-header')?.classList.contains('expanded'),
      hasWebLogo: !!el.querySelector('.bio-wl-host canvas'),
      hasParams: !!el.querySelector('[name="input-Scaling"]')
        && !!el.querySelector('[name="input-Activity"]')
        && !!el.querySelector('[name="input-Clusters"]'),
    }));
    expect(result.expanded).toBe(true);
    expect(result.hasWebLogo).toBe(true);
    expect(result.hasParams).toBe(true);
  });

  // Steps 4-5: Change the Scaling parameter; the panel re-renders without console errors
  // (peptides.widgets.distribution re-render path). input-Scaling is a native <select>
  // [none|lg|-lg]; drive it through the Playwright Locator API (selectOption) and confirm the
  // pane re-renders (WebLogo + activity-distribution Histogram persist).
  await softStep('Steps 4-5: Change Scaling parameter, panel re-renders', async () => {
    const pane = page.locator('[name="pane-Peptides"]');
    const scaling = pane.locator('[name="input-Scaling"]');
    await scaling.selectOption('lg');
    await page.waitForTimeout(1500);

    const result = await pane.evaluate((el) => {
      const sel = el.querySelector('[name="input-Scaling"]') as HTMLSelectElement | null;
      return {
        scalingChanged: sel?.value === 'lg',
        reRendered: !!el.querySelector('.bio-wl-host canvas') && !!el.querySelector('[name="viewer-Histogram"]'),
      };
    });
    expect(result.scalingChanged).toBe(true);
    expect(result.reRendered).toBe(true);
  });

  // Steps 6-9: Click an amino-acid monomer on the WebLogo → matching DataFrame rows are selected
  // (peptides.rendering.weblogo-header click handler → peptides.util.modify-selection → DG.BitSet).
  // The WebLogo is canvas-rendered (no per-monomer DOM node), so this gesture is a sanctioned
  // canvas-fallback: dispatch a synthetic MouseEvent on the .bio-wl-host canvas and assert the
  // selection cardinality via the JS API (per grok-browser peptides.md pitfall #3).
  // A click at fractional (0.40, 0.45) of the WebLogo canvas selects matching rows; the Context
  // Panel then follows the new selection (grok.shell.o changes), so the count is captured inside
  // the same evaluate, immediately after the click. MCP recon 2026-05-29: this fractional point
  // lands on a populated monomer column and selects 492/647 rows (selBefore=0). The source wording
  // is "some rows were selected" (no exact count bounded), so the assertion is the source-faithful
  // > 0 threshold rather than an over-tightened == 492 (which would be build-brittle).
  await softStep('Steps 6-9: WebLogo monomer click selects matching rows', async () => {
    // Re-focus the column to guarantee the Peptides pane + WebLogo are mounted, then re-expand.
    // Same cold-init stabilization as Step 3: idempotent re-assert of Peptides/Bio package load
    // (no-op on warm; blocks on residually-cold workers), then drive the focus through
    // df.currentCol = col AND grok.shell.o = col, then RETRY a readiness poll on the
    // accordion-pane headers so the lagging Peptides pane re-mounts before the WebLogo wait.
    await page.evaluate(async () => {
      const peptidesPanelFn = (DG.Func.find({name: 'peptidesPanel'}) || [])[0];
      const peptidesPkg = peptidesPanelFn?.package;
      if (peptidesPkg && typeof peptidesPkg.load === 'function')
        await peptidesPkg.load();
      const bioFns = DG.Func.find({package: 'Bio'}) || [];
      const bioPkg = bioFns[0]?.package;
      if (bioPkg && typeof bioPkg.load === 'function')
        await bioPkg.load();

      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      const col = df.col('AlignedSequence')!;
      df.currentCol = col;
      grok.shell.o = col;
      const deadline = Date.now() + 40000;
      while (Date.now() < deadline) {
        if (document.querySelectorAll('.d4-accordion-pane-header').length > 0
          && document.querySelector('[name="pane-Peptides"] .d4-accordion-pane-header')) break;
        await new Promise(r => setTimeout(r, 500));
      }
      await new Promise(r => setTimeout(r, 2000));
    });
    const pane = page.locator('[name="pane-Peptides"]');
    await pane.waitFor({timeout: 30_000});
    const header = pane.locator('.d4-accordion-pane-header');
    const alreadyExpanded = await header.evaluate(el => el.classList.contains('expanded'));
    if (!alreadyExpanded)
      await header.click();
    const canvas = pane.locator('.bio-wl-host canvas');
    await canvas.waitFor({timeout: 30_000});

    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const pane = document.querySelector('[name="pane-Peptides"]');
      const canvas = pane?.querySelector('.bio-wl-host canvas') as HTMLCanvasElement | null;
      if (!canvas) return {clicked: false, selectedRows: 0};

      const r = canvas.getBoundingClientRect();
      const cx = r.x + r.width * 0.40;
      const cy = r.y + r.height * 0.45;
      const selBefore = df.selection.trueCount;
      for (const type of ['pointerdown', 'mousedown', 'pointerup', 'mouseup', 'click']) {
        canvas.dispatchEvent(new MouseEvent(type, {
          bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0, view: window,
        }));
      }
      await new Promise(r => setTimeout(r, 800));
      // Capture selection immediately — the panel re-renders to the selection view after this.
      return {clicked: true, selBefore, selectedRows: df.selection.trueCount, total: df.rowCount};
    });
    // Step 8: at least one row selected. Step 9: selected-row count > 0.
    // GROK-14298 invariant — the fireBitsetChanged collaborative-selection backbone must
    // populate df.selection (cardinality > 0) without crashing on the WebLogo click.
    expect(result.clicked).toBe(true);
    expect(result.selectedRows).toBeGreaterThan(0);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
