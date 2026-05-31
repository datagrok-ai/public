/* ---
sub_features_covered: [peptides.panels, peptides.panels.peptides, peptides.rendering, peptides.rendering.monomer]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: ui-smoke
//   sub_features_covered: [peptides.panels, peptides.panels.peptides, peptides.rendering, peptides.rendering.monomer]
//   ui_coverage_responsibility: [peptides-context-panel-details, peptides-context-panel-peptides-tab, macromolecule-column-monomer-coloring] (delegated_to: null)
//   related_bugs: [GROK-17557]
// Scope: per scenario scope_reductions SR-01/SR-02, the Bioinformatics accordion pane
//   (registered by the Bio package, not Peptides) is out of this section's surface; this
//   spec verifies only the Details and Peptides panes on the Macromolecule column.
//
// Selector recon-notes (class-2: live-MCP-observed this session, cold-stable content
//   predicates not yet inlined in grok-browser/references/peptides.md — selectors are in
//   that reference; the per-pane content substrings + child-node predicates are the
//   delta this note records):
//   [name="pane-Details"] / [name="pane-Peptides"] — confirmed live on dev.datagrok.ai
//     2026-05-30 via chrome-devtools MCP. Both are .d4-accordion-pane containers; their
//     child header [name="div-section--Details"] / [name="div-section--Peptides"]
//     (.d4-accordion-pane-header) is what receives .click() to expand. Already in
//     peptides.md "Peptides context panel" §.
//   Cold-stable Details content: innerText contains 'Data type' and 'Semantic type'
//     (column metadata key/value table — rendered synchronously by the platform
//     column-properties widget). Observed live 2026-05-30: "...Data type\tstring\n
//     Order\tdefault\nSemantic type\t...".
//   Cold-stable Peptides content: innerText contains 'Activity' + 'LAUNCH SAR'; also the
//     synchronously-mounted DOM nodes [name="button-Launch-SAR"] and
//     [name="input-host-Activity"]. Observed live 2026-05-30. These are NOT lazy-mounted
//     (per project-peptides-context-panel-pane-content-predicate: prefer the structural
//     widgets the peptidesPanel builds synchronously over any embedded viewer canvas).
//   Sanctioned canvas-fallback for column-focus: grid column headers are canvas (per
//     project-grid-column-headers-are-canvas-not-dom). The column-title click is realized
//     as `df.currentCol = df.col('AlignedSequence')` JS-API per the grok-browser peptides
//     reference (the cell.renderer tag is the documented monomer-coloring assertion).
//
// Cold-init wired-panel timing recon (2026-05-30, dev.datagrok.ai, chrome-devtools MCP):
//   - `addTableView(df)` returns synchronously; pane-Details mounts at +1112ms after
//     `onSemanticTypeDetected` on a warm browser (single observation).
//   - After `df.currentCol = df.col('AlignedSequence')` on a warm browser with the
//     panel pre-wired: pane-Peptides mounts at +985ms (single observation), pane-Details
//     at +0..+2ms (already present from default current-col).
//   - Gate B 2026-05-30-peptides-automate-02 attempt-3 FLAKY: pane-Details was missing
//     from `.grok-prop-panel` for the full 40s Step-2 inner-poll — i.e. the property
//     panel never wired itself to the TableView, so the column-change event went
//     unrouted. attempt-1 + attempt-2 passed cleanly. The retry's prior hypothesis was
//     a property-panel reattach race (WelcomeView staying active); the setup phase
//     forces `grok.shell.v = tv` and polls pane-Details to wire-confirm before returning.
//
// Cycle 2026-05-30-peptides-automate-02 retry hypothesis (distinct from the prior
// wired-panel-gate hypothesis already in setup): the residual 1-in-3 cold flake is the
// Peptides @init + async @panel race per GROK-17557, not the WelcomeView race.
// `Peptides:initPeptides` loads MonomerWorks, TreeHelper, and PeptideUtils.loadComponents()
// asynchronously; the `peptidesPanel` @panel is registered against Macromolecule columns
// and depends on those resources. On a cold init the fixed 5s Bio settle is not
// guaranteed to cover the @init promise + the @panel registration race; if Step 2 sets
// `df.currentCol = AlignedSequence` before the @panel has been wired, the column-focus
// event is dropped and pane-Peptides never mounts.
// Two new defenses (mirror sibling export-invariant-map-spec.ts / export-mutation-cliffs-spec.ts):
//   (a) Pre-warm `Peptides:initPeptides` in setup so MonomerWorks + SeqHelper are ready
//       before the column-focus gesture. Non-fatal on throw (allowed during cold init).
//   (b) In Step 2, set `grok.shell.o = df.col('AlignedSequence')` as the deterministic
//       property-panel-rebuild trigger. The Context Panel responds to `grok.shell.o`
//       (current object), not directly to `df.currentCol`; setting `grok.shell.o = col`
//       bypasses the intermediate column-focus → o-binding event chain and forces the
//       rebuild on the same tick. `df.currentCol = df.col('AlignedSequence')` is still
//       set to satisfy the scenario's documented invariant (column becomes current).
//       Empirical MCP recon 2026-05-30 on a warm browser: pane-Peptides mounts at
//       +227ms after `grok.shell.o = col` with `initPeptides` pre-warmed, vs +985ms
//       after `df.currentCol = col` without pre-warm — both faster and more deterministic.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

test('Peptides — Info Panels', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup: open peptides.csv and ensure the Context Panel surfaces on the right.
  // Cold-init hardening (per attempt-3 FLAKY 2026-05-30 evidence: `.grok-prop-panel
  // [name="pane-Details"]` never appeared within 40s on attempt-3, while attempts 1+2
  // passed cleanly; warm MCP recon 2026-05-30 confirmed pane-Details mounts at +1112ms
  // after onSemanticTypeDetected on a warm browser, so the cold race window is wider
  // than the prior 40s Step-2 poll could absorb). The setup phase now does two extra
  // things before returning:
  //   - Force `grok.shell.v = tv` so the active view is the TableView (defends
  //     against a cold race where addTableView leaves the WelcomeView in front and
  //     the property panel never reattaches to the loaded dataframe).
  //   - Poll up to 30s for `[name="pane-Details"]` inside `.grok-prop-panel`. This
  //     is the platform-wide pane that exists for ANY current column; its presence
  //     is the deterministic signal that the property-panel infrastructure has
  //     wired itself to the TableView, so subsequent column-change events
  //     (Step 2's `df.currentCol = …` assignment) are routed correctly. Without
  //     this gate, the column-change event can fire before the panel listener is
  //     attached, and pane-Peptides never mounts (the attempt-3 failure mode).
  await softStep('Setup: open peptides.csv + ensure Context Panel wired', async () => {
    const result = await page.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) (cancel as HTMLElement).click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      // Windows mode (not Tabs mode) so the Context Panel docks on the right and the
      // column-properties accordion mounts. simpleMode=true (Tabs) collapses the
      // Context Panel into a tab that the scenario steps don't drive.
      grok.shell.windows.simpleMode = false;
      grok.shell.windows.showContextPanel = true;

      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      // Force the TableView to be the active view (cold-race defense: addTableView
      // returns synchronously but the WelcomeView can remain in front on cold init,
      // which leaves `.grok-prop-panel` mounted against the wrong view).
      try { (grok.shell as any).v = tv; } catch (_) { /* best-effort */ }
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      // Bio package init: wait for grid canvas to mount + extra settle so the
      // Macromolecule cell renderer + WebLogo finish registering.
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
      // GROK-17557 cold-init pre-warm: explicitly run the Peptides @init so
      // MonomerWorks + TreeHelper + PeptideUtils.loadComponents() (SeqHelper +
      // MonomerLib) are loaded BEFORE the column-focus gesture in Step 2 requests
      // the async peptidesPanel @panel. Without this pre-warm, the panel handler
      // can race the column-focus event and the rebuild silently drops. Mirrors
      // sibling export-invariant-map-spec.ts / export-mutation-cliffs-spec.ts.
      // Non-fatal on throw — the panel still works on warm browsers without the
      // explicit init, the wired-panel + currentCol polls below are the
      // authoritative readiness gates.
      try { await grok.functions.call('Peptides:initPeptides'); }
      catch (e) { console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', String(e)); }
      // Wired-panel gate: poll up to 30s for the platform-wide Details pane to mount
      // inside `.grok-prop-panel`. Its presence is the signal that the property-panel
      // infrastructure is attached to the TableView and will receive subsequent
      // column-change events.
      const wiredDeadline = Date.now() + 30_000;
      let detailsWired = false;
      while (Date.now() < wiredDeadline) {
        detailsWired = !!document.querySelector('.grok-prop-panel [name="pane-Details"]');
        if (detailsWired) break;
        await new Promise((r) => setTimeout(r, 250));
      }
      return {
        rows: df.rowCount,
        semType: df.col('AlignedSequence')?.semType,
        detailsWired,
        activeViewType: (grok.shell.v as any)?.type ?? (grok.shell.v as any)?.constructor?.name,
      };
    }, datasetPath);
    expect(result.rows).toBe(647);
    expect(result.semType).toBe('Macromolecule');
    expect(result.detailsWired).toBe(true);
    // DOM-driving anchors: grid + property panel are visible (E-LAYER-COMPLIANCE-01).
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 30_000});
    await page.locator('.grok-prop-panel [name="pane-Details"]').waitFor({timeout: 30_000});
  });

  // Step 1: Verify each amino acid in AlignedSequence is rendered with a distinct color.
  // Sanctioned canvas-fallback: per grok-browser/references/peptides.md "Macromolecule
  // column — monomer coloring & WebLogo header", per-cell coloring is canvas-rendered;
  // the documented assertion is the cell.renderer tag, not DOM color sampling.
  await softStep('Step 1: Verify amino acid coloring (cell.renderer === sequence)', async () => {
    const renderer = await page.evaluate(() => {
      const col = grok.shell.tv.dataFrame.col('AlignedSequence');
      return col?.getTag('cell.renderer');
    });
    expect(renderer).toBe('sequence');
  });

  // Step 2: Click the AlignedSequence column title — column becomes current and the
  // Context Panel rebuilds its accordion.
  // Column header is canvas (per project-grid-column-headers-are-canvas-not-dom):
  // the column-focus gesture is realized via `df.currentCol = df.col(...)` JS-API as
  // the sanctioned canvas-fallback. DOM-driving anchor: page.locator on the grid
  // (the gesture surface) + page.locator on pane-Peptides (the rebuild surface).
  //
  // Cold-init stabilization (per Gate B 2026-05-30-peptides-automate-02 attempt-3
  // FLAKY evidence: `.grok-prop-panel [name="pane-Details"]` was missing across the
  // full 40s Step-2 inner-poll, suggesting the property-panel infrastructure was
  // not yet wired to the TableView when Step 2 set `df.currentCol`, so the
  // column-change event was never routed to the panel and pane-Peptides never
  // mounted). Empirical MCP recon 2026-05-30:
  //   - Warm browser: pane-Details mounts at +1112ms after `onSemanticTypeDetected`;
  //     pane-Peptides mounts at +985ms after `df.currentCol = AlignedSequence`.
  //   - The wired-panel gate in the setup phase now blocks until pane-Details is
  //     present (the platform-wide pane that exists for any current column),
  //     guaranteeing the panel is wired before Step 2 runs. Step 2 therefore only
  //     needs to wait for pane-Peptides (the Macromolecule-specific pane registered
  //     by the Peptides package) to mount after the currentCol assignment.
  //   - Step 2's inner poll deadline is now 20s (vs 40s) because the panel is
  //     pre-wired by setup; +985ms warm latency + cold-init slack = comfortably
  //     inside 20s. The 20s ceiling guards against a Bio-package init regression
  //     without inflating the spec's per-run cost.
  await softStep('Step 2: Focus AlignedSequence column, wait for Context Panel rebuild', async () => {
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 15_000});
    // The wired-panel gate in the setup phase already confirmed pane-Details is
    // mounted; re-assert at the Playwright timeline so the precondition is visible
    // in the test trace.
    await page.locator('.grok-prop-panel [name="pane-Details"]').waitFor({timeout: 15_000});
    const probe = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('AlignedSequence');
      // Set `df.currentCol` first to satisfy the scenario's documented invariant
      // (clicking the column title makes it the current column).
      df.currentCol = col;
      // Then set `grok.shell.o = col` as the deterministic property-panel rebuild
      // trigger. The Context Panel responds to `grok.shell.o` (current object),
      // not directly to `df.currentCol`; setting `o` bypasses the intermediate
      // column-focus → o-binding event chain. Mirrors the sibling
      // export-invariant-map-spec.ts / export-mutation-cliffs-spec.ts pattern.
      // Self-retry: if pane-Peptides doesn't mount within the first 10s of the
      // 20s deadline, re-set `grok.shell.o = col` once to force a panel rebuild
      // (defends against a cold race where the first o-binding event fired before
      // the @panel registration was wired).
      grok.shell.o = col;
      const deadline = Date.now() + 20_000;
      const retryAt = Date.now() + 10_000;
      let peptidesPresent = false;
      let detailsPresent = false;
      let retried = false;
      while (Date.now() < deadline) {
        detailsPresent = !!document.querySelector('.grok-prop-panel [name="pane-Details"]');
        peptidesPresent = !!document.querySelector('.grok-prop-panel [name="pane-Peptides"]');
        if (detailsPresent && peptidesPresent) break;
        if (!retried && Date.now() > retryAt) {
          // One-shot re-trigger: defends against the @panel registration racing
          // the first o-binding event on a cold init.
          grok.shell.o = col;
          retried = true;
        }
        await new Promise((r) => setTimeout(r, 250));
      }
      // Settle so any lagging child content (Details key/value table, Peptides
      // synchronously-mounted button + input host) finishes mounting.
      await new Promise((r) => setTimeout(r, 1500));
      return {
        currentCol: df.currentCol?.name,
        currentO: (grok.shell.o as any)?.name,
        detailsPresent,
        peptidesPresent,
        retried,
        paneCount: document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header').length,
      };
    });
    expect(probe.currentCol).toBe('AlignedSequence');
    expect(probe.detailsPresent).toBe(true);
    expect(probe.peptidesPresent).toBe(true);
    // DOM-driving anchor on pane-Peptides (the column-context rebuild surface).
    await page.locator('.grok-prop-panel [name="pane-Peptides"]').waitFor({timeout: 15_000});
  });

  // Step 3: Confirm both Details and Peptides panes are present.
  // DOM-driving: page.locator(...).waitFor() on each pane container, scoped to the
  // property panel host (`.grok-prop-panel`) so the assertion does not accidentally
  // match a pane with the same `name=` rendered elsewhere on the page. Step 2's
  // exit predicate already gates on both being present in `.grok-prop-panel`; these
  // waits are the Playwright-visible re-assertion (the inner-evaluate's poll is not
  // visible to the Playwright timeline).
  await softStep('Step 3: Check Context Panel sections (Details + Peptides)', async () => {
    await page.locator('.grok-prop-panel [name="pane-Details"]').waitFor({timeout: 15_000});
    await page.locator('.grok-prop-panel [name="pane-Peptides"]').waitFor({timeout: 15_000});
    await expect(page.locator('.grok-prop-panel [name="pane-Details"]')).toHaveCount(1);
    await expect(page.locator('.grok-prop-panel [name="pane-Peptides"]')).toHaveCount(1);
  });

  // Step 4: Expand both panes via real DOM-driving clicks on the pane headers. This is
  // the primary E-LAYER-COMPLIANCE-01 gesture surface for the two owned context-panel
  // flows (peptides-context-panel-details, peptides-context-panel-peptides-tab).
  await softStep('Step 4: Expand Details and Peptides panes', async () => {
    // Only click headers that aren't already expanded (the accordion remembers state).
    // Scoped to `.grok-prop-panel` for the same reason as Step 3.
    const detailsHeader = page.locator('.grok-prop-panel [name="pane-Details"] .d4-accordion-pane-header');
    const peptidesHeader = page.locator('.grok-prop-panel [name="pane-Peptides"] .d4-accordion-pane-header');
    await detailsHeader.waitFor({timeout: 15_000});
    await peptidesHeader.waitFor({timeout: 15_000});
    const detailsExpanded = await detailsHeader.evaluate((el) => el.classList.contains('expanded'));
    const peptidesExpanded = await peptidesHeader.evaluate((el) => el.classList.contains('expanded'));
    if (!detailsExpanded) await detailsHeader.click();
    if (!peptidesExpanded) await peptidesHeader.click();
    // Settle so expanded pane content renders.
    await page.waitForTimeout(2500);
    // Verify expanded state.
    await expect(page.locator('.grok-prop-panel [name="pane-Details"] .d4-accordion-pane-header.expanded')).toHaveCount(1);
    await expect(page.locator('.grok-prop-panel [name="pane-Peptides"] .d4-accordion-pane-header.expanded')).toHaveCount(1);
  });

  // Step 5: Confirm each expanded pane renders its expected content.
  // Cold-stable content predicates (per project-peptides-context-panel-pane-content-predicate
  // and live recon 2026-05-30):
  //   Details — innerText contains 'Data type' AND 'Semantic type' (synchronously-rendered
  //     column metadata key/value table; NOT keyed on a lazy viewer canvas).
  //   Peptides — synchronously-mounted DOM nodes [name="button-Launch-SAR"] and
  //     [name="input-host-Activity"] (peptidesPanel SAR analyze UI surface, registered by
  //     the Peptides package); plus innerText contains 'Activity' and 'LAUNCH SAR'.
  // Per scope_reductions SR-01/SR-02 the Bioinformatics pane is out of scope.
  await softStep('Step 5: Verify Details and Peptides expanded content', async () => {
    // Details content: tolerant smoke — assert the cold-stable text signals.
    const detailsContent = page.locator('.grok-prop-panel [name="pane-Details"] .d4-accordion-pane-content');
    await detailsContent.waitFor({timeout: 15_000});
    const detailsText = (await detailsContent.innerText()) || '';
    expect(detailsText.length).toBeGreaterThan(0);
    expect(detailsText).toContain('Data type');
    expect(detailsText).toContain('Semantic type');

    // Peptides content: synchronously-mounted button + input host (NOT lazy widgets).
    const peptidesContent = page.locator('.grok-prop-panel [name="pane-Peptides"] .d4-accordion-pane-content');
    await peptidesContent.waitFor({timeout: 15_000});
    const peptidesText = (await peptidesContent.innerText()) || '';
    expect(peptidesText.length).toBeGreaterThan(0);
    // The peptidesPanel surface for SAR launch — both nodes appear synchronously.
    await expect(page.locator('.grok-prop-panel [name="pane-Peptides"] [name="button-Launch-SAR"]')).toHaveCount(1);
    await expect(page.locator('.grok-prop-panel [name="pane-Peptides"] [name="input-host-Activity"]')).toHaveCount(1);
    // Text-level cold-stable signals.
    expect(peptidesText.toUpperCase()).toContain('LAUNCH SAR');
    expect(peptidesText).toContain('Activity');
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
