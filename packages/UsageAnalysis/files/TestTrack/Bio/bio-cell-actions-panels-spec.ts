/* ---
sub_features_covered:
  - bio.actions.copy-as
  - bio.editors.get-region
  - bio.editors.split-to-monomers
  - bio.panels.composition-analysis
  - bio.panels.monomer-info
--- */
//   related_bugs: [] (no curated bug-library/bio.yaml entry has any of the
//   Hypothesis category: test-bug (per §"Hypothesis protocol" round 1) —
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
test.use(specTestOptions);
const DATASET_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';
test('Bio cell-context actions + Context Pane info panels + custom editors', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  // bio.md "GROK-18616 entry-path-class invariant" — readCsv triggers the
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
  }, DATASET_PATH);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch {  }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });
  await page.evaluate(async () => {
    const candidates = [
      ['Bio:getRegionTopMenu', 'Bio:getRegion', 'Bio:extractRegionTopMenu', 'Bio:extractRegion'],
      ['Bio:splitToMonomersTopMenu', 'Bio:splitToMonomers'],
      ['Bio:compositionAnalysisWidget'],
      ['Bio:monomerInfoPanel'],
      ['Bio:addCopyMenu'],
    ];
    const findAny = (names: string[]): boolean => {
      for (const n of names) {
        try {
          if ((grok as any).functions.find && (grok as any).functions.find(n)) return true;
        } catch {  }
      }
      return false;
    };
    const deadline = Date.now() + 15_000;
    while (Date.now() < deadline) {
      if (candidates.every(findAny)) return;
      await new Promise((r) => setTimeout(r, 300));
    }
    await new Promise((r) => setTimeout(r, 1500));
  });
  await page.waitForTimeout(2000);
  const setupProbe = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
    return {
      hasMacromoleculeCol: !!macro,
      macroName: macro?.name ?? null,
      macroSemType: macro?.semType ?? null,
      rowCount: df.rowCount,
    };
  });
  expect(setupProbe.hasMacromoleculeCol,
    'atlas bio.detector contract: readCsv path MUST classify a Macromolecule column synchronously').toBe(true);
  expect(setupProbe.macroSemType).toBe('Macromolecule');
  expect(setupProbe.rowCount,
    'scenario .md Setup: ≥ 5 rows so cell-level operations have non-degenerate input').toBeGreaterThanOrEqual(5);
  await softStep('Scenario 1 Step 1-2: addCopyMenu function registered + 4 notation paths reachable', async () => {
    // methods that do NOT exist on the real API (test-bug, fabricated method
    const result: {
      addCopyMenuRegistered: boolean;
      fastaOk: boolean;
      separatorOk: boolean;
      helmOk: boolean;
      bilnOk: boolean;
      fastaSample: string | null;
      errFasta: string | null;
      errSeparator: string | null;
      errHelm: string | null;
      errBiln: string | null;
    } = await page.evaluate(async () => {
      const g = (window as any).grok;
      let addCopyMenuRegistered = false;
      try {
        const find = (g as any).functions && (g as any).functions.find;
        for (const candidate of ['Bio:addCopyMenu', 'addCopyMenu']) {
          try { if (find && find(candidate)) { addCopyMenuRegistered = true; break; } } catch {  }
        }
      } catch {  }
      const df = g.shell.tv.dataFrame;
      const macroCol = Array.from({length: df.columns.length}, (_: unknown, i: number) =>
        df.columns.byIndex(i)).find((c: any) => c.semType === 'Macromolecule') as any;
      let fastaOk = false; let separatorOk = false; let helmOk = false; let bilnOk = false;
      let fastaSample: string | null = null;
      let errFasta: string | null = null;
      let errSeparator: string | null = null;
      let errHelm: string | null = null;
      let errBiln: string | null = null;
      try {
        const seqHelper = await (g as any).functions.call('Bio:getSeqHelper', {});
        const handler = seqHelper.getSeqHandler(macroCol);
        try {
          const srcSS = handler.getSplitted(0);
          const joiner = handler.getJoiner({notation: 'fasta'});
          const fasta = joiner(srcSS);
          fastaOk = typeof fasta === 'string' && fasta.length > 0;
          fastaSample = fastaOk ? fasta.slice(0, 32) : null;
        } catch (e: any) { errFasta = (e && (e.message || String(e))) || 'unknown'; }
        const defaultSeparator = '-';
        for (const tgt of ['separator', 'helm', 'biln']) {
          try {
            const srcSS = handler.getSplitted(0);
            const joiner = handler.getJoiner({notation: tgt, separator: tgt === 'separator' ? defaultSeparator : undefined});
            const converted = joiner(srcSS);
            const ok = typeof converted === 'string' && converted.length > 0;
            if (tgt === 'separator') separatorOk = ok;
            if (tgt === 'helm') helmOk = ok;
            if (tgt === 'biln') bilnOk = ok;
          } catch (e: any) {
            const msg = (e && (e.message || String(e))) || 'unknown';
            if (tgt === 'separator') errSeparator = msg;
            if (tgt === 'helm') errHelm = msg;
            if (tgt === 'biln') errBiln = msg;
          }
        }
      } catch {  }
      return {addCopyMenuRegistered, fastaOk, separatorOk, helmOk, bilnOk, fastaSample,
        errFasta, errSeparator, errHelm, errBiln};
    });
    if (!result.fastaOk && result.errFasta)
      console.warn(`Scenario 1 FASTA joiner error: ${result.errFasta}`);
    if (!result.separatorOk && result.errSeparator)
      console.warn(`Scenario 1 SEPARATOR joiner error: ${result.errSeparator}`);
    if (!result.helmOk && result.errHelm)
      console.warn(`Scenario 1 HELM joiner error: ${result.errHelm}`);
    if (!result.bilnOk && result.errBiln)
      console.warn(`Scenario 1 BILN joiner error: ${result.errBiln}`);
    // Primary assertion (scenario .md Step 3): the FASTA Copy-as path
    // (matching source notation) MUST produce a non-empty round-trip string
    // — this is the closest JS-API surrogate for "after clicking Copy as
    // FASTA, the clipboard contains the FASTA-form representation".
    expect(result.fastaOk,
      'addCopyMenu FASTA equivalent path: getJoiner({notation:"fasta"})(getSplitted(0)) MUST produce a non-empty string').toBe(true);
    expect(result.fastaSample).toBeTruthy();
    // Function-registry presence: addCopyMenu is registered (so the cell
    // menu would surface the four entries on real-OS right-click).
    expect(result.addCopyMenuRegistered,
      'addCopyMenu (package.ts#L1527) must be registered in the function table').toBe(true);
    // Cross-notation reachability — at least one of separator/helm/biln
    // conversions must succeed. Exact subset depends on the source sequence
    // alphabet (e.g. BILN may decline on a pure-FASTA-peptide source); the
    // integration-level assertion (per scenario .md Expected) is that the
    // cross-notation conversions are reachable, not that all four notations
    // accept every source unit.
    expect(result.separatorOk || result.helmOk || result.bilnOk,
      'at least one cross-notation conversion (SEPARATOR/HELM/BILN) must succeed').toBe(true);
  });
  // -------------------------------------------------------------------
  // Scenario 2 — Composition analysis + Monomer info panels on Context Pane.
  //
  // Trigger: set grok.shell.o to a SemanticValue wrapping the first
  // Macromolecule cell. This is the cell-level Context Pane entrypoint
  // — distinct from the column-level entrypoint
  // (grok.shell.o = df.col(...)) documented in bio.md:499-530 which is
  // covered by sibling scenarios.
  //
  // Asserted panels:
  //   • Composition analysis (package.ts#L403, semType Macromolecule)
  //     — MUST appear on a Macromolecule cell selection. Spec asserts
  //     the panel surfaces with non-empty content.
  //   • Monomer info (package.ts#L412, semType Monomer) — the decorator
  //     restricts this panel to semType 'Monomer' cells, not
  //     'Macromolecule'. Spec asserts presence-or-absence-with-rationale
  //     per scenario .md Step 4 Expected ("…or the cell-level summary for
  //     the parent Macromolecule cell — the panel must surface non-empty
  //     content").
  // -------------------------------------------------------------------
  await softStep('Scenario 2 Step 1: single-click a Macromolecule cell so it becomes the current cell', async () => {
    // Round-1 retry empirical correction: the prior author called
    // DG.SemanticValue.fromTableCell(g.shell.tv.grid.cell(name, idx)). The
    // grid.cell(...) accessor returns a GridCell; fromTableCell expects a
    // DataFrame Cell (per public/js-api/src/grid.ts:1376 — `static fromTableCell(cell: Cell)`).
    // The arg-type mismatch could surface as either a null SemanticValue or
    // a no-op subject set, in either case leaving the Context Pane WITHOUT
    // a cell-scope subject and the Composition pane never surfaces — the
    // most plausible root cause of B-RUN-PASS at Scenario 2 Step 2.
    //
    // Same-paradigm tactical fix: use the DataFrame Cell accessor
    // (df.cell(rowIdx, colName)) which returns the right type. Per
    // composition-analysis-widget.ts:18 the widget reads val.cell.column /
    // val.cell.rowIndex — both available on a DataFrame Cell.
    await page.evaluate(() => {
      const g = (window as any).grok;
      const df = g.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const macroCol = (cols.find((c: any) => c.semType === 'Macromolecule') as any);
      df.currentRowIdx = 0;
      df.currentCol = macroCol;
      // Set the cell-level subject so the Context Pane surfaces cell-scope
      // panels (the per-cell entrypoint that compositionAnalysisWidget +
      // monomerInfoPanel register against). Pass a DataFrame Cell (not a
      // GridCell) — fromTableCell signature expects a DataFrame Cell.
      const DG = (window as any).DG;
      const cell = df.cell(0, macroCol.name);
      g.shell.o = DG.SemanticValue.fromTableCell(cell);
    });
    // Give the property-panel framework a moment to populate. Context
    // Pane reconcile is asynchronous (panel widgets are invoked via the
    // platform's @panel registration pipeline after grok.shell.o changes).
    // The cold-start window can be wider than the prior 2s — bump to 4s
    // for resilience (same paradigm as composition-analysis-spec.ts Step 4
    // settle for the property-grid bind).
    await page.waitForTimeout(4000);
  });
  await softStep('Scenario 2 Step 2-3: locate the Composition analysis panel and verify non-empty content', async () => {
    // Poll up to 60s for the property panel to show a section whose
    // header text matches "Composition analysis" (the @grok.decorators.panel
    // name in package.ts#L403). The panel may surface inside an accordion
    // section that requires expand-on-click to render content. Round-1
    // retry tactical reinforcement: bumped poll window from 30s → 60s to
    // accommodate the cold-start @panel registration race that may have
    // contributed to the prior Gate B B-STAB-01 (the prior author already
    // gave Scenarios 3/4 a 60s dialog tolerance; symmetric for the
    // Context Pane panel surface).
    const result: {found: boolean; expandedHasContent: boolean; headerTexts: string[]} = await page.evaluate(async () => {
      const deadline = Date.now() + 60_000;
      let foundHeader: Element | null = null;
      const matchHeader = (h: Element) => {
        const t = (h.textContent || '').trim();
        return t === 'Composition analysis' || t.toLowerCase().startsWith('composition');
      };
      while (Date.now() < deadline) {
        const headers = Array.from(document.querySelectorAll(
          '.grok-prop-panel .d4-accordion-pane-header, .d4-accordion-pane-header'));
        const h = headers.find(matchHeader);
        if (h) { foundHeader = h; break; }
        await new Promise((r) => setTimeout(r, 300));
      }
      if (!foundHeader) {
        const headerTexts = Array.from(document.querySelectorAll(
          '.grok-prop-panel .d4-accordion-pane-header'))
          .map((h) => (h.textContent || '').trim())
          .filter((t) => t.length > 0);
        return {found: false, expandedHasContent: false, headerTexts};
      }
      // Expand the pane (click header) if it's collapsed.
      (foundHeader as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 1500));
      // Walk up to the .d4-accordion-pane ancestor and look inside its
      // content region for non-empty content (any rendered child element
      // — the composition table, a label, etc.).
      let pane: Element | null = foundHeader.parentElement;
      while (pane && !pane.classList?.contains('d4-accordion-pane'))
        pane = pane.parentElement;
      const content = pane?.querySelector(':scope > :not(.d4-accordion-pane-header)');
      const text = (content?.textContent || '').trim();
      const childCount = content ? content.children.length : 0;
      const expandedHasContent = !!content && (text.length > 0 || childCount > 0);
      const headerTexts: string[] = [];
      return {found: true, expandedHasContent, headerTexts};
    });
    // Diagnostic: print the headers observed when the panel was not found
    // — surfaces the actual property-panel section names so a retry can
    // refine the header-match predicate.
    if (!result.found && result.headerTexts.length > 0)
      console.warn(`Composition analysis pane not found; observed headers: ${result.headerTexts.join(', ')}`);
    expect(result.found,
      'Composition analysis pane MUST appear on a Macromolecule-cell Context Pane (package.ts#L403)').toBe(true);
    expect(result.expandedHasContent,
      'Composition analysis pane MUST render non-empty content on expand').toBe(true);
  });
  await softStep('Scenario 2 Step 4: locate the Monomer info panel (semType Monomer; see Selector recon-notes)', async () => {
    // monomerInfoPanel (package.ts#L412) is decorated with semType:
    // 'Monomer' — the param decorator restricts the panel to Monomer
    // cells (individual per-position monomer cells produced by Split to
    // Monomers), not the parent Macromolecule sequence cell. On a
    // Macromolecule cell selection the panel may or may not surface
    // depending on the platform's @panel registration filter on
    // SemanticValue.semType. The scenario .md Step 4 Expected
    // accommodates both cases ("…or the cell-level summary for the
    // parent Macromolecule cell — the panel must surface non-empty
    // content"). We assert non-empty content IF the panel surfaces; we
    // do NOT fail-hard when it does not (the scenario's primary
    // assertion is "no error balloon appears in either panel").
    const result: {found: boolean; hasContent: boolean} = await page.evaluate(async () => {
      const headers = Array.from(document.querySelectorAll(
        '.grok-prop-panel .d4-accordion-pane-header, .d4-accordion-pane-header'));
      const h = headers.find((el) => {
        const t = (el.textContent || '').trim();
        return t === 'Monomer' || t.toLowerCase().startsWith('monomer');
      });
      if (!h) return {found: false, hasContent: false};
      (h as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 1500));
      let pane: Element | null = h.parentElement;
      while (pane && !pane.classList?.contains('d4-accordion-pane'))
        pane = pane.parentElement;
      const content = pane?.querySelector(':scope > :not(.d4-accordion-pane-header)');
      const text = (content?.textContent || '').trim();
      const childCount = content ? content.children.length : 0;
      const hasContent = !!content && (text.length > 0 || childCount > 0);
      return {found: true, hasContent};
    });
    // Soft predicate: when the panel surfaces, it MUST be non-empty. When
    // it does not surface (semType Monomer filter on a Macromolecule
    // cell), the assertion is the no-error-balloon predicate at the
    // end of Scenario 2.
    if (result.found) {
      expect(result.hasContent,
        'When the Monomer pane surfaces on a Macromolecule-cell selection it MUST render non-empty content').toBe(true);
    }
  });
  await softStep('Scenario 2: no balloon error fired by either context-pane info panel', async () => {
    const balloonError = await page.evaluate(() => {
      const errors = document.querySelectorAll('.d4-balloon.error, .grok-balloon-error');
      return errors.length;
    });
    expect(balloonError, 'no error balloon must be emitted by either Context Pane panel').toBe(0);
  });
  // -------------------------------------------------------------------
  // Scenario 3 — Get Region editor (GetRegionEditor; package.ts#L213).
  //
  // Atlas sub_feature id: bio.editors.get-region.
  // Live menu label: "Extract Region..." (bio.md:334-348; Bio CLAUDE.md
  // top-menu table). The dialog title is "Get Region" — the dialog name
  // [name="dialog-Get-Region"] persists across the menu-label rename.
  // -------------------------------------------------------------------
  await softStep('Scenario 3 Step 1-2: click Bio > Calculate > Extract Region — GetRegionEditor dialog opens', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      const group = document.querySelector('[name="div-Bio---Calculate"]');
      if (!group) throw new Error('[name="div-Bio---Calculate"] group anchor not found');
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      const leaf = document.querySelector('[name="div-Bio---Calculate---Extract-Region..."]');
      if (!leaf) throw new Error('[name="div-Bio---Calculate---Extract-Region..."] leaf not found');
      (leaf as HTMLElement).click();
    });
    // Cold-start dialog tolerance (matches analyze-spec.ts L180).
    await page.locator('[name="dialog-Get-Region"]').waitFor({state: 'visible', timeout: 60_000});
  });
  await softStep('Scenario 3 Step 3: GetRegionEditor sequence-column selector is populated with the FASTA column', async () => {
    // The custom editor exposes a sequence-column input bound to the
    // active TableView's Macromolecule column. bio.md:344 names the
    // selector pair [name="input-host-Sequence"] / [name="input-Sequence"].
    // Asserts:
    //   • the sequence-column input is attached to the dialog DOM (input
    //     widgets render under [name="input-host-…"]);
    //   • the dialog is visible (visible-state predicate from waitFor above).
    const seqInputPresent = await page.locator(
      '[name="dialog-Get-Region"] [name="input-host-Sequence"]').count();
    expect(seqInputPresent,
      'GetRegionEditor MUST render a sequence-column selector input (bio.md:344)').toBeGreaterThan(0);
  });
  await softStep('Scenario 3 Step 4: Cancel via Escape closes the dialog with no balloon error', async () => {
    await page.keyboard.press('Escape');
    await page.locator('[name="dialog-Get-Region"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(balloonError,
      'Cancel must close GetRegionEditor with no error balloon').toBe(0);
  });
  // -------------------------------------------------------------------
  // Scenario 4 — Split to Monomers editor (SplitToMonomersEditor;
  // package.ts#L226).
  //
  // Atlas sub_feature id: bio.editors.split-to-monomers.
  // Menu label is stable: "Split to Monomers..." (bio.md:395-405).
  // -------------------------------------------------------------------
  await softStep('Scenario 4 Step 1: click Bio > Transform > Split to Monomers — SplitToMonomersEditor dialog opens', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      const group = document.querySelector('[name="div-Bio---Transform"]');
      if (!group) throw new Error('[name="div-Bio---Transform"] group anchor not found');
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      const leaf = document.querySelector('[name="div-Bio---Transform---Split-to-Monomers..."]');
      if (!leaf) throw new Error('[name="div-Bio---Transform---Split-to-Monomers..."] leaf not found');
      (leaf as HTMLElement).click();
    });
    await page.locator('[name="dialog-Split-to-Monomers"]').waitFor({state: 'visible', timeout: 60_000});
  });
  await softStep('Scenario 4 Step 2: SplitToMonomersEditor sequence-column selector is populated', async () => {
    const seqInputPresent = await page.locator(
      '[name="dialog-Split-to-Monomers"] [name="input-host-Sequence"]').count();
    expect(seqInputPresent,
      'SplitToMonomersEditor MUST render a sequence-column selector input (bio.md:403)').toBeGreaterThan(0);
  });
  await softStep('Scenario 4 Step 3: Cancel via Escape closes the dialog with no balloon error', async () => {
    await page.keyboard.press('Escape');
    await page.locator('[name="dialog-Split-to-Monomers"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(balloonError,
      'Cancel must close SplitToMonomersEditor with no error balloon').toBe(0);
  });
  // Cleanup: close all views (free TableView state).
  await page.evaluate(() => (window as any).grok.shell.closeAll());
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
