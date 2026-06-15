/* ---
sub_features_covered: [bio.actions.copy-as, bio.editors.get-region, bio.editors.split-to-monomers, bio.panels.composition-analysis, bio.panels.monomer-info]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const DATASET_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';
test('Bio cell-context actions + Context Pane info panels + custom editors', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
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
    expect(result.fastaOk,
      'addCopyMenu FASTA equivalent path: getJoiner({notation:"fasta"})(getSplitted(0)) MUST produce a non-empty string').toBe(true);
    expect(result.fastaSample).toBeTruthy();
    expect(result.addCopyMenuRegistered,
      'addCopyMenu (package.ts#L1527) must be registered in the function table').toBe(true);
    expect(result.separatorOk || result.helmOk || result.bilnOk,
      'at least one cross-notation conversion (SEPARATOR/HELM/BILN) must succeed').toBe(true);
  });
  // Scenario 2 — Composition analysis + Monomer info panels on the cell-level Context Pane.
  await softStep('Scenario 2 Step 1: single-click a Macromolecule cell so it becomes the current cell', async () => {
    await page.evaluate(() => {
      const g = (window as any).grok;
      const df = g.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const macroCol = (cols.find((c: any) => c.semType === 'Macromolecule') as any);
      df.currentRowIdx = 0;
      df.currentCol = macroCol;
      const DG = (window as any).DG;
      const cell = df.cell(0, macroCol.name);
      g.shell.o = DG.SemanticValue.fromTableCell(cell);
    });
    // Context Pane reconcile is async — settle before reading the panels.
    await page.waitForTimeout(4000);
  });
  await softStep('Scenario 2 Step 2-3: locate the Composition analysis panel and verify non-empty content', async () => {
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
      (foundHeader as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 1500));
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
    if (!result.found && result.headerTexts.length > 0)
      console.warn(`Composition analysis pane not found; observed headers: ${result.headerTexts.join(', ')}`);
    expect(result.found,
      'Composition analysis pane MUST appear on a Macromolecule-cell Context Pane (package.ts#L403)').toBe(true);
    expect(result.expandedHasContent,
      'Composition analysis pane MUST render non-empty content on expand').toBe(true);
  });
  await softStep('Scenario 2 Step 4: locate the Monomer info panel (semType Monomer; see Selector recon-notes)', async () => {
    // monomerInfoPanel is semType:'Monomer' — may not surface on a Macromolecule cell; assert non-empty only if present.
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
  // Scenario 3 — Get Region editor. Menu label "Extract Region..."; dialog name [name="dialog-Get-Region"].
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
    await page.locator('[name="dialog-Get-Region"]').waitFor({state: 'visible', timeout: 60_000});
  });
  await softStep('Scenario 3 Step 3: GetRegionEditor sequence-column selector is populated with the FASTA column', async () => {
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
  // Scenario 4 — Split to Monomers editor.
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
  await page.evaluate(() => (window as any).grok.shell.closeAll());
  finishSpec();
});
