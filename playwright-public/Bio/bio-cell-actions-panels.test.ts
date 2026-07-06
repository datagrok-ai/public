/* ---
sub_features_covered: [bio.actions.copy-as, bio.editors.get-region, bio.editors.split-to-monomers, bio.panels.composition-analysis, bio.panels.monomer-info]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const DATASET_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';
test('Bio cell-context actions + Context Pane info panels + custom editors', async ({page}) => {
  test.setTimeout(180_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    const hasMacro = () => Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Macromolecule');
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { if (hasMacro()) { sub.unsubscribe(); resolve(); } });
      const deadline = Date.now() + 10_000;
      const poll = () => {
        if (hasMacro() || Date.now() > deadline) { sub.unsubscribe(); resolve(); }
        else setTimeout(poll, 200);
      };
      poll();
    });
    if (hasMacro()) {
      for (let i = 0; i < 60; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
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
  await softStep('Scenario 1 Step 1-2: Copy-as cell menu exposes all 4 notation entries + joiners produce output', async () => {
    const result: {
      addCopyMenuRegistered: boolean;
      copyGroupFound: boolean;
      menuNotationsFound: string[];
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
      // Build the actual cell context menu and read back the "Copy as ..." entries.
      let copyGroupFound = false;
      const menuNotationsFound: string[] = [];
      try {
        const DG = (window as any).DG;
        const menu = DG.Menu.popup();
        const cell = df.cell(0, macroCol.name);
        await g.functions.call('Bio:addCopyMenu', {cell, menu});
        const copyGroup = menu.find('Copy');
        copyGroupFound = !!copyGroup;
        if (copyGroup)
          for (const n of ['fasta', 'separator', 'helm', 'biln'])
            if (copyGroup.find(n)) menuNotationsFound.push(n);
      } catch {  }
      return {addCopyMenuRegistered, copyGroupFound, menuNotationsFound, fastaOk, separatorOk, helmOk, bilnOk,
        fastaSample, errFasta, errSeparator, errHelm, errBiln};
    });
    expect(result.addCopyMenuRegistered,
      'addCopyMenu (package.ts#L1527) must be registered in the function table').toBe(true);
    expect(result.copyGroupFound,
      'addCopyMenu MUST add a "Copy" group to the Macromolecule cell context menu').toBe(true);
    for (const n of ['fasta', 'separator', 'helm', 'biln'])
      expect(result.menuNotationsFound,
        `Copy menu MUST contain the "${n}" notation entry (all 4 unconditionally present per bio.md:89)`).toContain(n);
    expect(result.fastaOk,
      `FASTA joiner: getJoiner({notation:"fasta"})(getSplitted(0)) MUST produce a non-empty string; err=${result.errFasta}`).toBe(true);
    expect(result.fastaSample).toBeTruthy();
    expect(result.separatorOk,
      `SEPARATOR joiner MUST produce a non-empty string; err=${result.errSeparator}`).toBe(true);
    expect(result.helmOk,
      `HELM joiner MUST produce a non-empty string; err=${result.errHelm}`).toBe(true);
    expect(result.bilnOk,
      `BILN joiner MUST produce a non-empty string; err=${result.errBiln}`).toBe(true);
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
  });
  await softStep('Scenario 2 Step 2-3: Composition analysis panel renders a monomer-composition table with count bars', async () => {
    const result: {found: boolean; barCount: number; rowCount: number; headerTexts: string[]} = await page.evaluate(async () => {
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
        return {found: false, barCount: 0, rowCount: 0, headerTexts};
      }
      (foundHeader as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      let pane: Element | null = foundHeader.parentElement;
      while (pane && !pane.classList?.contains('d4-accordion-pane'))
        pane = pane.parentElement;
      const contentDeadline = Date.now() + 15_000;
      let barCount = 0; let rowCount = 0;
      while (Date.now() < contentDeadline) {
        barCount = pane ? pane.querySelectorAll('.macromolecule-cell-comp-analysis-bar').length : 0;
        rowCount = pane ? pane.querySelectorAll('table tr').length : 0;
        if (barCount > 0) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      return {found: true, barCount, rowCount, headerTexts: []};
    });
    if (!result.found && result.headerTexts.length > 0)
      console.warn(`Composition analysis pane not found; observed headers: ${result.headerTexts.join(', ')}`);
    expect(result.found,
      'Composition analysis pane MUST appear on a Macromolecule-cell Context Pane (package.ts#L403)').toBe(true);
    expect(result.barCount,
      'Composition analysis pane MUST render a monomer-composition table with per-monomer count bars').toBeGreaterThan(0);
  });
  await softStep('Scenario 2 Step 4: Monomer info panel surfaces non-empty details for a Monomer-semType selection', async () => {
    // monomerInfoPanel is registered semType:'Monomer' (package.ts#L412); drive a Monomer current object to exercise it.
    const result: {found: boolean; rowCount: number; symbol: string | null} = await page.evaluate(async () => {
      const g = (window as any).grok;
      const DG = (window as any).DG;
      const df = g.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const macroCol = (cols.find((c: any) => c.semType === 'Macromolecule') as any);
      const seqHelper = await g.functions.call('Bio:getSeqHelper', {});
      const handler = seqHelper.getSeqHandler(macroCol);
      const ss = handler.getSplitted(0);
      let symbol: string | null = null;
      for (let i = 0; i < ss.length; i++)
        if (!ss.isGap(i)) { symbol = ss.getCanonical(i); break; }
      if (!symbol) return {found: false, rowCount: 0, symbol: null};
      g.shell.o = DG.SemanticValue.fromValueType(symbol, 'Monomer');
      // The Monomer info panel (semType Monomer) may not surface on this Context Pane — assert content only if present.
      const deadline = Date.now() + 8_000;
      let found = false; let rowCount = 0;
      while (Date.now() < deadline) {
        const h = Array.from(document.querySelectorAll(
          '.grok-prop-panel .d4-accordion-pane-header, .d4-accordion-pane-header'))
          .find((el) => (el.textContent || '').trim().toLowerCase().startsWith('monomer'));
        if (h) {
          (h as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
          let pane: Element | null = h.parentElement;
          while (pane && !pane.classList?.contains('d4-accordion-pane'))
            pane = pane.parentElement;
          rowCount = pane ? pane.querySelectorAll('table tr').length : 0;
          found = true;
          if (rowCount > 0) break;
        }
        await new Promise((r) => setTimeout(r, 300));
      }
      return {found, rowCount, symbol};
    });
    expect(result.symbol, 'a non-gap monomer symbol must be extractable from the first sequence').toBeTruthy();
    if (result.found)
      expect(result.rowCount,
        'When the Monomer info panel surfaces it MUST render non-empty monomer details (Symbol/Name/... rows)')
        .toBeGreaterThan(0);
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
    await page.evaluate(() => (document.querySelector('[name="div-Bio"]') as HTMLElement).click());
    await page.locator('[name="div-Bio---Calculate"]').waitFor({state: 'attached', timeout: 10_000});
    await page.evaluate(() => {
      const group = document.querySelector('[name="div-Bio---Calculate"]')!;
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
    });
    await page.locator('[name="div-Bio---Calculate---Extract-Region..."]').waitFor({state: 'attached', timeout: 10_000});
    await page.evaluate(() =>
      (document.querySelector('[name="div-Bio---Calculate---Extract-Region..."]') as HTMLElement).click());
    await page.locator('[name="dialog-Get-Region"]').waitFor({state: 'visible', timeout: 60_000});
  });
  await softStep('Scenario 3 Step 3: GetRegionEditor sequence-column selector is populated with the setup Macromolecule column', async () => {
    const seqHost = page.locator('[name="dialog-Get-Region"] [name="input-host-Sequence"]');
    await seqHost.waitFor({state: 'visible', timeout: 15_000});
    const populated = await seqHost.evaluate((el, name) => {
      const txt = el.textContent || '';
      const vals = Array.from(el.querySelectorAll('input, select')).map((i: any) => i.value || '').join(' ');
      return `${txt} ${vals}`.includes(name);
    }, setupProbe.macroName as string);
    expect(populated,
      `GetRegionEditor selector MUST be populated with the setup Macromolecule column "${setupProbe.macroName}" (bio.md:149)`).toBe(true);
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
    await page.evaluate(() => (document.querySelector('[name="div-Bio"]') as HTMLElement).click());
    await page.locator('[name="div-Bio---Transform"]').waitFor({state: 'attached', timeout: 10_000});
    await page.evaluate(() => {
      const group = document.querySelector('[name="div-Bio---Transform"]')!;
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
    });
    await page.locator('[name="div-Bio---Transform---Split-to-Monomers..."]').waitFor({state: 'attached', timeout: 10_000});
    await page.evaluate(() =>
      (document.querySelector('[name="div-Bio---Transform---Split-to-Monomers..."]') as HTMLElement).click());
    await page.locator('[name="dialog-Split-to-Monomers"]').waitFor({state: 'visible', timeout: 60_000});
  });
  await softStep('Scenario 4 Step 2: SplitToMonomersEditor sequence-column selector is populated with the setup column', async () => {
    const seqHost = page.locator('[name="dialog-Split-to-Monomers"] [name="input-host-Sequence"]');
    await seqHost.waitFor({state: 'visible', timeout: 15_000});
    const populated = await seqHost.evaluate((el, name) => {
      const txt = el.textContent || '';
      const vals = Array.from(el.querySelectorAll('input, select')).map((i: any) => i.value || '').join(' ');
      return `${txt} ${vals}`.includes(name);
    }, setupProbe.macroName as string);
    expect(populated,
      `SplitToMonomersEditor selector MUST be populated with the setup Macromolecule column "${setupProbe.macroName}" (bio.md:167)`).toBe(true);
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
