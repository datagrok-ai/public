/* ---
sub_features_covered: [chem.sketcher, chem.sketcher.ocl]
--- */
// Frontmatter extraction:
//   target_layer: playwright
//   coverage_type: regression
//   sub_features_covered: [chem.sketcher, chem.sketcher.ocl]
//   ui_coverage_responsibility: [chem-filter-substructure-open, chem-filter-substructure-apply,
//     chem-filter-substructure-clear, chem-filter-backend-switch-sync, chem-filter-sketcher-reopen]
//   related_bugs: [GROK-12581, GROK-12903, GROK-12905, GROK-14028]
//
// Cross-context sketcher bugs where the sketcher is wired into the FILTER PANEL (not the
// standalone widget). All steps are real UI gestures (filters.md reference):
//   open Filters panel → molecule column's substructure filter card has a `.sketch-link`
//   → click it → the SAME sketcher dialog opens (SMILES input + hamburger ≡ backend switcher).
//
//   Block A  apply substructure filter (type SMILES → Enter → OK → rows filtered)  baseline
//            then clear → input empties + rows restored ............................ GROK-14028
//   Block B  switch backend via the filter sketcher's hamburger → the global
//            currentSketcherType (what the column hamburger / Context Pane read) follows  GROK-12581 / GROK-12903
//   Block C  switch backend then back to OpenChemLib, close + reopen the filter
//            sketcher → it reopens ................................................. GROK-12905
//
// Recon (chrome-devtools MCP @ dev.datagrok.ai): filter panel `[name="viewer-Filters"]`,
// cards `.d4-filter`, structure card has `.sketch-link` (class `chem-sketch-link sketch-link`)
// that opens a `.d4-dialog` with `input[placeholder*="SMILES"]` + `.fa-bars`. Applying
// `c1ccccc1` on smiles-50.csv drops filter.trueCount 50→47. Backend switch in the dialog is
// in-place on one widget (stable, unlike fresh-widget churn). NB: getFiltersGroup() must run
// AFTER semType detection or the molecule column won't get the substructure filter (filters.md).
//
// Paired scenario: filter-panel-cross-context.md
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

// No `storageState`: `loginToDatagrok(page)` inside the test body is the canonical auth path
// per spec-login.ts (DATAGROK_AUTH_TOKEN injection by `grok test`). A `storageState: 'auth.json'`
// directive here would shadow that path — Playwright would resolve the file at module-load
// against `public/packages/UsageAnalysis/` and ENOENT the whole spec before any test() runs.
test.use(specTestOptions);

const FILTERS = '[name="viewer-Filters"]';

async function getCurrentSketcherType(page: import('@playwright/test').Page): Promise<string> {
  return page.evaluate(() => (window as any).DG.chem.currentSketcherType);
}

// Open (or reopen) the molecule column's filter sketcher dialog: click the `.sketch-link`
// (empty state) or the rendered canvas (filled state).
async function openFilterSketcher(page: import('@playwright/test').Page): Promise<void> {
  const link = page.locator(`${FILTERS} .d4-filter .sketch-link`);
  if (await link.count() > 0)
    await link.first().click();
  else
    await page.locator(`${FILTERS} .d4-filter .chem-external-sketcher-canvas`).first().click();
  await page.locator('.d4-dialog input[placeholder*="SMILES" i]').waitFor({timeout: 15000});
}

async function closeDialog(page: import('@playwright/test').Page): Promise<void> {
  const cancel = page.locator('.d4-dialog [name="button-CANCEL"], .d4-dialog .ui-btn-cancel');
  if (await cancel.count() > 0) await cancel.first().click();
  else {
    const ok = page.locator('.d4-dialog [name="button-OK"], .d4-dialog .ui-btn-ok');
    if (await ok.count() > 0) await ok.first().click();
  }
  await page.waitForTimeout(800);
}

// Switch the sketcher backend through the dialog's hamburger ≡ menu (real DOM clicks).
async function switchBackendInDialog(page: import('@playwright/test').Page, name: string): Promise<void> {
  await page.locator('.d4-dialog .fa-bars').first().click();
  await page.waitForTimeout(800);
  await page.locator('.d4-menu-item-label').filter({hasText: new RegExp(`^${name}$`)}).first().click();
  await page.waitForFunction((b) => (window as any).DG?.chem?.currentSketcherType === b, name, {timeout: 60000});
  await page.waitForTimeout(2500);
}

test('Chem: Filter Panel sketcher — apply / clear / backend-switch sync / reopen', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);
  await page.waitForTimeout(3000);

  let total = 0;
  await softStep('Open chem table + Filter Panel (substructure filter on molecule column)', async () => {
    const res = await page.evaluate(async () => {
      const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/smiles-50.csv');
      const tv = grok.shell.addTableView(df);
      // wait for semType detection
      let t0 = Date.now();
      while (Date.now() - t0 < 25000) {
        if (df.columns.toList().some((c: any) => c.semType === 'Molecule')) break;
        await sleep(300);
      }
      // Force Chem + sketcher backends to initialize in this fresh session — otherwise the
      // molecule column won't get the substructure (sketcher) filter when the panel opens.
      DG.chem.currentSketcherType = 'OpenChemLib';
      const probe = new DG.chem.Sketcher();
      const pd = ui.dialog('init').add(probe.root); pd.show();
      t0 = Date.now();
      while (Date.now() - t0 < 20000) { if (probe.sketcher?.isInitialized) break; await sleep(300); }
      pd.close();
      (window as any).__sk_err = [];
      const orig = console.error;
      console.error = function(...a: any[]) { (window as any).__sk_err.push(a.map((x: any) => String(x)).join(' ')); orig.apply(console, a as any); };
      // open filters; retry getFiltersGroup until the substructure sketch-link appears
      let found = false;
      for (let i = 0; i < 4 && !found; i++) {
        tv.getFiltersGroup();
        const w = Date.now();
        while (Date.now() - w < 9000) {
          if (document.querySelector('[name="viewer-Filters"] .d4-filter .sketch-link')) { found = true; break; }
          await sleep(400);
        }
      }
      return {total: df.rowCount, found};
    });
    expect(res.found, 'molecule column must get a substructure (sketch) filter in the panel').toBe(true);
    total = res.total;
    expect(total, 'dataset loaded').toBeGreaterThan(0);
  });

  await softStep('Block A: apply substructure filter (SMILES → Enter → OK) filters rows', async () => {
    await openFilterSketcher(page);
    const smiles = page.locator('.d4-dialog input[placeholder*="SMILES" i]');
    await smiles.click();
    await page.keyboard.type('c1ccccc1', {delay: 30});
    await page.keyboard.press('Enter');
    await page.waitForTimeout(2500);
    const ok = page.locator('.d4-dialog [name="button-OK"], .d4-dialog .ui-btn-ok');
    if (await ok.count() > 0) await ok.first().click();
    await page.waitForTimeout(3000);
    const trueCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    console.log(`[fp] applied benzene: total=${total} trueCount=${trueCount}`);
    expect(trueCount, 'substructure filter must reduce the visible rows').toBeLessThan(total);
    expect(trueCount, 'benzene matches most rows (not zero)').toBeGreaterThan(0);
  });

  await softStep('Block A: Reset on the Filter Panel clears input + restores rows (GROK-14028)', async () => {
    // GROK-14028 path per the scenario .md: the gesture is the Filter Panel's Reset
    // (the substructure card's `.chem-clear-sketcher-button`), NOT a manual input wipe —
    // the bug was that Reset failed to clear the persisted sketcher input. Manually emptying
    // the input would assert a tautology and never exercise the fix. The button is
    // hover-revealed, so click it via DOM (mirrors chem-grok-14028-spec.ts).
    const clicked = await page.evaluate(async () => {
      const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));
      let btn: HTMLElement | null = null;
      for (let i = 0; i < 30; i++) {
        btn = document.querySelector('[name="viewer-Filters"] .d4-filter .chem-clear-sketcher-button');
        if (btn) break;
        await sleep(500);
      }
      if (!btn) return false;
      btn.click();
      await sleep(2500);
      return true;
    });
    expect(clicked, 'GROK-14028: Reset (.chem-clear-sketcher-button) must be present on the Filter Panel after applying').toBe(true);
    const trueCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(trueCount, 'GROK-14028: rows restored after Reset on the Filter Panel').toBe(total);
    // Reopen the sketcher to confirm Reset also cleared its persisted input line.
    await openFilterSketcher(page);
    const inputVal = await page.locator('.d4-dialog input[placeholder*="SMILES" i]').inputValue();
    console.log(`[fp] after Reset: input="${inputVal}" trueCount=${trueCount}/${total}`);
    expect(inputVal, 'GROK-14028: Reset must clear the sketcher input line').toBe('');
    await closeDialog(page);
  });

  await softStep('Block B: switch backend in the filter sketcher → global propagates (GROK-12581/12903)', async () => {
    await openFilterSketcher(page);
    await switchBackendInDialog(page, 'Ketcher');
    expect(await getCurrentSketcherType(page),
      'switching in the Filter Panel sketcher must update the shared current sketcher (read by hamburger/Context Pane)').toBe('Ketcher');
    await switchBackendInDialog(page, 'ChemDraw');
    expect(await getCurrentSketcherType(page), 'second switch also propagates').toBe('ChemDraw');
  });

  await softStep('Block C: switch back to OpenChemLib, close, reopen filter sketcher (GROK-12905)', async () => {
    await switchBackendInDialog(page, 'OpenChemLib');
    expect(await getCurrentSketcherType(page)).toBe('OpenChemLib');
    await closeDialog(page);
    // reopen from the filter panel — must open again
    await openFilterSketcher(page);
    const reopened = await page.locator('.d4-dialog input[placeholder*="SMILES" i]').count();
    expect(reopened, 'GROK-12905: sketcher must reopen from the Filter Panel after backend churn').toBeGreaterThan(0);
    await closeDialog(page);
  });

  await softStep('No sketcher console errors during the filter-panel walk', async () => {
    const errs = await page.evaluate(() => ((window as any).__sk_err ?? []) as string[]);
    const sk = errs.filter((e) => /sketcher|search pattern cannot be set|substructure|chem-filter/i.test(e));
    expect(sk.length, `sketcher console errors: ${JSON.stringify(sk.slice(0, 5))}`).toBe(0);
  });

  await page.evaluate(() => {
    try { grok.userSettings.add('sketcher', 'selected', 'OpenChemLib'); } catch (e) {}
    grok.shell.closeAll();
  }).catch(() => {});

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
