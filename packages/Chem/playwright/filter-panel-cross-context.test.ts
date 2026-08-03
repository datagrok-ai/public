/* ---
sub_features_covered: [chem.sketcher, chem.sketcher.ocl]
--- */
// Sketcher wired into the FILTER PANEL: apply/clear substructure filter (GROK-14028),
// backend-switch sync to global currentSketcherType (GROK-12581/12903), close+reopen (GROK-12905).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const FILTERS = '[name="viewer-Filters"]';

async function getCurrentSketcherType(page: import('@playwright/test').Page): Promise<string> {
  return page.evaluate(() => (window as any).DG.chem.currentSketcherType);
}

// Open/reopen the molecule column's filter sketcher: click .sketch-link (empty) or canvas (filled).
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

// Switch the sketcher backend through the dialog's hamburger menu.
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
      let t0 = Date.now();
      while (Date.now() - t0 < 25000) {
        if (df.columns.toList().some((c: any) => c.semType === 'Molecule')) break;
        await sleep(300);
      }
      // Force Chem + sketcher init so the molecule column gets the substructure filter.
      DG.chem.currentSketcherType = 'OpenChemLib';
      const probe = new DG.chem.Sketcher();
      const pd = ui.dialog('init').add(probe.root); pd.show();
      t0 = Date.now();
      while (Date.now() - t0 < 20000) { if (probe.sketcher?.isInitialized) break; await sleep(300); }
      pd.close();
      (window as any).__sk_err = [];
      const orig = console.error;
      console.error = function(...a: any[]) { (window as any).__sk_err.push(a.map((x: any) => String(x)).join(' ')); orig.apply(console, a as any); };
      // retry getFiltersGroup until the substructure sketch-link appears
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
    // GROK-14028: exercise the Filter Panel's Reset (.chem-clear-sketcher-button, hover-revealed, click via DOM).
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
    // reopen to confirm Reset cleared the persisted input line
    await openFilterSketcher(page);
    const inputVal = await page.locator('.d4-dialog input[placeholder*="SMILES" i]').inputValue();
    console.log(`[fp] after Reset: input="${inputVal}" trueCount=${trueCount}/${total}`);
    expect(inputVal, 'GROK-14028: Reset must clear the sketcher input line').toBe('');
    await closeDialog(page);
  });

  await softStep('Block B: switch backend in the filter sketcher → global propagates (GROK-12581/12903)', async () => {
    // ChemDraw/Marvin are proprietary backends absent from the public build, so
    // this exercises bidirectional propagation across the two public backends:
    // OpenChemLib (Chem) and Ketcher (KetcherSketcher).
    await openFilterSketcher(page);
    await switchBackendInDialog(page, 'Ketcher');
    expect(await getCurrentSketcherType(page),
      'switching in the Filter Panel sketcher must update the shared current sketcher (read by hamburger/Context Pane)').toBe('Ketcher');
    await switchBackendInDialog(page, 'OpenChemLib');
    expect(await getCurrentSketcherType(page), 'second switch also propagates').toBe('OpenChemLib');
  });

  await softStep('Block C: switch back to OpenChemLib, close, reopen filter sketcher (GROK-12905)', async () => {
    await switchBackendInDialog(page, 'OpenChemLib');
    expect(await getCurrentSketcherType(page)).toBe('OpenChemLib');
    await closeDialog(page);
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

  finishSpec();
});
