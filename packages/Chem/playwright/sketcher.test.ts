/* ---
sub_features_covered: [chem.actions.copy-as, chem.actions.copy-molfile-v2000, chem.actions.copy-smiles, chem.sketcher, chem.sketcher.cell-editor, chem.sketcher.ocl]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Chem: Sketcher Favorites + Recent + Copy as SMILES/MOLBLOCK + input round-trip', async ({page}) => {
  test.setTimeout(120_000);

  const CYCLOHEXANE = 'C1CCCCC1';
  const ETHANOL = 'CCO';
  const FAV_KEYS = ['chem-molecule-favorites', 'smiles-favorites'];
  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));

  await page.context().grantPermissions(['clipboard-read', 'clipboard-write']);

  const smilesInput = page.locator('.d4-dialog input[placeholder*="SMILES" i]');

  const openHamburger = async (): Promise<void> => {
    const found = await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');
      const h = dialog?.querySelector('.fa-bars, [name="icon-font-icon-menu"], .grok-icon.fa-bars') as HTMLElement | null;
      if (!h) return false;
      h.click();
      return true;
    });
    expect(found, 'sketcher hamburger icon not found').toBe(true);
  };

  const waitForLabel = async (re: RegExp): Promise<void> => {
    await page.waitForFunction((src) => {
      const rx = new RegExp(src, 'i');
      return Array.from(document.querySelectorAll('.d4-menu-item-label')).some((l) => rx.test(l.textContent || ''));
    }, re.source, {timeout: 10_000});
  };

  const clickLabel = async (re: RegExp): Promise<void> => {
    await waitForLabel(re);
    await page.evaluate((src) => {
      const rx = new RegExp(src, 'i');
      const it = Array.from(document.querySelectorAll('.d4-menu-item-label')).find((l) => rx.test(l.textContent || ''));
      (it!.closest('.d4-menu-item') as HTMLElement).click();
    }, re.source);
  };

  const menuLabels = async (): Promise<string[]> =>
    page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-item-label')).map((l) => l.textContent!.trim()));

  const clipboard = async (): Promise<string> => (await page.evaluate(() => navigator.clipboard.readText())).trim();

  const dismissMenu = async (): Promise<void> => {
    await page.evaluate(() => {
      const t = document.querySelector('.d4-dialog .d4-dialog-title, .d4-dialog .d4-dialog-header, .d4-dialog') as HTMLElement | null;
      t?.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    });
  };

  const typeSmiles = async (smiles: string): Promise<void> => {
    await smilesInput.click();
    await smilesInput.fill(smiles);
    await smilesInput.press('Enter');
    await expect.poll(() => smilesInput.inputValue(), {timeout: 10_000}).toBe(smiles);
  };

  await loginToDatagrok(page);
  await page.waitForTimeout(3000);

  // navigator.clipboard is only defined in a secure context (https:// or localhost). CI runs
  // against a plain-http dev host, so the Copy/Paste round-trip steps (8/8b/9/9b) can't read the
  // clipboard there — grantPermissions above cannot conjure the API. Detect availability once and
  // skip only those steps; the menu-presence assertions (Step 3) still run everywhere.
  const clipboardOk = await page.evaluate(() => {
    try { return !!(navigator.clipboard && typeof navigator.clipboard.readText === 'function'); }
    catch { return false; }
  });

  await softStep('Step 1: Open smiles-50.csv + Molecule semType ready', async () => {
    await page.evaluate(async (favKeys) => {
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      favKeys.concat(['chem-molecule-recent', 'smiles-recent']).forEach((k) => localStorage.removeItem(k));
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/smiles-50.csv');
      grok.shell.addTableView(df);
      (window as any).__sk_errors = [];
      const orig = console.error;
      console.error = function(...args: any[]) {
        (window as any).__sk_errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
    }, FAV_KEYS);
    await waitForChemMenu(page);
  });

  await softStep('Step 2: Open sketcher cell editor (JS API)', async () => {
    await page.evaluate(async () => {
      const molColName = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule')?.name ?? 'canonical_smiles';
      const smiles = grok.shell.t.col(molColName).get(0);
      const sketcher = grok.chem.sketcher(() => null, smiles);
      ui.dialog('Sketcher').add(sketcher).show();
    });
    await page.locator('.d4-dialog').waitFor({timeout: 10_000});
    await smilesInput.waitFor({timeout: 30_000});
  });

  await softStep('Step 3: Hamburger menu exposes Favorites + Recent + Copy as SMILES/MOLBLOCK', async () => {
    await openHamburger();
    await waitForLabel(/Copy as SMILES/);
    const labels = await menuLabels();
    expect(labels.some((l) => /Favorites/i.test(l)), `Favorites missing: ${JSON.stringify(labels.slice(0, 20))}`).toBe(true);
    expect(labels.some((l) => /Recent/i.test(l)), `Recent missing: ${JSON.stringify(labels.slice(0, 20))}`).toBe(true);
    expect(labels.some((l) => /Copy as SMILES/i.test(l))).toBe(true);
    expect(labels.some((l) => /Copy as MOLBLOCK/i.test(l) || /Copy as Mol/i.test(l))).toBe(true);
  });

  const favCount = async (): Promise<number> => page.evaluate((keys) => {
    for (const k of keys) {
      const v = JSON.parse(localStorage.getItem(k) ?? '[]');
      if (Array.isArray(v) && v.length && typeof v[0] === 'string' && v[0].length) return v.length;
    }
    return 0;
  }, FAV_KEYS);

  await softStep('Step 4: Favorites > Add to Favorites registers the current molecule', async () => {
    // "Add to Favorites" persists this.getMolFile(), which stays empty until the sketcher holds a
    // molecule. The seed molecule loads into the Dart sketcher unreliably via this JS-API inplace
    // entry point, so first type a molecule (proven reliable in later steps) to make it "current",
    // then drive the menu item. Retry the (open menu -> hover Favorites -> Add to Favorites) cycle
    // until the Favorites collection in localStorage actually grows.
    await dismissMenu(); // close the menu left open by Step 3 so openHamburger opens a fresh one
    await typeSmiles(CYCLOHEXANE);
    await expect.poll(async () => {
      await openHamburger();
      await waitForLabel(/^Favorites/);
      await page.evaluate(() => {
        const fav = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find((l) => /^Favorites/i.test((l.textContent || '').trim()));
        (fav?.closest('.d4-menu-item') as HTMLElement | undefined)?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      });
      await clickLabel(/Add to Favorites/);
      return favCount();
    }, {timeout: 45_000, intervals: [1000, 2000, 3000]}).toBeGreaterThan(0);
  });

  await softStep('Step 5: Type C1CCCCC1 in SMILES input + Enter — molecule updates', async () => {
    await typeSmiles(CYCLOHEXANE);
  });

  await softStep('Step 6: Reopen menu — Favorites persists, Recent + Favorites groups render', async () => {
    await openHamburger();
    await waitForLabel(/Favorites/);
    const labels = await menuLabels();
    expect(labels.some((l) => /Favorites/i.test(l)), `Favorites group gone on reopen: ${JSON.stringify(labels.slice(0, 20))}`).toBe(true);
    expect(labels.some((l) => /Recent/i.test(l)), `Recent group gone on reopen: ${JSON.stringify(labels.slice(0, 20))}`).toBe(true);
    // Favorites survives the edit + reopen (persisted collection).
    expect(await favCount(), 'Favorites collection lost after edit + reopen').toBeGreaterThan(0);
    // Recent list *content* (the edited SMILES appearing) is only populated via the external-mode
    // cell-editor OK path, which the ui.dialog(inplace) JS-API entry point cannot drive — deferred.
    await dismissMenu();
  });

  await softStep('Step 8: Copy as SMILES places the current molecule SMILES on the clipboard', async () => {
    if (!clipboardOk) { console.warn('[SKIP] Step 8 — navigator.clipboard unavailable (insecure origin)'); return; }
    await openHamburger();
    await clickLabel(/Copy as SMILES/);
    await expect.poll(() => clipboard(), {timeout: 10_000}).toBe(CYCLOHEXANE);
  });

  await softStep('Step 8b: Paste the copied SMILES restores it into the molecular input', async () => {
    if (!clipboardOk) { console.warn('[SKIP] Step 8b — navigator.clipboard unavailable (insecure origin)'); return; }
    await typeSmiles(ETHANOL);
    await smilesInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Control+V');
    await smilesInput.press('Enter');
    await expect.poll(() => smilesInput.inputValue(), {timeout: 10_000}).toBe(CYCLOHEXANE);
  });

  await softStep('Step 9: Copy as MOLBLOCK places a V2000 molblock on the clipboard', async () => {
    if (!clipboardOk) { console.warn('[SKIP] Step 9 — navigator.clipboard unavailable (insecure origin)'); return; }
    await openHamburger();
    await clickLabel(/Copy as MOL/);
    await expect.poll(() => clipboard(), {timeout: 10_000}).toContain('V2000');
    const mol = await clipboard();
    expect(mol, `MOLBLOCK missing terminator: ${mol.slice(0, 60)}`).toContain('M  END');
  });

  await softStep('Step 9b: Paste the copied MOLBLOCK restores the pre-edit molecule', async () => {
    if (!clipboardOk) { console.warn('[SKIP] Step 9b — navigator.clipboard unavailable (insecure origin)'); return; }
    await typeSmiles(ETHANOL);
    await smilesInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Control+V');
    // MOLBLOCK paste is applied to the sketch canvas (input is not echoed) — verify via a Copy-as-SMILES read-back.
    await openHamburger();
    await clickLabel(/Copy as SMILES/);
    await expect.poll(() => clipboard(), {timeout: 10_000}).toBe(CYCLOHEXANE);
  });

  await softStep('Step 10: Close sketcher — no console errors fired', async () => {
    await page.evaluate(() => {
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) (cancel as HTMLElement).click();
      });
    });
    await page.locator('.d4-dialog').first().waitFor({state: 'detached', timeout: 10_000});
    const errs = await page.evaluate(() => ((window as any).__sk_errors ?? []) as string[]);
    const sketcherErrs = errs.filter((e) => /sketcher|chem\.sketcher|copy as/i.test(e));
    expect(sketcherErrs.length, `Console errors during sketcher walk: ${JSON.stringify(sketcherErrs.slice(0, 5))}`).toBe(0);
    const relErrors = pageErrors.filter((e) => /sketcher|chem/i.test(e));
    expect(relErrors.length, `Uncaught page errors during sketcher walk: ${JSON.stringify(relErrors.slice(0, 5))}`).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
