/* ---
sub_features_covered: [chem.analyze.elemental, chem.analyze.elemental.run, chem.analyze.elemental.top-menu]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

interface Variant { id: string; path: string; format: string; opener: 'csv' | 'table'; }
const variants: Variant[] = [
  {id: 'D1', path: 'System:AppData/Chem/tests/smiles-50.csv',                 format: 'smiles',   opener: 'csv'},
  {id: 'D2', path: 'System:AppData/Chem/mol1K.sdf',                    format: 'molV2000', opener: 'table'},
  {id: 'D3', path: 'System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf',  format: 'molV3000', opener: 'table'},
];

test('Chem: Elemental Analysis multi-format walk (smiles / molV2000 / molV3000)', async ({page}) => {
  test.setTimeout(180_000);

  await loginToDatagrok(page);

  await softStep('Setup: close all + selenium flags + hook console.error', async () => {
    await page.evaluate(() => {
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      (window as any).__ea_errors = [];
      const orig = console.error;
      console.error = function(...args: any[]) {
        (window as any).__ea_errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
    });
    await page.waitForTimeout(500);
  });

  for (const v of variants) {
    await softStep(`${v.id} (${v.format}): open ${v.path}`, async () => {
      await page.evaluate(async ({path}: {path: string; opener: string}) => {
        grok.shell.closeAll();
        const isSdf = path.toLowerCase().endsWith('.sdf');
        let df: any;
        if (isSdf) {
          await ((DG as any).Func.find({name: 'OpenFile'})[0])
            .prepare({fullPath: path}).call(undefined, undefined, {processed: false});
          for (let __i = 0; __i < 30 && !grok.shell.t; __i++) await new Promise(r => setTimeout(r, 200));
          df = grok.shell.t;
        } else {
          df = await grok.dapi.files.readCsv(path);
          grok.shell.addTableView(df);
        }
        (window as any).__ea_df = df;
        (window as any).__ea_preCount = df.columns.length;
        (window as any).__ea_preRowCount = df.rowCount;
      }, {path: v.path, opener: v.opener});
    });

    await softStep(`${v.id}: verify Molecule semType + Grid renders`, async () => {
      await waitForMolecule(page);
      await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 20000});
    });

    await softStep(`${v.id}: open Chem > Analyze > Elemental Analysis... dialog`, async () => {
      await page.evaluate(async () => {
        const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement | null;
        if (!chemMenu) throw new Error('[name="div-Chem"] not found in DOM (Chem top-menu missing)');
        chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
        const findEa = () => Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find(m => (m.textContent ?? '').trim() === 'Elemental Analysis...') as HTMLElement | undefined;
        let ea = findEa();
        for (let i = 0; i < 100 && !ea; i++) { await new Promise(r => setTimeout(r, 50)); ea = findEa(); }
        if (!ea) throw new Error('"Elemental Analysis..." sub-menu item not found');
        (ea.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      });
      await page.locator('.d4-dialog').waitFor({timeout: 15000});
    });

    await softStep(`${v.id}: toggle all per-element checkboxes ON`, async () => {
      const counts = await page.evaluate(() => {
        const dlg = document.querySelector('.d4-dialog');
        if (!dlg) return {total: 0, allChecked: false, checkedCount: 0};
        const checkboxes = Array.from(dlg.querySelectorAll('input[type="checkbox"]')) as HTMLInputElement[];
        for (const cb of checkboxes)
          if (!cb.checked) cb.click();
        const checkedCount = checkboxes.filter(cb => cb.checked).length;
        (window as any).__ea_checkedCount = checkedCount;
        return {total: checkboxes.length, allChecked: checkboxes.every(cb => cb.checked), checkedCount};
      });
      expect(counts.total, `${v.id}: per-element checkboxes must exist in dialog`).toBeGreaterThan(0);
      expect(counts.allChecked, `${v.id}: all per-element checkboxes must be ON (checked ${counts.checkedCount}/${counts.total})`).toBe(true);
    });

    await softStep(`${v.id}: click OK and verify per-element atom-count columns appended`, async () => {
      const preCount = await page.evaluate(() => (window as any).__ea_preCount as number);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      // Poll only for column append; read any balloon purely to enrich the failure message.
      const outcome = await page.waitForFunction((preCount: number) => {
        const df = (window as any).__ea_df;
        if (df.columns.length > preCount) return {appended: df.columns.length - preCount, balloon: null};
        const balloon = document.querySelector('.d4-balloon, .grok-balloon');
        if (balloon && balloon.textContent && balloon.textContent.length > 5)
          return {appended: 0, balloon: balloon.textContent.slice(0, 200)};
        return false;
      }, preCount, {timeout: 30000}).then(jsh => jsh.jsonValue()).catch(() => null) as {appended: number; balloon: string | null} | null;
      expect(outcome, `${v.id}: no per-element columns appended within 30s (balloon=${outcome && (outcome as any).balloon})`).toBeTruthy();
      expect(outcome!.appended, `${v.id}: Elemental Analysis must append per-element columns, not fail with balloon=${outcome!.balloon}`).toBeGreaterThan(0);

      const verify = await page.evaluate(() => {
        const df = (window as any).__ea_df;
        const preRow = (window as any).__ea_preRowCount as number;
        const preCol = (window as any).__ea_preCount as number;
        const appendedCols = df.columns.toList().slice(preCol);
        const first = appendedCols[0];
        const raw = first.getRawData ? first.getRawData() : null;
        const n = Math.min(df.rowCount, 20);
        let allInt = true;
        if (raw) for (let i = 0; i < n; i++) { const val = raw[i]; if (!Number.isInteger(val) || val < 0) { allInt = false; break; } }
        return {rowCount: df.rowCount, preRow, firstName: first.name, allInt: raw ? allInt : true};
      });
      expect(verify.rowCount, `${v.id}: grid row count must be unchanged after append`).toBe(verify.preRow);
      expect(verify.allInt, `${v.id}: appended column "${verify.firstName}" must hold non-negative integer atom counts`).toBe(true);
    });
  }

  const consoleErrors = await page.evaluate(() => ((window as any).__ea_errors as string[] | undefined) ?? []);
  // .md step 5: "No console errors fire" — gate on every console.error, allowlisting only explicit non-fatal/warning noise.
  const errors = consoleErrors.filter((e: string) => !/non-fatal|warning/i.test(e));
  expect(
    errors.length,
    `Console errors during Elemental Analysis multi-format walk: ${JSON.stringify(errors.slice(0, 5))}`,
  ).toBe(0);

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
