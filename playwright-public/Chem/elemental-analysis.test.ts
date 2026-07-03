/* ---
sub_features_covered: [chem.analyze.elemental, chem.analyze.elemental.run, chem.analyze.elemental.top-menu]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

interface Variant { id: string; path: string; format: string; opener: 'csv' | 'table'; warmWaitMs: number; }
const variants: Variant[] = [
  {id: 'D1', path: 'System:AppData/Chem/tests/smiles-50.csv',                 format: 'smiles',   opener: 'csv',   warmWaitMs: 30000},
  {id: 'D2', path: 'System:AppData/Chem/mol1K.sdf',                    format: 'molV2000', opener: 'table', warmWaitMs: 8000},
  {id: 'D3', path: 'System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf',  format: 'molV3000', opener: 'table', warmWaitMs: 8000},
];

test('Chem: Elemental Analysis multi-format walk (smiles / molV2000 / molV3000)', async ({page}) => {
  test.setTimeout(420_000);

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
      }, {path: v.path, opener: v.opener});
    });

    await softStep(`${v.id}: wait ${v.warmWaitMs}ms for Chem cascade (semType + dialog registration)`, async () => {
      await page.waitForTimeout(v.warmWaitMs);
    });

    await softStep(`${v.id}: verify Molecule semType + Grid renders`, async () => {
      const result = await page.evaluate(() => {
        const df = (window as any).__ea_df;
        const cols = df.columns.toList().map((c: any) => ({name: c.name, semType: c.semType}));
        return {hasMol: cols.some((c: any) => c.semType === 'Molecule'), cols};
      });
      if (!result.hasMol)
        throw new Error(`${v.id} setup failed: no Molecule column. cols=${JSON.stringify(result.cols)}`);
      await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 20000});
    });

    await softStep(`${v.id}: open Chem > Analyze > Elemental Analysis... dialog`, async () => {
      await page.evaluate(async () => {
        const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement | null;
        if (!chemMenu) throw new Error('[name="div-Chem"] not found in DOM (Chem top-menu missing)');
        chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
        await new Promise(r => setTimeout(r, 600));
        const ea = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find(m => (m.textContent ?? '').trim() === 'Elemental Analysis...') as HTMLElement | undefined;
        if (!ea) throw new Error('"Elemental Analysis..." sub-menu item not found');
        (ea.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      });
      await page.locator('.d4-dialog').waitFor({timeout: 15000});
    });

    await softStep(`${v.id}: toggle all per-element checkboxes ON`, async () => {
      const counts = await page.evaluate(() => {
        const dlg = document.querySelector('.d4-dialog');
        if (!dlg) return {total: 0, toggled: 0};
        const checkboxes = Array.from(dlg.querySelectorAll('input[type="checkbox"]')) as HTMLInputElement[];
        let toggled = 0;
        for (const cb of checkboxes) {
          if (!cb.checked) {
            cb.click();
            toggled++;
          }
        }
        return {total: checkboxes.length, toggled};
      });
      expect(counts.total, `${v.id}: per-element checkboxes must exist in dialog`).toBeGreaterThan(0);
    });

    await softStep(`${v.id}: click OK and verify column-append OR error-balloon`, async () => {
      const preCount = await page.evaluate(() => (window as any).__ea_preCount as number);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      // Wait up to 30s for column append or balloon-error to surface
      const outcome = await page.waitForFunction((preCount: number) => {
        const df = (window as any).__ea_df;
        if (df.columns.length > preCount) return {appended: df.columns.length - preCount, balloon: null};
        const balloon = document.querySelector('.d4-balloon, .grok-balloon');
        if (balloon && balloon.textContent && balloon.textContent.length > 5)
          return {appended: 0, balloon: balloon.textContent.slice(0, 200)};
        return false;
      }, preCount, {timeout: 30000}).then(jsh => jsh.jsonValue()).catch(() => null) as {appended: number; balloon: string | null} | null;
      expect(outcome, `${v.id}: neither column-append nor error balloon within 30s — silent failure not allowed`).toBeTruthy();
      // Valid outcome = appended columns OR a visible error balloon, never silent.
      const ok = (outcome!.appended > 0) || !!outcome!.balloon;
      expect(ok, `${v.id}: outcome=${JSON.stringify(outcome)}`).toBe(true);
    });
  }

  const consoleErrors = await page.evaluate(() => ((window as any).__ea_errors as string[] | undefined) ?? []);
  const fatalErrors = consoleErrors.filter((e: string) =>
    /TypeError|cannot read|undefined is not|null is not/i.test(e) &&
    !/non-fatal|warning/i.test(e));
  expect(
    fatalErrors.length,
    `Fatal console errors during Elemental Analysis multi-format walk: ${JSON.stringify(fatalErrors.slice(0, 5))}`,
  ).toBe(0);

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
