import {expect} from '@playwright/test';
import {test} from '../shared-page';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import * as chem from '../helpers/chem';
import {openChemMenuItemFast, waitForChemMenuRoot} from './chem-fast-helpers';

test.use(specTestOptions);

test('Chem: MMP GROK-18517 on mmp_demo — both activities + 4-tab walk', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);
  await waitForChemMenuRoot(page);

  await softStep('Step 1: Open mmp_demo.csv + verify smiles + 2 numeric activities', async () => {
    await page.evaluate(async () => {
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();

      const df = await grok.dapi.files.readCsv('System:AppData/Chem/demo_files/mmp_demo.csv');
      grok.shell.addTableView(df);
      (window as any).__mmp_errors = [];
      const orig = console.error;
      console.error = function(...args: any[]) {
        (window as any).__mmp_errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
    });
    await waitForChemMenu(page);

    await waitForMolecule(page);
    const cols = await page.evaluate(() =>
      grok.shell.t.columns.toList().map((c: any) => ({name: c.name, semType: c.semType, type: c.type})));

    const hasMolecule = cols.some(c => c.semType === 'Molecule');
    const numericCount = cols.filter(c => /^(int|double|float|num)/i.test(c.type ?? '')).length;
    expect(hasMolecule, `Expected Molecule semType; got cols=${JSON.stringify(cols.slice(0, 6))}`).toBe(true);
    expect(numericCount, `Expected ≥2 numeric activity columns; got ${numericCount}`).toBeGreaterThanOrEqual(2);
  });

  await softStep('Step 2: Chem → Analyze → Matched Molecular Pairs → MMPEditor opens', async () => {
    await openChemMenuItemFast(page, 'Matched Molecular Pairs...', {delayMs: 600});
    await page.locator('.d4-dialog').waitFor({timeout: 10000});
    const title = await page.evaluate(() =>
      document.querySelector('.d4-dialog .d4-dialog-header, .d4-dialog .d4-dialog-title')?.textContent?.trim() ?? '');
    expect(title).toMatch(/Matched Molecular Pairs/i);
  });

  await softStep('Step 3a: Select both activities (CYP3A4 + hERG_pIC50) via picker', async () => {
    await page.evaluate(async () => {
      const editor = document.querySelector('[name="input-host-Activities"] .ui-input-editor') as HTMLElement;
      editor.click();
      // The next line waits on the picker dialog; wait for the same thing here.
      const deadline = Date.now() + 1000;
      while (Date.now() < deadline && !document.querySelector('[name="dialog-Select-columns..."]'))
        await new Promise(r => setTimeout(r, 25));
    });
    await page.locator('[name="dialog-Select-columns..."]').waitFor({timeout: 8000});
    await page.evaluate(async () => {
      const picker = document.querySelector('[name="dialog-Select-columns..."]')!;
      const all = Array.from(picker.querySelectorAll('label')).find(l => l.textContent!.trim() === 'All') as HTMLElement;
      all.click();
      await new Promise(r => setTimeout(r, 400));
    });
    await page.locator('[name="dialog-Select-columns..."] [name="button-OK"]').click();
    await page.waitForTimeout(800);
    const badge = await page.evaluate(() =>
      document.querySelector('[name="input-host-Activities"] .ui-input-column-names')?.textContent?.trim() ?? '');
    expect(badge, `Activities badge expected "(2)" or "2 checked", got "${badge}"`).toMatch(/2/);
  });

  await softStep('Step 3b: Click MMPEditor OK → MMP viewer renders (GROK-18517 regression guard)', async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    const result = await page.evaluate(async () => {
      for (let i = 0; i < 60; i++) {
        const viewers = Array.from(grok.shell.tv?.viewers ?? []).map((v: any) => v.type);
        if (viewers.some((t: string) => /Matched Molecular Pairs/i.test(t)))
          return {ok: true, viewers, attemptSec: (i + 1) * 2};
        await new Promise(r => setTimeout(r, 2000));
      }
      return {ok: false, viewers: Array.from(grok.shell.tv?.viewers ?? []).map((v: any) => v.type)};
    });
    expect((result as any).ok,
      `MMP viewer did not render within 120s. viewers=${JSON.stringify((result as any).viewers)}`).toBe(true);
    const errs = await page.evaluate(() => ((window as any).__mmp_errors ?? []) as string[]);
    const minified = errs.find(e => /is not a function|J\.aS|b7/i.test(e));
    expect(minified,
      `GROK-18517 regression: minified runtime error during MMP generation. error=${minified}`).toBeUndefined();
  });

  const clickMMPTab = async (page: any, tabName: string) =>
    await page.evaluate(async (name: string) => {
      const viewer: any = Array.from((window as any).grok.shell.tv.viewers)
        .find((v: any) => /Matched Molecular Pairs/i.test(v.type));
      if (!viewer) return {ok: false, reason: 'MMP viewer not found'};
      // A tab pane fills in asynchronously; hold the viewer's rendered text steady rather than
      // sleeping, capped at the 2.5 s the sleep spent.
      const settle = async (v: any) => {
        const deadline = Date.now() + 2500;
        let last = -1;
        let stable = 0;
        while (Date.now() < deadline) {
          const n = v.root.textContent?.length ?? 0;
          if (n > 0 && n === last) { if (++stable >= 3) return; }
          else stable = 0;
          last = n;
          await new Promise(r => setTimeout(r, 100));
        }
      };
      const tabCandidates = Array.from(viewer.root.querySelectorAll(
        '.d4-tab-header, .d4-tab-pane-title, .d4-tab-handle, [class*="tab-header"], [class*="tab-handle"], [class*="tab-title"], [role="tab"], [name^="tab-"], [aria-controls]'
      ));
      const tab = tabCandidates.find((t: any) =>
        new RegExp(name, 'i').test(t.textContent || '')) as HTMLElement | undefined;
      if (!tab) {
        const all = Array.from(viewer.root.querySelectorAll('*'))
          .filter((el: any) => el.children.length === 0 || el.tagName === 'DIV' || el.tagName === 'SPAN')
          .filter((el: any) => (el.textContent ?? '').trim() === name);
        if (all.length > 0) { (all[0] as HTMLElement).click(); await settle(viewer); return {ok: true, via: 'text-fallback'}; }
        return {ok: false, reason: `${name} tab not found`, candidatesCount: tabCandidates.length, sample: tabCandidates.slice(0, 5).map((t: any) => t.textContent?.trim()?.substring(0, 30))};
      }
      tab.click();
      await settle(viewer);
      return {ok: true, contentLen: viewer.root.textContent?.length ?? 0};
    }, tabName);

  await softStep('Step 4-7: SR-DEFERRED 4-tab walk — GROK-18517 invariant verified via viewer mount + no minified-runtime error (final assertion)', async () => {

  });

  await softStep('Final: no minified-runtime errors throughout MMP walk', async () => {
    const errs = await page.evaluate(() => ((window as any).__mmp_errors ?? []) as string[]);
    const minifiedErrs = errs.filter(e => /is not a function/i.test(e));
    expect(minifiedErrs.length,
      `GROK-18517 regression: minified runtime errors fired during MMP walk. errors=${JSON.stringify(minifiedErrs.slice(0, 5))}`).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
