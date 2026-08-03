import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const datasets = [
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv', units: 'helm'},
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv', units: 'fasta'},
];
// GROK-15176: scan V3000 atom lines for heavy atoms illegally flagged MASS=1.
function heavyAtomIsotopeFlags(mol: string): string[] {
  if (!mol) return [];
  return mol.split(/\r?\n/)
    .filter((l) => /^M  V30 \d+ [A-Za-z]+ /.test(l))
    .filter((l) => / MASS=1\b/.test(l))
    .filter((l) => !/ [HD] /.test(l));
}
for (const ds of datasets) {
  test(`Bio Transform To Atomic Level + round-trip on ${ds.name}`, async ({page}) => {
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
    }, ds.path);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await page.evaluate(async () => {
      const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
      }
      await new Promise((r) => setTimeout(r, 3000));
    });
    await softStep(`${ds.name}: Macromolecule column has units=${ds.units}`, async () => {
      const info: {hasMacro: boolean, units: string | null} = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        return {hasMacro: !!macro, units: macro?.meta?.units ?? null};
      });
      expect(info.hasMacro).toBe(true);
      expect(info.units).toBe(ds.units);
    });
    await softStep(`${ds.name}: Bio | Transform | To Atomic Level adds a Molecule molblock column`, async () => {
      const beforeMolCount: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Molecule').length;
      });
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Transform"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector(
          '[name="div-Bio---Transform---To-Atomic-Level..."]') as HTMLElement).click();
      });
      await page.locator('[name="dialog-To-Atomic-Level"]').waitFor({timeout: 60_000});
      await page.locator('[name="dialog-To-Atomic-Level"] [name="button-OK"]').click();
      await page.waitForFunction((base) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Molecule').length > base;
      }, beforeMolCount, {timeout: 120_000});
      const info: {hasMol: boolean, units: string | null, name: string | null, firstCell: string | null} =
        await page.evaluate(() => {
          const df = grok.shell.tv.dataFrame;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const mol: any = cols.reverse().find((c: any) => c.semType === 'Molecule');
          let firstCell: string | null = null;
          if (mol) {
            for (let r = 0; r < Math.min(df.rowCount, 50); r++) {
              const v = mol.get(r);
              if (v != null && String(v).length > 0) { firstCell = String(v); break; }
            }
          }
          return {hasMol: !!mol, units: mol?.meta?.units ?? null, name: mol?.name ?? null, firstCell};
        });
      expect(info.hasMol).toBe(true);
      expect(info.units).toBe('molblock');
      expect(info.firstCell).not.toBeNull();
      expect(info.firstCell!).toMatch(/M\s+V30\s+BEGIN\s+CTAB/);
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-To-Atomic-Level"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });
    // Scenario 2 — Column Context Panel "To Atomic Level" action opens the same dialog.
    await softStep(`${ds.name}: Context Panel column-action "To Atomic Level" opens the same dialog`, async () => {
      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        if (macro) grok.shell.o = macro;
      });
      await page.waitForTimeout(800);
      await page.evaluate(() => {
        for (const h of Array.from(document.querySelectorAll('.d4-accordion-pane-header'))) {
          if ((h as HTMLElement).textContent === 'Actions')
            h.dispatchEvent(new MouseEvent('click', {bubbles: true}));
        }
      });
      await page.waitForTimeout(800);
      const clicked: boolean = await page.evaluate(() => {
        const labels = Array.from(document.querySelectorAll('label.d4-link-action'));
        const link = labels.find((l) => (l as HTMLElement).textContent?.trim()
          .toLowerCase().startsWith('to atomic level'));
        if (link) {
          (link as HTMLElement).click();
          return true;
        }
        return false;
      });
      expect(clicked).toBe(true);
      await page.locator('[name="dialog-To-Atomic-Level"]').waitFor({timeout: 30_000});
      const cancelBtn = page.locator('[name="dialog-To-Atomic-Level"] [name="button-CANCEL"]');
      if (await cancelBtn.count() > 0)
        await cancelBtn.first().click();
      else
        await page.keyboard.press('Escape');
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-To-Atomic-Level"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });
    await softStep(`${ds.name}: Bio | Transform | Molecules to HELM round-trip adds a Macromolecule column`, async () => {
      const beforeMacroCount: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Macromolecule').length;
      });
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Transform"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector(
          '[name="div-Bio---Transform---Molecules-to-HELM..."]') as HTMLElement).click();
      });
      await page.locator('[name="dialog-Molecules-to-HELM"]').waitFor({timeout: 60_000});
      await page.locator('[name="dialog-Molecules-to-HELM"] [name="button-OK"]').click();
      await page.waitForFunction((base) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Macromolecule').length > base;
      }, beforeMacroCount, {timeout: 180_000});
      const info: {newMacroCount: number, lastUnits: string | null, lastName: string | null} =
        await page.evaluate(() => {
          const df = grok.shell.tv.dataFrame;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const macros = cols.filter((c: any) => c.semType === 'Macromolecule');
          const last: any = macros[macros.length - 1];
          return {
            newMacroCount: macros.length,
            lastUnits: last?.meta?.units ?? null,
            lastName: last?.name ?? null,
          };
        });
      expect(info.newMacroCount).toBeGreaterThan(beforeMacroCount);
      expect(info.lastUnits).toBe('helm');
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-Molecules-to-HELM"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });
    // Scenario 5 — GROK-15176: to-atomic-level must not flag heavy atoms with MASS=1 (breaks PubChem standardization).
    await softStep(`${ds.name}: API wrappers seq2atomic + toAtomicLevelSingleSeq produce V3K molfiles (GROK-15176 isotope guard)`, async () => {
      const out: {molLinear: string | null, molNonLinear: string | null, linearErr: string | null, nonLinearErr: string | null} =
        await page.evaluate(async () => {
          const result: any = {molLinear: null, molNonLinear: null, linearErr: null, nonLinearErr: null};
          try {
            const r = await (grok as any).functions.call('Bio:toAtomicLevelSingleSeq', {sequence: 'ACDEFGHIK'});
            result.molLinear = typeof r === 'string' ? r : (r?.molfile ?? r?.result ?? null);
          } catch (e: any) {
            result.linearErr = String(e?.message ?? e);
          }
          try {
            const r = await (grok as any).functions.call('Bio:seq2atomic', {
              seq: 'PEPTIDE1{A.C.D.E.F.G.H.I.K}$$$$V2.0', nonlinear: true});
            result.molNonLinear = typeof r === 'string' ? r : (r?.molfile ?? r?.result ?? null);
          } catch (e: any) {
            result.nonLinearErr = String(e?.message ?? e);
          }
          return result;
        });
      expect(out.linearErr).toBeNull();
      expect(out.nonLinearErr).toBeNull();
      expect(out.molLinear).not.toBeNull();
      expect(out.molNonLinear).not.toBeNull();
      expect(out.molLinear!.length).toBeGreaterThan(0);
      expect(out.molNonLinear!.length).toBeGreaterThan(0);
      expect(out.molLinear!).toMatch(/V3000/);
      expect(out.molNonLinear!).toMatch(/V3000/);
      expect(out.molLinear!).toMatch(/M\s+V30\s+BEGIN\s+CTAB/);
      expect(out.molNonLinear!).toMatch(/M\s+V30\s+BEGIN\s+CTAB/);
      // GROK-15176 invariant: no heavy atom carries MASS=1.
      const offendersLinear = heavyAtomIsotopeFlags(out.molLinear!);
      const offendersNonLinear = heavyAtomIsotopeFlags(out.molNonLinear!);
      expect(offendersLinear, `GROK-15176 regression on linear wrapper: ${offendersLinear.length} heavy-atom MASS=1 flag(s). Offending atom lines:\n${offendersLinear.join('\n')}`).toEqual([]);
      expect(offendersNonLinear, `GROK-15176 regression on non-linear wrapper: ${offendersNonLinear.length} heavy-atom MASS=1 flag(s). Offending atom lines:\n${offendersNonLinear.join('\n')}`).toEqual([]);
    });
    finishSpec();
  });
}
