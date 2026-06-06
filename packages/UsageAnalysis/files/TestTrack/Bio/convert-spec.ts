import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv', units: 'fasta'},
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv', units: 'helm'},
  {name: 'MSA', path: 'System:AppData/Bio/tests/filter_MSA.csv', units: 'separator'},
];
for (const ds of datasets) {
  test(`Bio Convert matrix on ${ds.name}`, async ({page}) => {
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
    await softStep(`${ds.name}: Dataset opens with Macromolecule sequence column (units=${ds.units})`, async () => {
      const info: {hasMacro: boolean, units: string | null} = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        return {hasMacro: !!macro, units: macro?.meta?.units ?? null};
      });
      expect(info.hasMacro).toBe(true);
      expect(info.units).toBe(ds.units);
    });
    await softStep(`${ds.name}: Calculate > Extract Region adds a Macromolecule sub-region column`, async () => {
      const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Calculate"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector('[name="div-Bio---Calculate---Extract-Region..."]') as HTMLElement).click();
      });
      // Dialog name is `dialog-Get-Region` (named after API getRegion, not the menu label).
      await page.locator('[name="dialog-Get-Region"]').waitFor({timeout: 60_000});
      await page.locator('[name="dialog-Get-Region"] [name="button-OK"]').click();
      await page.waitForFunction(
        (b) => grok.shell.tv.dataFrame.columns.length > b, before, {timeout: 30_000});
      const info: {hasRegion: boolean, units: string | null, name: string | null} = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const region: any = cols.reverse().find((c: any) =>
          c.semType === 'Macromolecule' && /:\s*\(\d+-\d+\)/.test(c.name));
        return {hasRegion: !!region, units: region?.meta?.units ?? null, name: region?.name ?? null};
      });
      expect(info.hasRegion).toBe(true);
      expect(info.units).toBe(ds.units);
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-Get-Region"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });
    await softStep(`${ds.name}: Transform > Convert Sequence Notation adds a new Macromolecule column`, async () => {
      const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Transform"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector(
          '[name="div-Bio---Transform---Convert-Sequence-Notation..."]') as HTMLElement).click();
      });
      await page.locator('[name="dialog-Convert-Sequence-Notation"]').waitFor({timeout: 60_000});
      await page.locator('[name="dialog-Convert-Sequence-Notation"] [name="button-OK"]').click();
      await page.waitForFunction(
        (b) => grok.shell.tv.dataFrame.columns.length > b, before, {timeout: 60_000});
      const info: {newMacroCount: number, lastUnits: string | null, allUnits: string[]} =
        await page.evaluate(() => {
          const df = grok.shell.tv.dataFrame;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const macroCols = cols.filter((c: any) => c.semType === 'Macromolecule');
          const last: any = macroCols[macroCols.length - 1];
          return {
            newMacroCount: macroCols.length,
            lastUnits: last?.meta?.units ?? null,
            allUnits: macroCols.map((c: any) => c.meta?.units ?? null).filter((u: any) => u !== null),
          };
        });
      expect(info.newMacroCount).toBeGreaterThanOrEqual(2);
      expect(info.lastUnits).not.toBe(ds.units);
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-Convert-Sequence-Notation"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });
    await softStep(`${ds.name}: Transform > To Atomic Level adds a Molecule molblock column`, async () => {
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
      const info: {hasMol: boolean, units: string | null} = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const mol: any = cols.reverse().find((c: any) => c.semType === 'Molecule');
        return {hasMol: !!mol, units: mol?.meta?.units ?? null};
      });
      expect(info.hasMol).toBe(true);
      expect(info.units).toBe('molblock');
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-To-Atomic-Level"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });
    await softStep(`${ds.name}: Transform > Split to Monomers adds per-position Monomer columns`, async () => {
      const beforeMonCount: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Monomer').length;
      });
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Transform"]')!
          .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector(
          '[name="div-Bio---Transform---Split-to-Monomers..."]') as HTMLElement).click();
      });
      await page.locator('[name="dialog-Split-to-Monomers"]').waitFor({timeout: 60_000});
      await page.locator('[name="dialog-Split-to-Monomers"] [name="button-OK"]').click();
      await page.waitForFunction((base) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Monomer').length > base;
      }, beforeMonCount, {timeout: 60_000});
      const monCount: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Monomer').length;
      });
      expect(monCount).toBeGreaterThan(beforeMonCount);
      expect(monCount).toBeGreaterThan(0);
    });
    finishSpec();
  });
}
