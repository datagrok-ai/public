/* ---
sub_features_covered: [bio.calculate.get-region, bio.calculate.get-region.top-menu, bio.detector, bio.transform.convert-notation, bio.transform.convert-notation.top-menu, bio.transform.split-to-monomers, bio.transform.to-atomic-level]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv', units: 'fasta', monomers: 39},
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv', units: 'helm', monomers: 17},
  {name: 'MSA', path: 'System:AppData/Bio/tests/filter_MSA.csv', units: 'separator', monomers: 17},
];
// Opens the Bio top-menu and drills to a leaf, polling for each menu node instead of fixed sleeps.
async function openBioMenuItem(page: Page, submenu: string, leaf: string): Promise<void> {
  await page.evaluate(async ({submenu, leaf}) => {
    const waitNode = async (sel: string): Promise<HTMLElement> => {
      for (let i = 0; i < 100; i++) {
        const el = document.querySelector(sel) as HTMLElement | null;
        if (el) return el;
        await new Promise((r) => setTimeout(r, 50));
      }
      throw new Error(`menu node not found: ${sel}`);
    };
    (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
    (await waitNode(`[name="div-Bio---${submenu}"]`))
      .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    (await waitNode(`[name="div-Bio---${submenu}---${leaf}"]`)).click();
  }, {submenu, leaf});
}
for (const ds of datasets) {
  test(`Bio Convert matrix on ${ds.name}`, async ({page}) => {
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
      const hasMacromolecule = () => Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
        .some((c: any) => c.semType === 'Macromolecule');
      for (let i = 0; i < 100 && !hasMacromolecule(); i++)
        await new Promise((r) => setTimeout(r, 100));
      if (hasMacromolecule())
        for (let i = 0; i < 60; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise((r) => setTimeout(r, 200));
        }
    }, ds.path);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await page.evaluate(async () => {
      const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
      for (let i = 0; i < 30; i++) {
        for (const fn of probes) {
          try { await (grok as any).functions.call(fn, {}); return; } catch { /* not ready yet */ }
        }
        await new Promise((r) => setTimeout(r, 100));
      }
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
      await openBioMenuItem(page, 'Calculate', 'Extract-Region...');
      // Canonical editor: the platform hosts the widget and names the dialog after the function's
      // friendly name 'Get Sequence Region' (Bio/package.ts:473) -> [name="dialog-Get-Sequence-Region"].
      await page.locator('[name="dialog-Get-Sequence-Region"]').waitFor({timeout: 60_000});
      await page.locator('[name="dialog-Get-Sequence-Region"] [name="button-OK"]').click();
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
        () => document.querySelectorAll('[name="dialog-Get-Sequence-Region"]').length === 0,
        null, {timeout: 15_000}).catch(() => {});
    });
    await softStep(`${ds.name}: Transform > Convert Sequence Notation adds a new Macromolecule column`, async () => {
      const beforeMacro: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Macromolecule').length;
      });
      const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await openBioMenuItem(page, 'Transform', 'Convert-Sequence-Notation...');
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
      // A new Macromolecule column was added by THIS step (not the earlier Extract-Region step),
      // in a different notation than the source. Exact target notation is env/default-dependent
      // (see convert.md unresolved_ambiguities: polytool-vs-transform-convert), so not pinned to a literal.
      expect(info.newMacroCount).toBeGreaterThan(beforeMacro);
      expect(info.lastUnits).not.toBeNull();
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
      await openBioMenuItem(page, 'Transform', 'To-Atomic-Level...');
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
      await openBioMenuItem(page, 'Transform', 'Split-to-Monomers...');
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
      // N per-position Monomer columns = alignment width of the source (convert.md matrix / convert-run.md):
      // FASTA -> 39, HELM -> 17, MSA -> 17. A broken split emitting a single column must fail here.
      expect(monCount - beforeMonCount).toBe(ds.monomers);
    });
    finishSpec();
  });
}
