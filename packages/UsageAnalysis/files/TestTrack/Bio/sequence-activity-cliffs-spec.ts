import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import * as bio from '../helpers/bio';
test.use(specTestOptions);
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/samples/FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/samples/HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/samples/MSA.csv'},
];
for (const ds of datasets) {
  test(`Bio Sequence Activity Cliffs on ${ds.name}`, async ({page}) => {
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
    await page.evaluate(async () => {
      const candidates = ['Bio:activityCliffsTopMenu', 'Bio:activityCliffs', 'Bio:macromoleculeActivityCliffs'];
      const findAny = (names: string[]): boolean => {
        for (const n of names) {
          try {
            if ((grok as any).functions.find && (grok as any).functions.find(n)) return true;
          } catch { /* try next */ }
        }
        return false;
      };
      const deadline = Date.now() + 15_000;
      while (Date.now() < deadline) {
        if (findAny(candidates)) return;
        await new Promise((r) => setTimeout(r, 300));
      }
      await new Promise((r) => setTimeout(r, 1500));
    });
    await page.waitForTimeout(2000);
    await softStep(`${ds.name}: Open Bio > Analyze > Activity Cliffs (defaults)`, async () => {
      await bio.openBioAnalyze(page, 'div-Bio---Analyze---Activity-Cliffs...');
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
      const title = await page.locator('.d4-dialog .d4-dialog-title').textContent();
      expect(title?.trim()).toBe('Activity Cliffs');
    });
    await softStep(`${ds.name}: Run with default parameters — ScatterPlot + embeddings appended`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 240_000});
      const result = await page.evaluate(() => ({
        hasScatter: Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        scatterCount: Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length,
        cols: grok.shell.tv.dataFrame.columns.length,
        hasParamsTag: !!grok.shell.tv.dataFrame.getTag('seqActivityCliffsParams'),
      }));
      expect(result.hasScatter).toBe(true);
      expect(result.scatterCount).toBe(1);
      expect(result.hasParamsTag).toBe(true);
    });
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });
    await softStep(`${ds.name}: Re-open Bio > Analyze > Activity Cliffs`, async () => {
      await bio.openBioAnalyze(page, 'div-Bio---Analyze---Activity-Cliffs...');
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
    });
    await softStep(`${ds.name}: Change Method to t-SNE and Similarity to Levenshtein`, async () => {
      await page.evaluate(async () => {
        const dlg = document.querySelector('.d4-dialog')!;
        const methodSel = (dlg.querySelector('[name="input-Method"]') as HTMLSelectElement)
          ?? (dlg.querySelector('[name="input-host-Method"] select') as HTMLSelectElement);
        const simSel = (dlg.querySelector('[name="input-Similarity"]') as HTMLSelectElement)
          ?? (dlg.querySelector('[name="input-host-Similarity"] select') as HTMLSelectElement);
        methodSel.value = 't-SNE';
        methodSel.dispatchEvent(new Event('input', {bubbles: true}));
        methodSel.dispatchEvent(new Event('change', {bubbles: true}));
        simSel.value = 'Levenshtein';
        simSel.dispatchEvent(new Event('input', {bubbles: true}));
        simSel.dispatchEvent(new Event('change', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
      });
      const verified = await page.evaluate(() => {
        const dlg = document.querySelector('.d4-dialog')!;
        const methodSel = (dlg.querySelector('[name="input-Method"]') as HTMLSelectElement)
          ?? (dlg.querySelector('[name="input-host-Method"] select') as HTMLSelectElement);
        const simSel = (dlg.querySelector('[name="input-Similarity"]') as HTMLSelectElement)
          ?? (dlg.querySelector('[name="input-host-Similarity"] select') as HTMLSelectElement);
        return {method: methodSel.value, sim: simSel.value};
      });
      expect(verified.method).toBe('t-SNE');
      expect(verified.sim).toBe('Levenshtein');
    });
    await softStep(`${ds.name}: Run with edited parameters — second cliff result docks`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 240_000});
      const result = await page.evaluate(() => ({
        hasScatter: Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        cols: grok.shell.tv.dataFrame.columns.length,
      }));
      expect(result.hasScatter).toBe(true);
      expect(result.cols).toBeGreaterThan(baseCols);
    });
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });
    finishSpec();
  });
}
