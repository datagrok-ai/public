/* ---
sub_features_covered: [bio.analyze.sequence-space, bio.analyze.sequence-space.editor, bio.analyze.sequence-space.top-menu, bio.analyze.sequence-space.transform]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import * as bio from '../helpers/bio';
test.use(specTestOptions);
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/tests/filter_FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/tests/filter_HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/tests/filter_MSA.csv'},
];
for (const ds of datasets) {
  test(`Bio Sequence Space on ${ds.name}`, async ({page}) => {
    test.setTimeout(180_000);
    stepErrors.length = 0;
    await loginToDatagrok(page);
    const detected = await page.evaluate(async (path) => {
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = true;
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      const hasMacro = (): boolean => Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
        .some((c: any) => c.semType === 'Macromolecule');
      const deadline = Date.now() + 20_000;
      while (Date.now() < deadline) {
        if (hasMacro()) return true;
        await new Promise((r) => setTimeout(r, 200));
      }
      return hasMacro();
    }, ds.path);
    expect(detected, `${ds.name}: Macromolecule column detected`).toBe(true);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await page.evaluate(async () => {
      const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
      const deadline = Date.now() + 30_000;
      while (Date.now() < deadline) {
        for (const fn of probes) {
          try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
        }
        await new Promise((r) => setTimeout(r, 300));
      }
    });
    await page.waitForFunction(() => {
      const f = (grok as any).functions.find;
      try { return !!(f && (f('Bio:sequenceSpaceTopMenu') || f('Bio:sequenceSpace'))); } catch { return false; }
    }, undefined, {timeout: 30_000});
    await softStep(`${ds.name}: Open Bio > Analyze > Sequence Space (defaults)`, async () => {
      await bio.openBioAnalyze(page, 'div-Bio---Analyze---Sequence-Space...');
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
      const title = await page.locator('.d4-dialog .d4-dialog-title').textContent();
      expect(title?.trim()).toBe('Sequence Space');
    });
    await softStep(`${ds.name}: Run with default parameters — ScatterPlot + embeddings appended`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 90_000});
      const result = await page.evaluate(() => ({
        hasScatter: Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        scatterCount: Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length,
        cols: grok.shell.tv.dataFrame.columns.length,
        colNames: Array.from({length: grok.shell.tv.dataFrame.columns.length},
          (_, i) => grok.shell.tv.dataFrame.columns.byIndex(i).name),
      }));
      expect(result.hasScatter).toBe(true);
      expect(result.scatterCount).toBe(1);
      expect(result.cols).toBeGreaterThan(baseCols);
      expect(result.colNames).toContain('Embed_X_1');
      expect(result.colNames).toContain('Embed_Y_1');
      expect(result.colNames.some((n: string) => /^Cluster \(DBSCAN\)/.test(n))).toBe(true);
    });
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });
    await softStep(`${ds.name}: Re-open Bio > Analyze > Sequence Space`, async () => {
      await bio.openBioAnalyze(page, 'div-Bio---Analyze---Sequence-Space...');
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
      });
      await page.waitForFunction(() => {
        const dlg = document.querySelector('.d4-dialog')!;
        const m = (dlg.querySelector('[name="input-Method"]') as HTMLSelectElement)
          ?? (dlg.querySelector('[name="input-host-Method"] select') as HTMLSelectElement);
        const s = (dlg.querySelector('[name="input-Similarity"]') as HTMLSelectElement)
          ?? (dlg.querySelector('[name="input-host-Similarity"] select') as HTMLSelectElement);
        return m?.value === 't-SNE' && s?.value === 'Levenshtein';
      }, undefined, {timeout: 15_000});
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
    await softStep(`${ds.name}: Run with edited parameters — second embedding result docks`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 90_000});
      const result = await page.evaluate(() => ({
        hasScatter: Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        cols: grok.shell.tv.dataFrame.columns.length,
        colNames: Array.from({length: grok.shell.tv.dataFrame.columns.length},
          (_, i) => grok.shell.tv.dataFrame.columns.byIndex(i).name),
      }));
      expect(result.hasScatter).toBe(true);
      expect(result.cols).toBeGreaterThan(baseCols);
      const embedXCount = result.colNames.filter((n: string) => /^Embed_X_\d+$/.test(n)).length;
      const embedYCount = result.colNames.filter((n: string) => /^Embed_Y_\d+$/.test(n)).length;
      const clusterCount = result.colNames.filter((n: string) => /^Cluster \(DBSCAN\)/.test(n)).length;
      expect(embedXCount).toBeGreaterThanOrEqual(2);
      expect(embedYCount).toBeGreaterThanOrEqual(2);
      expect(clusterCount).toBeGreaterThanOrEqual(2);
    });
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });
    finishSpec();
  });
}
