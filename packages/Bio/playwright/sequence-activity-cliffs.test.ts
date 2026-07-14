/* ---
sub_features_covered: [bio.analyze.activity-cliffs, bio.analyze.activity-cliffs.editor, bio.analyze.activity-cliffs.init, bio.analyze.activity-cliffs.top-menu, bio.analyze.activity-cliffs.transform]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import * as bio from '@datagrok-libraries/test/src/playwright/bio';
test.use(specTestOptions);
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/samples/FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/samples/HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/samples/MSA.csv'},
];
for (const ds of datasets) {
  test(`Bio Sequence Activity Cliffs on ${ds.name}`, async ({page}) => {
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
    }, ds.path);
    await page.waitForFunction(() =>
      !!grok.shell.t && grok.shell.t.columns.toList().some((c: any) => c.semType === 'Macromolecule'),
      null, {timeout: 30_000});
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
    await page.waitForFunction((names) => {
      const g = grok as any;
      return names.some((n: string) => { try { return !!(g.functions.find && g.functions.find(n)); } catch { return false; } });
    }, ['Bio:activityCliffsTopMenu', 'Bio:activityCliffs', 'Bio:macromoleculeActivityCliffs'], {timeout: 30_000});
    await softStep(`${ds.name}: Open Bio > Analyze > Activity Cliffs (defaults)`, async () => {
      await bio.openBioAnalyze(page, 'div-Bio---Analyze---Activity-Cliffs...');
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 60_000});
      const title = await page.locator('.d4-dialog .d4-dialog-title').textContent();
      // The Bio editor is registered as 'Sequence Activity Cliffs' (packages/Bio/src/package.ts:527);
      // tolerate the 'Sequence ' prefix rather than pinning the exact string.
      expect(title?.trim()).toContain('Activity Cliffs');
    });
    await softStep(`${ds.name}: Run with default parameters — ScatterPlot + embeddings appended`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base &&
          Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 120_000});
      const result = await page.evaluate(() => ({
        scatterCount: Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length,
        scatterTitles: Array.from((grok.shell.tv as any).viewers)
          .filter((v: any) => v.type === 'Scatter plot').map((v: any) => v.props?.title ?? ''),
        cols: grok.shell.tv.dataFrame.columns.length,
        paramsTag: grok.shell.tv.dataFrame.getTag('seqActivityCliffsParams'),
      }));
      expect(result.scatterCount).toBe(1);
      expect(result.scatterTitles).toContain('Activity cliffs');
      expect(result.cols).toBeGreaterThanOrEqual(baseCols + 2);
      expect(result.paramsTag, 'seqActivityCliffsParams tag missing').toBeTruthy();
      const params = JSON.parse(result.paramsTag as string);
      expect(Object.keys(params).length).toBeGreaterThan(0);
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
      const before = await page.evaluate(() => ({
        cols: grok.shell.tv.dataFrame.columns.length,
        scatterCount: Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length,
      }));
      expect(before.scatterCount).toBe(1);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      // GROK-19150: assert only that a distinct second viewer DOCKS (count=2); the
      // multi-instance click-routing isolation invariant is delegated to bio-grok-19150-spec.ts.
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base.cols &&
          Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot').length === 2,
        before, {timeout: 120_000});
      const result = await page.evaluate(() => {
        const scatters = Array.from((grok.shell.tv as any).viewers).filter((v: any) => v.type === 'Scatter plot');
        return {
          scatterCount: scatters.length,
          distinct: scatters.length === 2 && scatters[0] !== scatters[1],
          cols: grok.shell.tv.dataFrame.columns.length,
        };
      });
      expect(result.scatterCount).toBe(2);
      expect(result.distinct).toBe(true);
      expect(result.cols).toBeGreaterThan(before.cols);
    });
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers))
        if ((v as any).type !== 'Grid') (v as any).close();
    });
    finishSpec();
  });
}
