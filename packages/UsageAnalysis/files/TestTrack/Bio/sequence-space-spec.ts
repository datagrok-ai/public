import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/samples/FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/samples/HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/samples/MSA.csv'},
];

for (const ds of datasets) {
  test(`Bio Sequence Space on ${ds.name}`, async ({page}) => {
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
      const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
      if (hasBioChem) {
        for (let i = 0; i < 60; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        await new Promise((r) => setTimeout(r, 5000));
      }
    }, ds.path);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

    // Step 1 (open dataset) handled in setup above.

    await softStep(`${ds.name}: Open Bio > Analyze > Sequence Space (defaults)`, async () => {
      // Scenario text says "Bio > Search > Sequence Space" but the actual menu
      // path on the live platform is "Bio > Analyze > Sequence Space...".
      await page.locator('[name="div-Bio"]').hover();
      await page.locator('[name="div-Bio---Analyze"]').hover();
      await page.locator('[name="div-Bio---Analyze---Sequence-Space..."]').click();
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 15000});
      const title = await page.locator('.d4-dialog .d4-dialog-title').textContent();
      expect(title?.trim()).toBe('Sequence Space');
    });

    await softStep(`${ds.name}: Run with default parameters`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base
          && Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 240_000});
      const hasScatter = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'));
      expect(hasScatter).toBe(true);
    });

    // Close the scatter plot from the first run so we can clearly observe the second one.
    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers)) {
        if ((v as any).type !== 'Grid') (v as any).close();
      }
    });

    await softStep(`${ds.name}: Re-open Bio > Analyze > Sequence Space`, async () => {
      await page.locator('[name="div-Bio"]').hover();
      await page.locator('[name="div-Bio---Analyze"]').hover();
      await page.locator('[name="div-Bio---Analyze---Sequence-Space..."]').click();
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({timeout: 15000});
    });

    await softStep(`${ds.name}: Change Method to t-SNE and Similarity to Levenshtein`, async () => {
      // During MCP run we verified the named hosts: input-host-Method (UMAP/t-SNE),
      // input-host-Similarity (Hamming/Levenshtein/Monomer chemical distance/Needlemann-Wunsch).
      await page.evaluate(async () => {
        const dlg = document.querySelector('.d4-dialog')!;
        const methodSel = dlg.querySelector('[name="input-host-Method"] select') as HTMLSelectElement;
        const simSel = dlg.querySelector('[name="input-host-Similarity"] select') as HTMLSelectElement;
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
        return {
          method: (dlg.querySelector('[name="input-host-Method"] select') as HTMLSelectElement).value,
          sim: (dlg.querySelector('[name="input-host-Similarity"] select') as HTMLSelectElement).value,
        };
      });
      expect(verified.method).toBe('t-SNE');
      expect(verified.sim).toBe('Levenshtein');
    });

    await softStep(`${ds.name}: Run with edited parameters`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base
          && Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 240_000});
      const hasScatter = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'));
      expect(hasScatter).toBe(true);
    });

    if (stepErrors.length > 0) {
      const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
    }
  });
}
