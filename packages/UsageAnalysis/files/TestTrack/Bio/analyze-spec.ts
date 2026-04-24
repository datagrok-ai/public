import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/samples/FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/samples/HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/samples/MSA.csv'},
];

for (const ds of datasets) {
  test(`Bio Analyze on ${ds.name}`, async ({page}) => {
    test.setTimeout(600_000);
    stepErrors.length = 0;

    await loginToDatagrok(page);

    await page.evaluate(async (path) => {
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = true;
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(() => resolve(), 3000);
      });
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
      if (hasBioChem) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        await new Promise((r) => setTimeout(r, 5000));
      }
    }, ds.path);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

    await softStep(`${ds.name}: Sequence Space with defaults`, async () => {
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]')!.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Sequence-Space..."]') as HTMLElement).click();
      });
      await page.locator('[name="button-OK"]').waitFor({timeout: 15000});
      await page.locator('[name="button-OK"]').click();
      await page.waitForFunction(() => {
        const viewers = Array.from((grok.shell.tv as any).viewers);
        return viewers.some((v: any) => v.type === 'Scatter plot');
      }, null, {timeout: 180000});
      const hasScatter = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'));
      expect(hasScatter).toBe(true);
    });

    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers)) { if ((v as any).type !== 'Grid') (v as any).close(); }
    });

    await softStep(`${ds.name}: Activity Cliffs with defaults`, async () => {
      const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]')!.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Activity-Cliffs..."]') as HTMLElement).click();
      });
      await page.locator('[name="button-OK"]').waitFor({timeout: 15000});
      await page.locator('[name="button-OK"]').click();
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base
          && Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'Scatter plot'),
        baseCols, {timeout: 180000});
      const cols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      expect(cols).toBeGreaterThan(baseCols);
    });

    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers)) { if ((v as any).type !== 'Grid') (v as any).close(); }
    });

    await softStep(`${ds.name}: Composition`, async () => {
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]')!.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Composition"]') as HTMLElement).click();
      });
      await page.waitForFunction(() => {
        const viewers = Array.from((grok.shell.tv as any).viewers);
        return viewers.some((v: any) => v.type === 'WebLogo');
      }, null, {timeout: 60000});
      const hasWebLogo = await page.evaluate(() =>
        Array.from((grok.shell.tv as any).viewers).some((v: any) => v.type === 'WebLogo'));
      expect(hasWebLogo).toBe(true);
    });

    // Step 4 only applies to Composition viewer (WebLogo). Open Properties via right-click
    // because WebLogo viewer does not render title bar icons even with selenium class.
    if (ds.name === 'FASTA') {
      await softStep(`${ds.name}: Composition viewer properties on context panel`, async () => {
        await page.evaluate(async () => {
          const weblogo = document.querySelector('[name="viewer-WebLogo"]') as HTMLElement;
          const canvas = weblogo.querySelector('canvas') as HTMLCanvasElement;
          const rect = canvas.getBoundingClientRect();
          canvas.dispatchEvent(new MouseEvent('contextmenu', {
            bubbles: true, cancelable: true, button: 2,
            clientX: rect.left + 5, clientY: rect.top + 5,
          }));
          await new Promise((r) => setTimeout(r, 500));
          const labels = Array.from(document.querySelectorAll('.d4-menu-item-label'));
          for (const l of labels) {
            if (l.textContent?.trim() === 'Properties...') {
              (l.closest('.d4-menu-item') as HTMLElement).click();
              break;
            }
          }
          await new Promise((r) => setTimeout(r, 800));
        });
      });
    }

    await page.evaluate(() => {
      for (const v of Array.from((grok.shell.tv as any).viewers)) { if ((v as any).type !== 'Grid') (v as any).close(); }
    });

    await softStep(`${ds.name}: Sequence Space with t-SNE + Levenshtein`, async () => {
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]')!.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 300));
        (document.querySelector('[name="div-Bio---Analyze---Sequence-Space..."]') as HTMLElement).click();
      });
      await page.locator('[name="button-OK"]').waitFor({timeout: 15000});
      await page.evaluate(async () => {
        const dialog = document.querySelector('.d4-dialog')!;
        const selects = Array.from(dialog.querySelectorAll('select')) as HTMLSelectElement[];
        // selects[2] = Method, selects[3] = Similarity (verified during MCP run)
        selects[2].value = 't-SNE';
        selects[2].dispatchEvent(new Event('change', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 200));
        selects[3].value = 'Levenshtein';
        selects[3].dispatchEvent(new Event('change', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 200));
      });
      await page.locator('[name="button-OK"]').click();
      await page.waitForFunction(() => {
        const viewers = Array.from((grok.shell.tv as any).viewers);
        return viewers.some((v: any) => v.type === 'Scatter plot');
      }, null, {timeout: 180000});
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
