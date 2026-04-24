import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/samples/FASTA.csv', seqCol: 'Sequence'},
  {name: 'HELM', path: 'System:AppData/Bio/samples/HELM.csv', seqCol: 'HELM'},
  {name: 'MSA', path: 'System:AppData/Bio/samples/MSA.csv', seqCol: 'MSA'},
];

for (const ds of datasets) {
  test(`Bio Convert on ${ds.name}`, async ({page}) => {
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

    await softStep(`${ds.name}: Dataset loaded with Macromolecule sequence column`, async () => {
      const hasMacro: boolean = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.some((c: any) => c.semType === 'Macromolecule');
      });
      expect(hasMacro).toBe(true);
    });

    // Bio > Calculate > Extract Region... (dialog title: "Get Region")
    await softStep(`${ds.name}: Calculate > Get Region adds a region column`, async () => {
      const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Calculate"]')!
          .dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector('[name="div-Bio---Calculate---Extract-Region..."]') as HTMLElement).click();
      });
      await page.locator('[name="dialog-Get-Region"]').waitFor({timeout: 15000});
      await page.locator('[name="button-OK"]').click();
      await page.waitForFunction(
        (b) => grok.shell.tv.dataFrame.columns.length > b, before, {timeout: 20000});
      const added: boolean = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const last = df.columns.byIndex(df.columns.length - 1);
        return last.semType === 'Macromolecule' && /:.*\(\d+-\d+\)/.test(last.name);
      });
      expect(added).toBe(true);
    });

    // Bio > PolyTool > Convert... — scenario intent (step: "Check a new column"):
    // opening the dialog and confirming it should add a new column to the dataframe.
    // The dialog itself carries name=`dialog-To-Atomic-Level` on the current build
    // (shared widget), so we assert the scenario's intent — a new column appears.
    await softStep(`${ds.name}: PolyTool > Convert adds a new column`, async () => {
      const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 500));
        document.querySelector('[name="div-Bio---PolyTool"]')!
          .dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 500));
        (document.querySelector('[name="div-Bio---PolyTool---Convert..."]') as HTMLElement).click();
      });
      await page.locator('.d4-dialog').first().waitFor({timeout: 15000});
      await page.waitForTimeout(400);
      await page.evaluate(() => {
        const dlg = document.querySelector('.d4-dialog');
        (dlg?.querySelector('[name="button-OK"]') as HTMLElement)?.click();
      });
      await page.waitForFunction(
        (b) => grok.shell.tv.dataFrame.columns.length > b, before, {timeout: 120000});
      const after: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      expect(after).toBeGreaterThan(before);
      // The PolyTool Convert dialog uses name=`dialog-To-Atomic-Level` (shared
      // widget). Wait for it to fully detach so the next softStep's
      // `dialog-To-Atomic-Level` locator doesn't match two dialogs at once.
      await page.waitForFunction(
        () => document.querySelectorAll('[name="dialog-To-Atomic-Level"]').length === 0,
        null, {timeout: 15000}).catch(() => {});
    });

    // Bio > Transform > To Atomic Level...
    await softStep(`${ds.name}: Transform > To Atomic Level adds a Molecule column`, async () => {
      const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Transform"]')!
          .dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector('[name="div-Bio---Transform---To-Atomic-Level..."]') as HTMLElement).click();
      });
      await page.locator('[name="dialog-To-Atomic-Level"]').waitFor({timeout: 15000});
      await page.locator('[name="dialog-To-Atomic-Level"] [name="button-OK"]').click();
      await page.waitForFunction(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.some((c: any) => c.semType === 'Molecule');
      }, null, {timeout: 120000});
      const info: {hasMol: boolean, units: string | null} = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const mol: any = cols.find((c: any) => c.semType === 'Molecule');
        return {hasMol: !!mol, units: mol?.meta?.units ?? null};
      });
      void before;
      expect(info.hasMol).toBe(true);
      expect(info.units).toBe('molblock');
    });

    // Bio > Transform > Split to Monomers...
    await softStep(`${ds.name}: Transform > Split to Monomers adds Monomer columns`, async () => {
      const before: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 400));
        document.querySelector('[name="div-Bio---Transform"]')!
          .dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
        (document.querySelector('[name="div-Bio---Transform---Split-to-Monomers..."]') as HTMLElement).click();
      });
      await page.locator('[name="dialog-Split-to-Monomers"]').waitFor({timeout: 15000});
      await page.locator('[name="dialog-Split-to-Monomers"] [name="button-OK"]').click();
      await page.waitForFunction(
        (b) => grok.shell.tv.dataFrame.columns.length > b, before, {timeout: 30000});
      const monCount: number = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.filter((c: any) => c.semType === 'Monomer').length;
      });
      expect(monCount).toBeGreaterThan(0);
    });

    if (stepErrors.length > 0) {
      const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
    }
  });
}
