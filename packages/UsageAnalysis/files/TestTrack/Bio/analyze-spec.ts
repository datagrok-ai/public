import {test, expect} from '@playwright/test';
import {specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasets = [
  {name: 'FASTA', path: 'System:AppData/Bio/samples/FASTA.csv'},
  {name: 'HELM', path: 'System:AppData/Bio/samples/HELM.csv'},
  {name: 'MSA', path: 'System:AppData/Bio/samples/MSA.csv'},
];

for (const ds of datasets) {
  test(`Bio Analyze on ${ds.name}`, async ({page}) => {
    stepErrors.length = 0;

    // Phase 1: Navigate
    await page.goto(baseUrl);
    await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell, {timeout: 15000});

    // Phase 2: Open dataset
    await page.evaluate(async (path) => {
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
        .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
      if (hasBioChem) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise(r => setTimeout(r, 200));
        }
        await new Promise(r => setTimeout(r, 5000));
      }
    }, ds.path);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

    // Step 1: Bio > Analyze > Sequence Space
    await softStep('Sequence Space with defaults', async () => {
      await page.evaluate(async () => {
        document.querySelector('[name="div-Bio"]').click();
        await new Promise(r => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]').dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise(r => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze---Sequence-Space..."]').click();
      });
      await page.locator('[name="button-OK"]').waitFor({timeout: 10000});
      await page.locator('[name="button-OK"]').click();
      await page.waitForFunction(() => {
        const viewers = Array.from(grok.shell.tv.viewers);
        return viewers.some(v => v.type === 'Scatter plot');
      }, {}, {timeout: 90000});
    });

    // Close non-grid viewers
    await page.evaluate(() => {
      for (const v of Array.from(grok.shell.tv.viewers)) { if (v.type !== 'Grid') v.close(); }
    });

    // Step 2: Bio > Analyze > Activity Cliffs
    await softStep('Activity Cliffs with defaults', async () => {
      await page.evaluate(async () => {
        document.querySelector('[name="div-Bio"]').click();
        await new Promise(r => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]').dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise(r => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze---Activity-Cliffs..."]').click();
      });
      await page.locator('[name="button-OK"]').waitFor({timeout: 10000});
      await page.locator('[name="button-OK"]').click();
      await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.length > 4, {}, {timeout: 90000});
    });

    // Close non-grid viewers
    await page.evaluate(() => {
      for (const v of Array.from(grok.shell.tv.viewers)) { if (v.type !== 'Grid') v.close(); }
    });

    // Step 3: Bio > Analyze > Composition
    await softStep('Composition analysis', async () => {
      await page.evaluate(async () => {
        document.querySelector('[name="div-Bio"]').click();
        await new Promise(r => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]').dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise(r => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze---Composition"]').click();
      });
      await page.waitForFunction(() => {
        const viewers = Array.from(grok.shell.tv.viewers);
        return viewers.some(v => v.type === 'WebLogo');
      }, {}, {timeout: 30000});
    });

    // Step 4: Check Composition viewer properties
    await softStep('Composition viewer properties', async () => {
      await page.evaluate(async () => {
        const weblogo = document.querySelector('[name="viewer-WebLogo"]');
        weblogo.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2, clientX: 500, clientY: 500}));
        await new Promise(r => setTimeout(r, 500));
        const items = document.querySelectorAll('.d4-menu-item-label');
        for (const item of items) {
          if (item.textContent.trim() === 'Properties...') {
            item.closest('.d4-menu-item').click();
            break;
          }
        }
      });
      await page.waitForTimeout(500);
    });

    // Step 5: Re-run Sequence Space with changed params
    await softStep('Sequence Space with t-SNE + Levenshtein', async () => {
      await page.evaluate(() => {
        for (const v of Array.from(grok.shell.tv.viewers)) { if (v.type !== 'Grid') v.close(); }
      });
      await page.evaluate(async () => {
        document.querySelector('[name="div-Bio"]').click();
        await new Promise(r => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze"]').dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise(r => setTimeout(r, 300));
        document.querySelector('[name="div-Bio---Analyze---Sequence-Space..."]').click();
      });
      await page.locator('select:near(:text("Method"))').first().selectOption('t-SNE');
      await page.locator('select:near(:text("Similarity"))').first().selectOption('Levenshtein');
      await page.locator('[name="button-OK"]').click();
      await page.waitForFunction(() => {
        const viewers = Array.from(grok.shell.tv.viewers);
        return viewers.some(v => v.type === 'Scatter plot');
      }, {}, {timeout: 90000});
    });

    if (stepErrors.length > 0) {
      const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
    }
  });
}
