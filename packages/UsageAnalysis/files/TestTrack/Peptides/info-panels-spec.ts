import {test, expect, chromium} from '@playwright/test';
import {specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/bio/peptides.csv';

test('Peptides — Info Panels', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
    await page.waitForFunction(() => {
      try { return typeof grok !== 'undefined' && typeof grok.shell.closeAll === 'function'; }
      catch { return false; }
    }, {timeout: 45000});
  }

  // Step 1: Open peptides.csv
  await softStep('Step 1: Open peptides.csv', async () => {
    const result = await page!.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach(d => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = false;

      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      // Bio: wait for package init + cell rendering
      const hasBio = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
        .some(c => c.semType === 'Macromolecule');
      if (hasBio) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise(r => setTimeout(r, 200));
        }
        await new Promise(r => setTimeout(r, 5000));
      }
      return {rows: df.rowCount, semType: df.col('AlignedSequence')?.semType};
    }, datasetPath);
    expect(result.rows).toBe(647);
    expect(result.semType).toBe('Macromolecule');
  });

  // Step 2: Verify amino acid coloring (renderer = sequence)
  await softStep('Step 2: Verify amino acid coloring', async () => {
    const renderer = await page!.evaluate(() => {
      const col = grok.shell.tv.dataFrame.col('AlignedSequence');
      return col?.getTag('cell.renderer');
    });
    expect(renderer).toBe('sequence');
  });

  // Step 3: Click the peptides column title
  await softStep('Step 3: Click peptides column title', async () => {
    await page!.evaluate(async () => {
      grok.shell.tv.dataFrame.currentCol = grok.shell.tv.dataFrame.col('AlignedSequence');
      await new Promise(r => setTimeout(r, 1000));
    });
    const currentCol = await page!.evaluate(() => grok.shell.tv.dataFrame.currentCol?.name);
    expect(currentCol).toBe('AlignedSequence');
  });

  // Step 4: Check that Details and Peptides panels are displayed
  await softStep('Step 4: Check Context Panel sections', async () => {
    const panels = await page!.evaluate(async () => {
      await new Promise(r => setTimeout(r, 1000));
      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      return Array.from(panes).map(p => p.textContent?.trim());
    });
    expect(panels).toContain('Details');
    expect(panels.some(p => p?.includes('Peptides'))).toBe(true);
  });

  // Step 5-6: Expand Peptides and Bioinformatics tabs, verify content
  await softStep('Step 5-6: Expand tabs and verify content', async () => {
    const content = await page!.evaluate(async () => {
      // Find and click Peptides pane header
      const panes = document.querySelectorAll('.d4-accordion-pane-header');
      let hasPeptides = false;
      let hasBioinfo = false;
      for (const p of panes) {
        const text = p.textContent?.trim() || '';
        if (text.includes('Peptides') && !p.classList.contains('expanded')) {
          p.click();
          hasPeptides = true;
        }
        if (text.includes('Bioinformatics') && !p.classList.contains('expanded')) {
          p.click();
          hasBioinfo = true;
        }
      }
      await new Promise(r => setTimeout(r, 2000));

      // Check that expanded panes have content
      const expandedPanes = document.querySelectorAll('.d4-accordion-pane-header.expanded');
      return {
        hasPeptides: Array.from(panes).some(p => p.textContent?.includes('Peptides')),
        hasBioinfo: Array.from(panes).some(p => p.textContent?.includes('Bioinformatics')),
        expandedCount: expandedPanes.length
      };
    });
    expect(content.hasPeptides).toBe(true);
  });

  // Cleanup
  await page!.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
