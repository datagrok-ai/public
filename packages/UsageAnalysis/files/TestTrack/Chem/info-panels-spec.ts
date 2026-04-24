import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Chem: Info Panels on smiles.csv', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup: open smiles.csv, wait for chem
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 500));
    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    await new Promise(r => setTimeout(r, 5000));
    grok.shell.windows.showContextPanel = true;
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Step 1-2: Verify smiles.csv loaded', async () => {
    const info = await page.evaluate(() => ({
      rows: grok.shell.t.rowCount,
      molCol: grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule')?.name,
    }));
    expect(info.rows).toBe(1000);
    expect(info.molCol).toBe('canonical_smiles');
  });

  await softStep('Step 3: Select canonical_smiles column → context panel shows column panes', async () => {
    const panes = await page.evaluate(async () => {
      grok.shell.windows.showContextPanel = true;
      grok.shell.o = grok.shell.t.col('canonical_smiles');
      await new Promise(r => setTimeout(r, 2500));
      return Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .map(h => h.textContent!.trim());
    });
    expect(panes).toContain('Chemistry');
    expect(panes).toContain('Details');
  });

  await softStep('Step 4: Expand all column info panels → no errors', async () => {
    const result = await page.evaluate(async () => {
      const targets = ['Details', 'Filter', 'Colors', 'Style', 'Settings', 'Plots', 'Advanced', 'Sticky meta', 'Chemistry'];
      const warns = (grok.shell.warnings || []).length;
      const balloons = Array.from(document.querySelectorAll('.grok-balloon-error')).length;
      for (const name of targets) {
        const pane = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
          .filter(h => h.textContent!.trim() === name).pop();
        if (pane && !pane.classList.contains('expanded')) {
          (pane as HTMLElement).click();
          await new Promise(r => setTimeout(r, 600));
        }
      }
      await new Promise(r => setTimeout(r, 1000));
      return {
        warnDelta: (grok.shell.warnings || []).length - warns,
        balloonDelta: Array.from(document.querySelectorAll('.grok-balloon-error')).length - balloons,
      };
    });
    expect(result.warnDelta).toBe(0);
    expect(result.balloonDelta).toBe(0);
  });

  await softStep('Step 5: Click first molecule → molecule context with no errors', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.windows.showContextPanel = true;
      await new Promise(r => setTimeout(r, 500));
      const warns = (grok.shell.warnings || []).length;
      const balloons = Array.from(document.querySelectorAll('.grok-balloon-error')).length;
      // Trigger molecule selection via grid cell
      const grid = grok.shell.tv.grid;
      grid.dataFrame.currentCell = grid.dataFrame.cell(0, 'canonical_smiles');
      const sv = DG.SemanticValue.fromTableCell(grok.shell.t.cell(0, 'canonical_smiles'));
      grok.shell.o = sv;
      // Wait for panels to populate
      for (let i = 0; i < 20; i++) {
        await new Promise(r => setTimeout(r, 500));
        const hs = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
          .map(h => h.textContent!.trim());
        if (hs.includes('Structure') || hs.includes('Chemistry')) break;
      }
      const headers = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .map(h => h.textContent!.trim());
      const molPanes = ['Chemistry', 'Biology', 'Databases', 'Structure'];
      for (const name of molPanes) {
        const pane = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
          .filter(h => h.textContent!.trim() === name).pop();
        if (pane && !pane.classList.contains('expanded')) {
          (pane as HTMLElement).click();
          await new Promise(r => setTimeout(r, 800));
        }
      }
      await new Promise(r => setTimeout(r, 1500));
      return {
        headers,
        warnDelta: (grok.shell.warnings || []).length - warns,
        balloonDelta: Array.from(document.querySelectorAll('.grok-balloon-error')).length - balloons,
      };
    });
    expect(result.headers.some((h: string) => h === 'Chemistry' || h === 'Structure')).toBe(true);
    expect(result.warnDelta).toBe(0);
    expect(result.balloonDelta).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
