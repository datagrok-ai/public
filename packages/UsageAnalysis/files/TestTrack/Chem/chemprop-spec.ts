import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Chem: ChemProp training', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await softStep('Open mol1K.sdf via grok.data.files.openTable', async () => {
    await page.evaluate(async () => {
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      const df = await (grok.data as any).files.openTable('System:AppData/Chem/mol1K.sdf');
      grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 5000);
      });
      await new Promise(r => setTimeout(r, 6000));
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
    const info = await page.evaluate(() => ({
      rows: grok.shell.t.rowCount,
      cols: grok.shell.t.columns.names(),
      molCol: grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule')?.name,
    }));
    expect(info.rows).toBeGreaterThan(100);
    expect(info.molCol).toBeTruthy();
  });

  await softStep('Open ML → Models → Train Model view', async () => {
    await page.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement;
      ml.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 700));
      const train = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Train Model...') as HTMLElement;
      if (train) (train.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 3500));
    });
    const viewName = await page.evaluate(() => grok.shell.v?.name);
    expect(viewName).toBe('Predictive model');
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
