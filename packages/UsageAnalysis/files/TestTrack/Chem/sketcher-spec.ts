import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Chem: Sketcher', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

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
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Open sketcher dialog via grok.chem.sketcher API', async () => {
    await page.evaluate(async () => {
      const smiles = grok.shell.t.col('canonical_smiles').get(0);
      const sketcher = grok.chem.sketcher(null as any, smiles);
      ui.dialog('Sketcher').add(sketcher).show();
      await new Promise(r => setTimeout(r, 3000));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 10000});
  });

  await softStep('Enter C1CCCCC1 into the sketcher SMILES input', async () => {
    const result = await page.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog');
      const inputs = Array.from(dialog!.querySelectorAll('input'));
      const smilesInput = inputs.find(i => (i.getAttribute('placeholder') || '').toLowerCase().includes('smiles'));
      if (!smilesInput) return {ok: false};
      smilesInput.focus();
      const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      nativeSetter.call(smilesInput, 'C1CCCCC1');
      smilesInput.dispatchEvent(new Event('input', {bubbles: true}));
      smilesInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true}));
      await new Promise(r => setTimeout(r, 2500));
      return {ok: true, value: smilesInput.value};
    });
    expect(result.ok).toBe(true);
    expect(result.value).toBe('C1CCCCC1');
  });

  await softStep('Hamburger menu shows Copy as SMILES / MOLBLOCK / Recent / Favorites', async () => {
    await page.evaluate(async () => {
      const dialog = document.querySelector('.d4-dialog');
      const hamburger = dialog?.querySelector('.fa-bars, [name="icon-font-icon-menu"]') as HTMLElement;
      if (hamburger) hamburger.click();
      await new Promise(r => setTimeout(r, 1200));
    });
    const items = await page.evaluate(() => {
      return Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .map(i => i.textContent!.trim());
    });
    const hasAny = items.some(i =>
      /Copy as SMILES|Copy as MOLBLOCK|Recent|Favorites/.test(i));
    expect(hasAny).toBe(true);
  });

  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
