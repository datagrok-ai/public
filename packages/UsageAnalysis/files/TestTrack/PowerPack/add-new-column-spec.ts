import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('PowerPack: Add new columns (Demog)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    // @ts-ignore
    grok.shell.settings.showFiltersIconsConstantly = true;
    // @ts-ignore
    grok.shell.windows.simpleMode = true;
    // @ts-ignore
    grok.shell.closeAll();
    // @ts-ignore
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    // @ts-ignore
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Open Demog dataset (HEIGHT/WEIGHT present)', async () => {
    const ok = await page.evaluate(() => {
      // @ts-ignore
      const df = grok.shell.tv.dataFrame;
      const names = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
      return names.includes('HEIGHT') && names.includes('WEIGHT');
    });
    expect(ok).toBe(true);
  });

  await softStep('Press the Add new column icon; dialog opens', async () => {
    await page.locator('[name="icon-add-new-column"]').click();
    await page.locator('.d4-dialog [name="button-Add-New-Column---OK"]').waitFor({timeout: 10000});
    await page.locator('.d4-dialog input.ui-input-addnewcolumn-name').waitFor({timeout: 10000});
    await page.locator('.d4-dialog .cm-content').waitFor({timeout: 10000});
  });

  await softStep('UI check: no overflow in contents, resize handles present', async () => {
    const res = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog') as HTMLElement;
      const contents = dlg.querySelector('.d4-dialog-contents') as HTMLElement;
      const overflowX = contents.scrollWidth > contents.clientWidth;
      const overflowY = contents.scrollHeight > contents.clientHeight;
      const resizers = dlg.querySelectorAll('[class*="resizer"]').length;
      return {overflowX, overflowY, resizers};
    });
    expect(res.overflowX).toBe(false);
    expect(res.overflowY).toBe(false);
    expect(res.resizers).toBeGreaterThan(0);
  });

  await softStep('Set name "New" and formula Round(${HEIGHT} + ${WEIGHT})', async () => {
    const nameInput = page.locator('.d4-dialog input.ui-input-addnewcolumn-name');
    await nameInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('New');
    await page.locator('.d4-dialog .cm-content').click();
    await page.keyboard.type('Round(${HEIGHT} + ${WEIGHT})');
    await page.waitForTimeout(500);
    const okEnabled = await page.evaluate(() =>
      document.querySelector('[name="button-Add-New-Column---OK"]')?.classList.contains('enabled'));
    expect(okEnabled).toBe(true);
  });

  await softStep('Press OK; "New" column is added', async () => {
    await page.locator('[name="button-Add-New-Column---OK"]').click();
    await page.waitForFunction(() => {
      // @ts-ignore
      const df = grok.shell.tv.dataFrame;
      const names = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
      return names.includes('New') && document.querySelectorAll('.d4-dialog').length === 0;
    }, {timeout: 15000});
  });

  await softStep('Reopen dialog, click history icon, select most recent activity (autofill)', async () => {
    await page.locator('[name="icon-add-new-column"]').click();
    await page.locator('.d4-dialog [name="icon-history"]').waitFor({timeout: 10000});
    await page.locator('.d4-dialog input.ui-input-addnewcolumn-name').waitFor({timeout: 10000});
    await page.locator('.d4-dialog [name="icon-history"]').click();
    const item = page.locator('.d4-menu-popup .d4-menu-item').first();
    await item.waitFor({timeout: 5000});
    await item.click();
    await page.waitForTimeout(1500);
    const filled = await page.evaluate(() => {
      const input = document.querySelector('.d4-dialog input.ui-input-addnewcolumn-name') as HTMLInputElement;
      const cm = document.querySelector('.d4-dialog .cm-content') as HTMLElement;
      return {name: input?.value, formula: cm?.textContent};
    });
    expect(filled.name).toBe('New');
    expect(filled.formula).toContain('Round');
    expect(filled.formula).toContain('HEIGHT');
    expect(filled.formula).toContain('WEIGHT');
    await page.locator('[name="button-Add-New-Column---CANCEL"]').click();
  });

  if (stepErrors.length > 0)
    throw new Error('Step errors:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
