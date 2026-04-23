import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('AddNewColumn: autocomplete (demog)', async ({page}) => {
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

  await softStep('Open demog, open Add New Column', async () => {
    await page.evaluate(() => (document.querySelector('[name="icon-add-new-column"]') as HTMLElement).click());
    await page.locator('.d4-dialog [name="button-Add-New-Column---OK"]').waitFor({timeout: 10000});
  });

  await softStep('Type "a" → autocomplete tooltip appears', async () => {
    await page.locator('.d4-dialog .cm-content').click();
    await page.keyboard.type('a');
    await page.waitForSelector('.cm-tooltip-autocomplete', {timeout: 3000});
    const items = await page.$$eval('.cm-tooltip-autocomplete li', nodes => nodes.slice(0, 8).map(n => (n as HTMLElement).textContent));
    expect(items.length).toBeGreaterThan(0);
    expect(items).toContain('Abs');
  });

  await softStep('Select function via Enter → inserted as "Abs(num)"', async () => {
    await page.keyboard.press('Enter');
    await page.waitForTimeout(400);
    const content = await page.evaluate(() => document.querySelector('.d4-dialog .cm-content')?.textContent);
    expect(content).toBe('Abs(num)');
  });

  await softStep('Clear, type "a", select via ArrowDown+Enter → inserted as "Acos(num)"', async () => {
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.keyboard.type('a');
    await page.waitForSelector('.cm-tooltip-autocomplete', {timeout: 3000});
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(300);
    const content = await page.evaluate(() => document.querySelector('.d4-dialog .cm-content')?.textContent);
    expect(content).toBe('Acos(num)');
  });

  await softStep('Clear field; Ctrl+Space → autocomplete tooltip', async () => {
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.keyboard.press('Control+Space');
    await page.waitForSelector('.cm-tooltip-autocomplete', {timeout: 3000});
    const count = await page.$$eval('.cm-tooltip-autocomplete li', nodes => nodes.length);
    expect(count).toBeGreaterThan(0);
  });

  await softStep('Press $ → tooltip lists columns', async () => {
    await page.keyboard.press('Escape');
    await page.evaluate(() => (document.querySelector('.d4-dialog .cm-content') as HTMLElement).focus());
    await page.keyboard.type('$');
    await page.waitForSelector('.cm-tooltip-autocomplete', {timeout: 3000});
    const items = await page.$$eval('.cm-tooltip-autocomplete li', nodes => nodes.map(n => (n as HTMLElement).textContent));
    expect(items).toContain('HEIGHT');
    expect(items).toContain('WEIGHT');
    expect(items).toContain('AGE');
  });

  await page.locator('[name="button-Add-New-Column---CANCEL"]').click();

  if (stepErrors.length > 0)
    throw new Error('Step errors:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
