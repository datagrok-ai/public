import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('AddNewColumn: hints — hover shows function signature', async ({page}) => {
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

  await softStep('Open Add New Column dialog', async () => {
    await page.evaluate(() => (document.querySelector('[name="icon-add-new-column"]') as HTMLElement).click());
    await page.locator('.d4-dialog [name="button-Add-New-Column---OK"]').waitFor({timeout: 10000});
  });

  await softStep('Insert a function (Abs) into the formula field', async () => {
    await page.locator('.d4-dialog .cm-content').click();
    await page.keyboard.type('a');
    await page.waitForSelector('.cm-tooltip-autocomplete', {timeout: 3000});
    await page.keyboard.press('Enter');
    await page.waitForTimeout(400);
    const content = await page.evaluate(() => document.querySelector('.d4-dialog .cm-content')?.textContent);
    expect(content).toBe('Abs(num)');
  });

  await softStep('Hover function name → tooltip with signature appears', async () => {
    const sigText = await page.evaluate(async () => {
      const line = document.querySelector('.d4-dialog .cm-content .cm-line') as HTMLElement;
      const range = document.createRange();
      range.selectNodeContents(line);
      const r = range.getClientRects()[0];
      const x = r.x + 10, y = r.y + r.height/2;
      for (const type of ['mouseover', 'mouseenter', 'mousemove'])
        line.dispatchEvent(new MouseEvent(type, {bubbles: true, cancelable: true, clientX: x, clientY: y, view: window}));
      await new Promise(res => setTimeout(res, 1500));
      const tip = document.querySelector('.cm-tooltip-hover');
      return tip?.textContent ?? null;
    });
    expect(sigText).not.toBeNull();
    expect(sigText).toContain('Abs');
    expect(sigText).toContain('num');
  });

  await page.locator('[name="button-Add-New-Column---CANCEL"]').click();

  if (stepErrors.length > 0)
    throw new Error('Step errors:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
