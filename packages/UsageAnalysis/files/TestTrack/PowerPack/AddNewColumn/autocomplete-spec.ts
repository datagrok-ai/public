import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e?.message ?? String(e)}); }
}

test('AddNewColumn: autocomplete (demog)', async ({page}) => {
  test.setTimeout(300_000);

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

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
