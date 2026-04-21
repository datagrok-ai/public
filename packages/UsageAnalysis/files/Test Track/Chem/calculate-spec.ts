import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  storageState: process.env.DATAGROK_STORAGE_STATE,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';
const datasetPath = 'System:DemoFiles/chem/smiles.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Chem: Calculate Descriptors', async ({page}) => {
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
  await page.locator('[name="Browse"]').waitFor({timeout: 60000});

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 500));
    const df = await grok.dapi.files.readCsv(path);
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
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Open Chem → Calculate → Descriptors dialog', async () => {
    await page.locator('[name="div-Chem"]').click();
    await page.waitForTimeout(500);
    await page.locator('.d4-menu-item-label', {hasText: /^Calculate$/}).first().hover();
    await page.waitForTimeout(500);
    await page.locator('.d4-menu-item-label', {hasText: /^Descriptors\.\.\.$/}).first().click();
    await page.locator('.d4-dialog').waitFor({timeout: 10000});
  });

  await softStep('Accept defaults → OK → new descriptor columns appended', async () => {
    const colsBefore = await page.evaluate(() => grok.shell.t.columns.length);
    // Expand a descriptor group first so at least one is selected, then OK.
    await page.evaluate(async () => {
      // Check any tree node (e.g., first descriptor category like "RDKit descriptors") to enable OK
      const tree = document.querySelector('.d4-dialog .d4-tree-view-root');
      if (tree) {
        const firstGroup = tree.querySelector('.d4-tree-view-checkbox');
        if (firstGroup) (firstGroup as HTMLElement).click();
      }
      await new Promise(r => setTimeout(r, 500));
    });
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(20000);
    const colsAfter = await page.evaluate(() => grok.shell.t.columns.length);
    expect(colsAfter).toBeGreaterThan(colsBefore);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
