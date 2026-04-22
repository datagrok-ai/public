import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  storageState: process.env.DATAGROK_STORAGE_STATE,
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Chem: Similarity Search', async ({page}) => {
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

  await softStep('Chem → Search → Similarity Search → viewer appears', async () => {
    await page.evaluate(async () => {
      const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const sim = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(m => m.textContent!.trim() === 'Similarity Search...') as HTMLElement;
      (sim.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 6000));
    });
    const hasSim = await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).some((v: any) => /Similarity/i.test(v.type || '')));
    expect(hasSim).toBe(true);
  });

  await softStep('Modify viewer options (fingerprint/limit/metric/cutoff) without error', async () => {
    const results = await page.evaluate(async () => {
      const simViewer: any = Array.from(grok.shell.tv.viewers).find((v: any) => /Similarity/i.test(v.type || ''));
      if (!simViewer) return {error: 'viewer not found'};
      const res: Record<string, boolean> = {};
      const beforeErr = (grok.shell.warnings || []).length;
      try { simViewer.setOptions({fingerprint: 'Pattern'}); await new Promise(r => setTimeout(r, 1500)); res.fingerprint = true; } catch (e) { res.fingerprint = false; }
      try { simViewer.setOptions({limit: 5}); await new Promise(r => setTimeout(r, 1500)); res.limit = true; } catch (e) { res.limit = false; }
      try { simViewer.setOptions({distanceMetric: 'Dice'}); await new Promise(r => setTimeout(r, 1500)); res.metric = true; } catch (e) { res.metric = false; }
      try { simViewer.setOptions({cutoff: 1.0}); await new Promise(r => setTimeout(r, 1500)); res.cutoff = true; } catch (e) { res.cutoff = false; }
      simViewer.setOptions({fingerprint: 'Morgan', limit: 12, distanceMetric: 'Tanimoto', cutoff: 0.01});
      (res as any).errDelta = (grok.shell.warnings || []).length - beforeErr;
      return res;
    });
    expect(results.fingerprint).toBe(true);
    expect(results.limit).toBe(true);
    expect(results.metric).toBe(true);
    expect(results.cutoff).toBe(true);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
