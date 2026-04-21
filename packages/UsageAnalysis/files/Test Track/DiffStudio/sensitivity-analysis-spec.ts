import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
});

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e.message ?? String(e)}); }
}

test('DiffStudio Sensitivity Analysis scenario', async ({page}) => {
  test.setTimeout(300_000);

  const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
  const login = process.env.DATAGROK_LOGIN ?? 'admin';
  const password = process.env.DATAGROK_PASSWORD ?? 'admin';

  // Login
  await page.goto(baseUrl);
  await page.waitForTimeout(2000);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  // Baseline setup: selenium class, Tabs mode, close everything.
  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
  });

  // Setup: open Bioreactor via DiffStudio:demoBioreactor
  await softStep('Setup — Open Bioreactor (DiffStudio:demoBioreactor)', async () => {
    await page.evaluate(async () => {
      const func = DG.Func.find({package: 'DiffStudio', name: 'demoBioreactor'})[0];
      await func.apply({});
    });
    await page.waitForFunction(() => grok.shell.v?.name === 'Bioreactor', null, {timeout: 30000});
    const viewName = await page.evaluate(() => grok.shell.v?.name);
    expect(viewName).toBe('Bioreactor');
  });

  // Step 1: Edit toggle off; form with inputs opens
  await softStep('Step 1 — Form with inputs is visible on Bioreactor', async () => {
    const info = await page.evaluate(() => {
      const inputs = document.querySelectorAll('[name^="input-host-"]').length;
      return {inputs, viewName: grok.shell.v?.name};
    });
    // 2b observed 23 inputs; Edit-toggle state is not verifiable via a stable selector.
    expect(info.viewName).toBe('Bioreactor');
    expect(info.inputs).toBeGreaterThan(0);
  });

  // Step 2: Click Sensitivity icon → "Bioreactor - comparison" view opens
  await softStep('Step 2 — Click Sensitivity icon, comparison view opens', async () => {
    await page.locator('.diff-studio-ribbon-sa-icon').first().click();
    await page.waitForFunction(
      () => (grok.shell.v?.name ?? '').includes('comparison'),
      null, {timeout: 30000});
    const viewName = await page.evaluate(() => grok.shell.v?.name);
    expect(viewName).toContain('comparison');
  });

  // Step 3: Select Parameters (Process mode, FFox, KKox, FFred)
  // The parameter switchers live in Dart-side widgets and were not scripted during 2b —
  // instead verify that the SA form exposes input hosts for the required params.
  await softStep('Step 3 — Parameter input hosts (FFox, KKox, FFred) present', async () => {
    const hosts = await page.evaluate(() => {
      const want = ['input-host-FFox', 'input-host-KKox', 'input-host-FFred'];
      return want.map(n => ({name: n, found: !!document.querySelector(`[name="${n}"]`)}));
    });
    for (const h of hosts)
      expect(h.found, `${h.name} should exist on SA form`).toBe(true);
  });

  // Step 4: Click Run icon → 4 viewers open
  // 2b observed: Run click did not start analysis (no view, no balloon, no error).
  // Try page.locator with force-click; poll tv.viewers for 90s. Fail AMBIGUOUS if nothing happens —
  // this may require a parameter switcher to be enabled first (2b could not toggle the switchers).
  await softStep('Step 4 — Click Run; expect 4 viewers (AMBIGUOUS — switchers may be required)', async () => {
    await page.locator('.grok-icon.fal.fa-play').first().click({force: true});

    let viewerCount = 0;
    for (let i = 0; i < 45; i++) {
      await page.waitForTimeout(2000);
      viewerCount = await page.evaluate(() => {
        const tv = grok.shell.tv;
        if (!tv) return 0;
        return Array.from(tv.viewers).length;
      });
      if (viewerCount >= 4) break;
    }

    if (viewerCount < 4) {
      throw new Error(
        `AMBIGUOUS: Run click on .grok-icon.fal.fa-play did not produce 4 viewers ` +
        `(got ${viewerCount}) after 90s. The SA Run may silently no-op when no parameter ` +
        `switcher is enabled — enabling at least one FFox/KKox/FFred switch is likely required ` +
        `before the play handler produces results. No balloon/warning was emitted to indicate this.`);
    }
    expect(viewerCount).toBeGreaterThanOrEqual(4);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
