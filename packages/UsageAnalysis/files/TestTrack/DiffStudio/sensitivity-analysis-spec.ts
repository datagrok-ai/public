import {test, expect} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

test('DiffStudio Sensitivity Analysis: Bioreactor, Process mode cascade, Run, 4 viewers', async ({page}) => {
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
  await page.waitForTimeout(2000);

  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  });

  await softStep('Step 1: Open Diff Studio + Bioreactor; Edit toggle OFF', async () => {
    await page.evaluate(async () => {
      const f = DG.Func.find({name: 'runDiffStudio', package: 'DiffStudio'})[0];
      const call = f.prepare();
      await call.call();
      grok.shell.addView(call.getOutputParamValue());
    });
    await page.waitForFunction(() => grok.shell.v?.name === 'Diff Studio', null, {timeout: 30000});
    await page.waitForTimeout(2000);
    const card = page.locator('.diff-studio-hub-card', {hasText: 'Bioreactor'}).first();
    await card.waitFor({timeout: 15000});
    await card.dblclick();
    await page.waitForFunction(() => grok.shell.v?.name === 'Bioreactor', null, {timeout: 30000});
    await page.waitForTimeout(2000);
    // Make sure Edit is OFF so form with inputs is visible
    const editOn = await page.evaluate(() =>
      document.querySelector('.d4-ribbon-item .ui-input-switch')?.classList.contains('ui-input-switch-on'));
    if (editOn) {
      await page.evaluate(() => {
        const editItem = Array.from(document.querySelectorAll('.d4-ribbon-item'))
          .find(el => el.textContent?.trim() === 'Edit') as HTMLElement;
        (editItem?.querySelector('.ui-input-bool-switch .ui-input-editor') as HTMLElement)?.click();
      });
      await page.waitForTimeout(1500);
    }
    const hasProcessMode = await page.locator('[name="input-host-Process-mode"]').count();
    expect(hasProcessMode).toBeGreaterThan(0);
  });

  await softStep('Step 2: Open Sensitivity Analysis view', async () => {
    await page.evaluate(() => {
      const item = Array.from(document.querySelectorAll('.d4-ribbon-item'))
        .find(el => el.textContent?.trim() === 'Sensitivity') as HTMLElement;
      const icon = item?.querySelector('.diff-studio-ribbon-sa-icon') as HTMLElement;
      const rect = icon.getBoundingClientRect();
      icon.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: rect.left+5, clientY: rect.top+5}));
      icon.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: rect.left+5, clientY: rect.top+5}));
      icon.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: rect.left+5, clientY: rect.top+5}));
    });
    await page.waitForFunction(() =>
      grok.shell.v?.name?.includes('comparison') || document.body.innerText.includes('Sensitivity Analysis'),
      null, {timeout: 30000});
    await page.waitForTimeout(2000);
    const sa = await page.evaluate(() => ({
      viewName: grok.shell.v?.name,
      hasSAText: document.body.innerText.includes('Sensitivity Analysis'),
      hasProcessMode: !!document.querySelector('[name="input-host-Process-mode"]'),
    }));
    expect(sa.hasSAText).toBe(true);
    expect(sa.hasProcessMode).toBe(true);
  });

  await softStep('Step 3: Modify Process mode; FFox/KKox/others cascade-update', async () => {
    const before = await page.evaluate(() => ({
      ffox: (document.querySelector('[name="input-host-FFox"] input.ui-input-editor') as HTMLInputElement)?.value,
      kkox: (document.querySelector('[name="input-host-KKox"] input.ui-input-editor') as HTMLInputElement)?.value,
      mea: (document.querySelector('[name="input-host-MEAthiol"] input.ui-input-editor') as HTMLInputElement)?.value,
      gas: (document.querySelector('[name="input-host-Gas"] input.ui-input-editor') as HTMLInputElement)?.value,
    }));
    await page.evaluate(() => {
      const select = document.querySelector('[name="input-host-Process-mode"] select') as HTMLSelectElement;
      select.value = 'Mode 1';
      select.dispatchEvent(new Event('input', {bubbles: true}));
      select.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.waitForTimeout(2500);
    const after = await page.evaluate(() => ({
      ffox: (document.querySelector('[name="input-host-FFox"] input.ui-input-editor') as HTMLInputElement)?.value,
      kkox: (document.querySelector('[name="input-host-KKox"] input.ui-input-editor') as HTMLInputElement)?.value,
      mea: (document.querySelector('[name="input-host-MEAthiol"] input.ui-input-editor') as HTMLInputElement)?.value,
      gas: (document.querySelector('[name="input-host-Gas"] input.ui-input-editor') as HTMLInputElement)?.value,
    }));
    const changed = Object.keys(before).filter(k => (before as any)[k] !== (after as any)[k]);
    expect(changed.length).toBeGreaterThanOrEqual(3);
  });

  await softStep('Step 4: Run Sensitivity Analysis; 4 viewers open', async () => {
    // Enable 2 input switchers so SA has varied inputs to analyze
    await page.evaluate(() => {
      (document.querySelector('[name="input-host-FFox"] .ui-input-switch') as HTMLElement)?.click();
      (document.querySelector('[name="input-host-KKox"] .ui-input-switch') as HTMLElement)?.click();
    });
    await page.waitForTimeout(1500);
    // Click Run (▶) on ribbon — fa-play icon
    await page.evaluate(() => {
      const runIcon = Array.from(document.querySelectorAll('.d4-ribbon-item i'))
        .find(e => e.className.includes('fa-play')) as HTMLElement;
      runIcon?.click();
    });
    // Wait for viewers to materialize
    await page.waitForTimeout(15000);
    const viewerCount = await page.evaluate(() => document.querySelectorAll('.d4-viewer').length);
    expect(viewerCount).toBeGreaterThanOrEqual(4);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
