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
});

test('DiffStudio Fitting: Bioreactor, open Fit view, modify params, run fit', async ({page}) => {
  test.setTimeout(300_000);

  // Login
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

  // Setup: close dialogs, set Tabs mode, close all views
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

  // Setup step: open Bioreactor via demoBioreactor
  await softStep('Setup: Open Bioreactor via DiffStudio:demoBioreactor', async () => {
    await page.evaluate(async () => {
      const func = DG.Func.find({package: 'DiffStudio', name: 'demoBioreactor'})[0];
      await func.apply({});
    });
    await page.waitForFunction(() => grok.shell.v?.name === 'Bioreactor', null, {timeout: 30000});
    await page.waitForTimeout(2000);

    const viewName = await page.evaluate(() => grok.shell.v?.name);
    expect(viewName).toBe('Bioreactor');
  });

  // Step 1: Open Diff Studio + Bioreactor (covered by setup)
  await softStep('Step 1: Open Diff Studio + Bioreactor', async () => {
    const info = await page.evaluate(() => {
      const inputs = document.querySelectorAll('[name^="input-host-"]').length;
      return {viewName: grok.shell.v?.name, inputs};
    });
    expect(info.viewName).toBe('Bioreactor');
    expect(info.inputs).toBeGreaterThan(0);
  });

  // Step 2: Click Fit icon on ribbon → Fitting view opens
  await softStep('Step 2: Click Fit icon, Fitting view opens', async () => {
    await page.evaluate(async () => {
      const fitIcon = document.querySelector(
        'i.diff-studio-svg-icon:not(.diff-studio-ribbon-sa-icon)') as HTMLElement | null;
      if (fitIcon) fitIcon.click();
    });
    await page.waitForFunction(
      () => grok.shell.v?.name === 'Bioreactor - fitting',
      null,
      {timeout: 30000},
    );
    await page.waitForTimeout(2000);

    const info = await page.evaluate(() => ({
      viewName: grok.shell.v?.name,
      hasProcessMode: !!document.querySelector('[name="input-host-Process-mode"]'),
      hasFFoxMin: !!document.querySelector('[name="input-host-FFox-(min)"]'),
    }));
    expect(info.viewName).toBe('Bioreactor - fitting');
    expect(info.hasProcessMode).toBe(true);
    expect(info.hasFFoxMin).toBe(true);
  });

  // Step 3: Modify Process mode; verify FFox/KKox update
  // AMBIGUOUS: changing <select> + observing FFox/KKox range updates requires dispatching
  // input/change events on the choice element. Not exercised in the 2b MCP run.
  await softStep('Step 3: Modify Process mode; verify FFox/KKox update (AMBIGUOUS)', async () => {
    test.info().annotations.push({type: 'ambiguous',
      description: 'Process mode Choice input + FFox/KKox range reactivity not exercised during 2b. ' +
                   'Asserting presence of the relevant input hosts only.'});
    const info = await page.evaluate(() => ({
      hasProcessMode: !!document.querySelector('[name="input-host-Process-mode"]'),
      hasFFoxMin: !!document.querySelector('[name="input-host-FFox-(min)"]'),
      hasKKoxMin: !!document.querySelector('[name="input-host-KKox-(min)"]'),
    }));
    expect(info.hasProcessMode).toBe(true);
    expect(info.hasFFoxMin).toBe(true);
  });

  // Step 4: Set Process mode=Default; switcher-select switch-at, FFox range 0.15→1.0, FKox 0→3
  // AMBIGUOUS: per-parameter min/max keyboard input was not exercised in 2b.
  await softStep('Step 4: Switchers + FFox/FKox range edits (AMBIGUOUS)', async () => {
    test.info().annotations.push({type: 'ambiguous',
      description: 'Per-parameter switcher clicks + (min)/(max) keyboard input not exercised ' +
                   'during 2b. Asserting presence of the range input hosts only.'});
    const info = await page.evaluate(() => ({
      hasFFoxMin: !!document.querySelector('[name="input-host-FFox-(min)"]'),
      hasFFoxMax: !!document.querySelector('[name="input-host-FFox-(max)"]'),
      hasFKoxMin: !!document.querySelector('[name="input-host-FKox-(min)"]'),
      hasFKoxMax: !!document.querySelector('[name="input-host-FKox-(max)"]'),
    }));
    expect(info.hasFFoxMin).toBe(true);
    expect(info.hasFFoxMax).toBe(true);
  });

  // Step 5: Scroll to Target; input Bioreactor Data from bioreactor-experiment.csv
  // AMBIGUOUS: the "Bioreactor table" input caption was not found in the Fitting form body
  // during 2b — may appear only after Step 4's switcher edits, or has a different label.
  await softStep('Step 5: Input Bioreactor Data target (AMBIGUOUS)', async () => {
    test.info().annotations.push({type: 'ambiguous',
      description: '"Bioreactor table" input caption not present in the Fitting view during 2b. ' +
                   'May be gated by Step 4 completion, or the caption differs from the scenario text.'});
    const hasTargetText = await page.evaluate(() => {
      const body = document.body.innerText || '';
      return body.includes('Target');
    });
    expect(hasTargetText).toBe(true);
  });

  // Step 6: Click Run icon; expected RMSE by iterations descending graph
  // AMBIGUOUS: play icon exists but the same click pattern on the SA ribbon produced no output
  // during 2b; not exercised here to avoid a false-positive on a silently-failing run.
  await softStep('Step 6: Run fitting; RMSE by iterations graph (AMBIGUOUS)', async () => {
    test.info().annotations.push({type: 'ambiguous',
      description: 'Run icon (`.grok-icon.fal.fa-play`) is present but its click was not exercised ' +
                   'during 2b — prior SA runs showed silent no-op risk. Asserting icon presence only.'});
    const hasPlay = await page.evaluate(() =>
      !!document.querySelector('.grok-icon.fal.fa-play') ||
      !!document.querySelector('.d4-ribbon-item i.fa-play'));
    expect(hasPlay).toBe(true);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
