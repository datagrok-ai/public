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

test('DiffStudio Files & Sharing: pk.ivp preview, modify Step/Count, URL reload', async ({page, context}) => {
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

  await softStep('Step 1: Open pk.ivp file preview', async () => {
    await page.evaluate(async () => {
      const previewFn = DG.Func.find({name: 'previewIvp', package: 'DiffStudio'})[0];
      const info = await grok.dapi.files.list('System:AppData/DiffStudio/library/', false, 'pk.ivp');
      const file = info[0];
      const call = previewFn.prepare({file: file});
      await call.call();
      const view = call.getOutputParamValue();
      if (view) grok.shell.addView(view);
    });
    await page.waitForFunction(() => grok.shell.v?.name === 'PK', null, {timeout: 30000});
    await page.waitForTimeout(2000);
    const inputs = await page.locator('[name^="input-host-"]').count();
    expect(inputs).toBeGreaterThan(5);
  });

  await softStep('Step 2: Set Step=0.1 (slider+textbox) and Count=4 (clicker)', async () => {
    const stepInp = page.locator('[name="input-host-step"] input.ui-input-editor');
    await stepInp.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('0.1');
    await page.keyboard.press('Tab');
    await page.waitForTimeout(1500);

    // Count via clicker — force click (icons only visible on hover)
    for (let i = 0; i < 3; i++) {
      await page.locator('[name="input-host-count"] [name="icon-plus"]').click({force: true});
      await page.waitForTimeout(400);
    }
    await page.waitForTimeout(1500);

    const stepVal = await stepInp.inputValue();
    const countVal = await page.locator('[name="input-host-count"] input.ui-input-editor').inputValue();
    expect(parseFloat(stepVal)).toBeCloseTo(0.1);
    expect(countVal).toBe('4');
    const url = page.url();
    expect(url).toContain('step=0.1');
    expect(url).toContain('count=4');
  });

  await softStep('Step 3: Open URL in a new tab; model loads with same inputs', async () => {
    const url = page.url();
    const newTab = await context.newPage();
    await newTab.goto(url);
    // login if needed on new tab (shared context — session cookie carries)
    await newTab.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'PK', null, {timeout: 60000});
    await newTab.waitForTimeout(3000);
    const state = await newTab.evaluate(() => ({
      viewName: grok.shell.v?.name,
      stepVal: (document.querySelector('[name="input-host-step"] input.ui-input-editor') as HTMLInputElement)?.value,
      countVal: (document.querySelector('[name="input-host-count"] input.ui-input-editor') as HTMLInputElement)?.value,
    }));
    expect(state.viewName).toBe('PK');
    expect(parseFloat(state.stepVal)).toBeCloseTo(0.1);
    expect(state.countVal).toBe('4');
    await newTab.close();
  });

  // Step 4: REMARK — not a test but a consistency check. With count<=3 expect no Multiaxis/Facet.
  await softStep('Step 4: REMARK — with 1-3 curves no Multiaxis/Facet; with >=4 curves they appear', async () => {
    // Currently count=4 — check what's there
    const tabsAt4 = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.tab-handle-text'))
        .map(t => t.textContent?.trim()).filter(Boolean));
    // Then drop to 1 and re-check
    const countInp = page.locator('[name="input-host-count"] input.ui-input-editor');
    await countInp.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('1');
    await page.keyboard.press('Tab');
    await page.waitForTimeout(2500);
    const tabsAt1 = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.tab-handle-text'))
        .map(t => t.textContent?.trim()).filter(Boolean));
    test.info().annotations.push({type: 'remark',
      description: `tabsAt4=${tabsAt4.join(',')} tabsAt1=${tabsAt1.join(',')}`});
    // REMARK only — assertion is permissive
    expect(tabsAt1).not.toContain('Multiaxis');
    expect(tabsAt1).not.toContain('Facet');
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
