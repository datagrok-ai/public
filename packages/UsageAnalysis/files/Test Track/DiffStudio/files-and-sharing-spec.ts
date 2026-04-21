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

test('DiffStudio Files & Sharing: Open pk.ivp, modify inputs, share URL', async ({page}) => {
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

  // Step 1: Navigate to Browse > Files > App Data > Diff Studio > Library; click pk.ivp; preview opens
  await softStep('Step 1: Open pk.ivp preview via /files/ URL', async () => {
    await page.goto(`${baseUrl}/files/System.AppData/DiffStudio/library/pk.ivp`);
    await page.waitForTimeout(2000);

    // Wait for the file view to appear (name contains 'pk.ivp')
    await page.waitForFunction(() => {
      const v = grok.shell.v;
      return !!v && typeof v.name === 'string' && v.name.toLowerCase().includes('pk.ivp');
    }, null, {timeout: 30000});

    const info = await page.evaluate(() => {
      const v = grok.shell.v;
      const inputHosts = document.querySelectorAll('.ui-input-root').length;
      const bodyText = document.body.innerText || '';
      return {
        viewName: v?.name ?? null,
        viewType: v?.type ?? null,
        inputHosts,
        hasMultiaxis: bodyText.includes('Multiaxis'),
        hasFacet: bodyText.includes('Facet'),
      };
    });

    test.info().annotations.push({type: 'info', description:
      `view=${info.viewName} type=${info.viewType} .ui-input-root count=${info.inputHosts} ` +
      `multiaxis=${info.hasMultiaxis} facet=${info.hasFacet}`});

    // Assert a view did open and contains pk.ivp in its name. Do NOT assert interactive inputs —
    // on dev, the /files/ URL opens a non-interactive file preview (0 inputs observed in 2b).
    expect(info.viewName).toBeTruthy();
    expect(info.viewName!.toLowerCase()).toContain('pk.ivp');
  });

  // Step 2: Set Step to 0.1 via slider; Count to 4 via clicker; view updates
  // 2b run log: FAIL — input-host-Step and input-host-Count are NOT present in this preview DOM.
  await softStep('Step 2: Modify Step=0.1 and Count=4 (FAIL per 2b — inputs absent in file preview)', async () => {
    const result = await page.evaluate(() => {
      const stepHost = document.querySelector('[name="input-host-Step"]');
      const countHost = document.querySelector('[name="input-host-Count"]');
      const anyInputHosts = document.querySelectorAll('[name^="input-host-"]').length;
      const inputRoots = document.querySelectorAll('.ui-input-root').length;
      return {
        stepHostFound: !!stepHost,
        countHostFound: !!countHost,
        anyInputHosts,
        inputRoots,
      };
    });

    test.info().annotations.push({type: 'ambiguous', description:
      `input-host-Step found=${result.stepHostFound}, input-host-Count found=${result.countHostFound}, ` +
      `total input-host-* count=${result.anyInputHosts}, .ui-input-root count=${result.inputRoots}. ` +
      'The /files/ URL opens a non-interactive file preview — Step/Count inputs are not exposed there.'});

    // This step is FAIL per the 2b log — assert so the Playwright run records it as such.
    expect(result.stepHostFound,
      'input-host-Step not present: /files/ URL opens a raw file preview, not an interactive DiffStudio model').toBe(true);
    expect(result.countHostFound,
      'input-host-Count not present: /files/ URL opens a raw file preview, not an interactive DiffStudio model').toBe(true);
  });

  // Step 3: Copy URL; paste in new tab; same inputs appear
  // 2b run log: AMBIGUOUS — depends on step 2 producing a parameterized URL; step 2 blocked.
  await softStep('Step 3: Share URL (AMBIGUOUS per 2b — blocked by step 2)', async () => {
    test.info().annotations.push({type: 'ambiguous', description:
      'Step 3 cannot be exercised because step 2 is blocked on dev — the /files/ URL preview ' +
      'does not expose interactive Step/Count inputs, so there is no parameterized URL to copy ' +
      'and validate in a new tab. Would be meaningful only after step 2 produces a URL with ' +
      'query-string inputs.'});

    const currentUrl = page.url();
    test.info().annotations.push({type: 'info', description:
      `current URL after step 2: ${currentUrl}`});
    // Soft existence check — we can at least confirm a URL was set to the pk.ivp file path.
    expect(currentUrl.toLowerCase()).toContain('pk.ivp');
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
