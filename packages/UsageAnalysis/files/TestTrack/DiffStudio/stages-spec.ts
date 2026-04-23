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

test('DiffStudio Stages (Acid Production): Load, Multiaxis+Facet, modify input, tooltips', async ({page}) => {
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

  await softStep('Step 1: Open Diff Studio + Acid Production from Library', async () => {
    await page.evaluate(async () => {
      const f = DG.Func.find({name: 'runDiffStudio', package: 'DiffStudio'})[0];
      const call = f.prepare();
      await call.call();
      grok.shell.addView(call.getOutputParamValue());
    });
    await page.waitForFunction(() => grok.shell.v?.name === 'Diff Studio', null, {timeout: 30000});
    await page.waitForTimeout(2000);
    const card = page.locator('.diff-studio-hub-card', {hasText: 'Acid Production'}).first();
    await card.waitFor({timeout: 15000});
    await card.dblclick();
    // The library card opens a view named 'GA-production' (not 'Acid Production')
    await page.waitForFunction(() =>
      grok.shell.v?.name === 'GA-production' || grok.shell.v?.name === 'Acid Production',
      null, {timeout: 30000});
    await page.waitForTimeout(2000);
    const inputs = await page.locator('[name^="input-host-"]').count();
    expect(inputs).toBeGreaterThan(10);
  });

  await softStep('Step 2: Multiaxis and Facet tabs updated', async () => {
    const tabs = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.tab-handle')).map(t => t.textContent?.trim()));
    expect(tabs).toEqual(expect.arrayContaining(['Multiaxis', 'Facet']));
  });

  await softStep('Step 3: Modify 1-st stage input; URL updates live', async () => {
    const inp = page.locator('[name="input-host-1-st-stage"] input.ui-input-editor');
    const before = await inp.inputValue();
    await inp.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('50');
    await page.keyboard.press('Tab');
    await page.waitForTimeout(2000);
    const after = await inp.inputValue();
    expect(after).toBe('50');
    expect(after).not.toBe(before);
    expect(page.url()).toContain('1-ststage=50');
  });

  await softStep('Step 4: Tooltips on inputs (1-st stage, biomass, glucose)', async () => {
    const tooltips = await page.evaluate(async () => {
      const names = ['1-st-stage', 'biomass', 'glucose'];
      const results: Record<string, string> = {};
      for (const n of names) {
        const host = document.querySelector(`[name="input-host-${n}"]`);
        const label = host?.querySelector('label') as HTMLElement;
        const rect = label?.getBoundingClientRect();
        label?.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: rect.left + 5, clientY: rect.top + 5}));
        label?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise(r => setTimeout(r, 900));
        results[n] = document.querySelector('.d4-tooltip')?.textContent?.trim() ?? '';
        label?.dispatchEvent(new MouseEvent('mouseleave', {bubbles: true}));
        await new Promise(r => setTimeout(r, 300));
      }
      return results;
    });
    expect(tooltips['1-st-stage'].length).toBeGreaterThan(0);
    expect(tooltips.biomass.length).toBeGreaterThan(0);
    expect(tooltips.glucose.length).toBeGreaterThan(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
