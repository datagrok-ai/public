import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('DiffStudio Scripting: Edit toggle, </> JS view, Run, Save with //tags: model, Model Hub', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);
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

  await softStep('Setup: Open Diff Studio + Bioreactor', async () => {
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
    await page.waitForTimeout(3000);
  });

  await softStep('Step 1: Turn on Edit toggle; click </> to open JS script view', async () => {
    await page.evaluate(() => {
      const editRibbonItem = Array.from(document.querySelectorAll('.d4-ribbon-item'))
        .find(el => el.textContent?.trim() === 'Edit') as HTMLElement;
      const editor = editRibbonItem?.querySelector('.ui-input-bool-switch .ui-input-editor') as HTMLElement;
      editor?.click();
    });
    await page.waitForTimeout(2000);
    const switchOn = await page.evaluate(() =>
      document.querySelector('.d4-ribbon-item .ui-input-switch')?.classList.contains('ui-input-switch-on'));
    expect(switchOn).toBe(true);

    await page.evaluate(() => {
      const span = Array.from(document.querySelectorAll('.d4-ribbon-name'))
        .find(el => el.textContent?.trim() === '</>') as HTMLElement;
      const rect = span.getBoundingClientRect();
      span.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: rect.left+5, clientY: rect.top+5}));
      span.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: rect.left+5, clientY: rect.top+5}));
      span.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: rect.left+5, clientY: rect.top+5}));
    });
    await page.waitForFunction(() => grok.shell.v?.type === 'ScriptView', null, {timeout: 15000});
    await page.waitForTimeout(1500);
    const info = await page.evaluate(() => ({
      type: grok.shell.v?.type,
      hasCM: document.querySelectorAll('.CodeMirror').length,
    }));
    expect(info.type).toBe('ScriptView');
    expect(info.hasCM).toBeGreaterThan(0);
  });

  await softStep('Step 2: Run the script; adjust Final input; live update', async () => {
    await page.locator('[name="icon-play"]').click();
    await page.waitForSelector('[name="input-host-Final"] input', {timeout: 30000});
    await page.waitForTimeout(2000);
    const finalInput = page.locator('[name="input-host-Final"] input');
    const before = await finalInput.inputValue();
    await finalInput.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('500');
    await page.keyboard.press('Tab');
    await page.waitForTimeout(2500);
    const after = await finalInput.inputValue();
    expect(after).toBe('500');
    expect(after).not.toBe(before);
    const missingFacet = await page.evaluate(() =>
      !document.body.innerText.includes('Facet') && !document.querySelector('[name="input-host-Process-mode"]'));
    expect(missingFacet).toBe(true);
  });

  await softStep('Step 3: Add //tags: model; Save the script', async () => {
    await page.evaluate(() => {
      const scriptView = Array.from(grok.shell.views).find(v => v.type === 'ScriptView');
      if (scriptView) grok.shell.v = scriptView;
    });
    await page.waitForTimeout(2000);
    const tagAdded = await page.evaluate(() => {
      const cm = (document.querySelector('.CodeMirror') as any)?.CodeMirror;
      if (!cm) return false;
      const text = cm.getValue();
      if (text.includes('//tags: model')) return true;
      cm.setValue(text.replace('//language: javascript', '//language: javascript\n//tags: model'));
      return cm.getValue().includes('//tags: model');
    });
    expect(tagAdded).toBe(true);
    await page.locator('[name="button-Save"]').click();
    await page.waitForTimeout(3000);
  });

  await softStep('Step 4: Access Model in Model Hub (Compute2:modelCatalog)', async () => {
    await page.evaluate(async () => {
      const f = DG.Func.find({package: 'Compute2', name: 'modelCatalog'})[0];
      const call = f.prepare();
      await call.call();
      const view = call.getOutputParamValue();
      if (view) grok.shell.addView(view);
    });
    await page.waitForFunction(() => grok.shell.v?.name === 'Model Hub', null, {timeout: 30000});
    await page.waitForTimeout(3000);
    const hasBioreactor = await page.evaluate(() =>
      !!Array.from(document.querySelectorAll('.d4-list-item')).find(el => el.textContent?.trim() === 'Bioreactor'));
    expect(hasBioreactor).toBe(true);
  });

  await softStep('Step 5: Interact with saved model; adjust Final', async () => {
    await page.evaluate(async () => {
      const items = Array.from(document.querySelectorAll('.d4-list-item'))
        .filter(el => el.textContent?.trim() === 'Bioreactor') as HTMLElement[];
      const card = items[items.length - 1];
      card.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
    });
    await page.waitForSelector('[name="input-host-Final"] input', {timeout: 30000});
    await page.waitForTimeout(2000);
    const inp = page.locator('[name="input-host-Final"] input');
    await inp.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('800');
    await page.keyboard.press('Tab');
    await page.waitForTimeout(2500);
    const after = await inp.inputValue();
    expect(after).toBe('800');
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
