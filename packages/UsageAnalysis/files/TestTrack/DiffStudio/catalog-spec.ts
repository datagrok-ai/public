import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('DiffStudio Catalog: PK-PD load, save to Model Hub, refresh, run, modify dose', async ({page}) => {
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

  await softStep('Step 1: Open PK-PD from Diff Studio Library', async () => {
    await page.evaluate(async () => {
      const f = DG.Func.find({name: 'runDiffStudio', package: 'DiffStudio'})[0];
      const call = f.prepare();
      await call.call();
      grok.shell.addView(call.getOutputParamValue());
    });
    await page.waitForFunction(() => grok.shell.v?.name === 'Diff Studio', null, {timeout: 30000});
    await page.waitForTimeout(2000);
    const card = page.locator('.diff-studio-hub-card', {hasText: 'PK-PD'}).first();
    await card.waitFor({timeout: 15000});
    await card.dblclick();
    await page.waitForFunction(() => grok.shell.v?.name === 'PK-PD', null, {timeout: 30000});
    await page.waitForTimeout(2000);
  });

  await softStep('Step 2: Click Save to Model Hub icon', async () => {
    await page.locator('.diff-studio-ribbon-save-to-model-catalog-icon').click();
    // Wait for balloon / notification
    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('.d4-balloon, .grok-notification'))
        .some(b => b.textContent?.includes('Saved to Library') || b.textContent?.includes('PK-PD.ivp')),
      null, {timeout: 15000}).catch(() => {});
    await page.waitForTimeout(2000);
  });

  await softStep('Step 3: Open Model Hub (Apps > Model Hub)', async () => {
    await page.evaluate(async () => {
      const f = DG.Func.find({package: 'Compute2', name: 'modelCatalog'})[0];
      const call = f.prepare();
      await call.call();
      const view = call.getOutputParamValue();
      if (view) grok.shell.addView(view);
    });
    await page.waitForFunction(() => grok.shell.v?.name === 'Model Hub', null, {timeout: 30000});
    await page.waitForTimeout(3000);
    const pkpdCount = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-list-item'))
        .filter(e => e.textContent?.trim() === 'PK-PD').length);
    expect(pkpdCount).toBeGreaterThan(0);
  });

  await softStep('Step 4: Click Refresh icon; catalog reloads', async () => {
    await page.evaluate(() => {
      const refresh = Array.from(document.querySelectorAll('.grok-icon'))
        .find(e => e.className.includes('fa-sync')) as HTMLElement;
      refresh?.click();
    });
    await page.waitForTimeout(3000);
    const pkpdCount = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-list-item'))
        .filter(e => e.textContent?.trim() === 'PK-PD').length);
    expect(pkpdCount).toBeGreaterThan(0);
  });

  await softStep('Step 5: Run PK-PD from catalog', async () => {
    await page.evaluate(() => {
      const items = Array.from(document.querySelectorAll('.d4-list-item'))
        .filter(e => e.textContent?.trim() === 'PK-PD') as HTMLElement[];
      const last = items[items.length - 1];
      last?.click();
      last?.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
    });
    await page.waitForFunction(() => grok.shell.v?.name === 'PK-PD', null, {timeout: 30000});
    await page.waitForTimeout(3000);
    const info = await page.evaluate(() => ({
      hasCount: !!document.querySelector('[name="input-host-count"]'),
      hasDose: !!document.querySelector('[name="input-host-dose"]'),
      hasMultiaxis: document.body.innerText.includes('Multiaxis'),
    }));
    expect(info.hasCount).toBe(true);
    expect(info.hasDose).toBe(true);
    expect(info.hasMultiaxis).toBe(true);
  });

  await softStep('Step 6: Modify dose input; verify live update and URL', async () => {
    const inp = page.locator('[name="input-host-dose"] input.ui-input-editor');
    const before = await inp.inputValue();
    await inp.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('5000');
    await page.keyboard.press('Tab');
    await page.waitForTimeout(2500);
    const after = await inp.inputValue();
    expect(after).toBe('5000');
    expect(after).not.toBe(before);
    expect(page.url()).toContain('dose=5000');
  });

  // Cleanup: delete the saved PK-PD file from the user's catalog
  await page.evaluate(async () => {
    try {
      const files = await grok.dapi.files.list('System:AppData/DiffStudio/', true, 'PK-PD.ivp');
      for (const f of files) {
        if (f.path?.includes('library/PK-PD')) continue;
        await grok.dapi.files.delete(f);
      }
    } catch (e) { /* ignore */ }
    grok.shell.closeAll();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
