import {test, expect} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Connections / Sparql', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  // Pre-clean any leftover test_sparql from prior runs.
  await page.evaluate(async () => {
    try {
      const c = await (window as any).grok.dapi.connections.filter('name = "test_sparql"').first();
      if (c) await (window as any).grok.dapi.connections.delete(c);
    } catch (_) {}
  });

  await page.goto(`${baseUrl}/connect?browse=connections`, {waitUntil: 'networkidle', timeout: 30000});
  await page.waitForTimeout(3000);

  await softStep('Step 1: Browse > Databases, press ellipsis to list providers without connections', async () => {
    await expect(page.locator('text=Databases').first()).toBeVisible({timeout: 10000});
    // Click the ellipsis *inside* the Databases group (page has multiple ellipses).
    await page.evaluate(async () => {
      const databasesNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((l) => l.textContent?.trim() === 'Databases')?.closest('.d4-tree-view-group') as HTMLElement | undefined;
      const ellipsis = databasesNode?.querySelector('[name="icon-ellipsis-h"]') as HTMLElement | null;
      ellipsis?.click();
      await new Promise((r) => setTimeout(r, 1500));
    });
    const sparqlVisible = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .some((l) => l.textContent?.trim() === 'Sparql'),
    );
    expect(sparqlVisible).toBe(true);
  });

  await softStep('Step 2: Right-click Sparql, choose New connection...', async () => {
    await page.evaluate(async () => {
      const sparql = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((l) => l.textContent?.trim() === 'Sparql') as HTMLElement | undefined;
      sparql?.scrollIntoView({block: 'center'});
      sparql?.closest('.d4-tree-view-node')?.dispatchEvent(
        new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}),
      );
      await new Promise((r) => setTimeout(r, 500));
      const item = document.querySelector('[name="div-New-connection..."]') as HTMLElement | null;
      item?.click();
      await new Promise((r) => setTimeout(r, 1000));
    });
    await expect(page.locator('.d4-dialog')).toBeVisible({timeout: 5000});
    await expect(page.locator('[name="input-Name"]')).toBeVisible();
  });

  await softStep('Step 3: Enter test_sparql to the Name field', async () => {
    await page.locator('[name="input-Name"]').click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('test_sparql');
    const v = await page.locator('[name="input-Name"]').inputValue();
    expect(v).toBe('test_sparql');
  });

  await softStep('Step 4: Endpoint + Requires Server=true, leave Prefixes empty', async () => {
    await page.locator('[name="input-Endpoint"]').click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('http://data.ontotext.com/repositories/data-last');

    await page.evaluate(() => {
      const cb = document.querySelector('[name="input-Requires-Server"]') as HTMLInputElement | null;
      if (cb && !cb.checked) cb.click();
    });

    const state = await page.evaluate(() => ({
      endpoint: (document.querySelector('[name="input-Endpoint"]') as HTMLInputElement)?.value,
      requiresServer: (document.querySelector('[name="input-Requires-Server"]') as HTMLInputElement)?.checked,
      prefixes: (document.querySelector('[name="Prefixes"]') as HTMLTextAreaElement)?.value ?? '',
    }));
    expect(state.endpoint).toBe('http://data.ontotext.com/repositories/data-last');
    expect(state.requiresServer).toBe(true);
    expect(state.prefixes).toBe('');
  });

  await softStep('Step 5: Click TEST — expect "connected successfully"', async () => {
    await page.locator('.d4-dialog [name="button-TEST"]').click();
    const balloonText = await page.evaluate(async () => {
      for (let i = 0; i < 60; i++) {
        const b = Array.from(document.querySelectorAll('.d4-balloon, .grok-balloon'))
          .find((el) => (el as HTMLElement).offsetParent !== null) as HTMLElement | undefined;
        if (b) return b.innerText.trim();
        await new Promise((r) => setTimeout(r, 500));
      }
      return '';
    });
    // Asserts scenario intent — fails on dev where DNS for data.ontotext.com is blocked.
    expect(balloonText).toContain('connected successfully');
  });

  await softStep('Step 6: Click OK — connection saved', async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForTimeout(1500);
    const saved = await page.evaluate(async () => {
      const c = await (window as any).grok.dapi.connections.filter('name = "test_sparql"').first();
      return !!c;
    });
    expect(saved).toBe(true);
  });

  await softStep('Step 7: Right-click test_sparql, Delete, confirm DELETE', async () => {
    // Refresh + re-expose Sparql + expand it so the freshly-saved connection shows.
    await page.evaluate(async () => {
      (document.querySelector('[name="icon-sync"], [name="icon-refresh"]') as HTMLElement | null)?.click();
      await new Promise((r) => setTimeout(r, 1500));

      const labels = () => Array.from(document.querySelectorAll('.d4-tree-view-group-label'));
      let sparql = labels().find((l) => l.textContent?.trim() === 'Sparql') as HTMLElement | undefined;
      if (!sparql) {
        const databasesNode = labels().find((l) => l.textContent?.trim() === 'Databases')
          ?.closest('.d4-tree-view-group') as HTMLElement | undefined;
        (databasesNode?.querySelector('[name="icon-ellipsis-h"]') as HTMLElement | null)?.click();
        await new Promise((r) => setTimeout(r, 1500));
        sparql = labels().find((l) => l.textContent?.trim() === 'Sparql') as HTMLElement | undefined;
      }
      const sparqlNode = sparql?.closest('.d4-tree-view-node') as HTMLElement | null;
      // Click triangle, then dblclick label as a fallback to ensure expansion.
      (sparqlNode?.querySelector('.d4-tree-view-tri') as HTMLElement | null)?.click();
      await new Promise((r) => setTimeout(r, 800));
      sparql?.dispatchEvent(new MouseEvent('dblclick', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 2500));
    });

    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label, .d4-link-label, label'))
        .some((l) => l.textContent?.trim() === 'test_sparql'),
    {timeout: 15000});

    await page.evaluate(async () => {
      const conn = Array.from(document.querySelectorAll(
        '.d4-tree-view-group-label, .d4-tree-view-item-label, .d4-link-label, label',
      )).find((l) => l.textContent?.trim() === 'test_sparql') as HTMLElement | undefined;
      conn?.scrollIntoView({block: 'center'});
      const target = (conn?.closest('.d4-tree-view-node') ?? conn) as HTMLElement | null;
      target?.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
      await new Promise((r) => setTimeout(r, 600));
      (document.querySelector('[name="div-Delete..."]') as HTMLElement | null)?.click();
      await new Promise((r) => setTimeout(r, 1000));
      (document.querySelector('.d4-dialog [name="button-DELETE"]') as HTMLElement | null)?.click();
      await new Promise((r) => setTimeout(r, 2000));
    });

    const stillExists = await page.evaluate(async () => {
      const c = await (window as any).grok.dapi.connections.filter('name = "test_sparql"').first()
        .catch(() => null);
      return !!c;
    });
    expect(stillExists).toBe(false);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
