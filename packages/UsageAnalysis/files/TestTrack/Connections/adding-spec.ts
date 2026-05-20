import {test, expect} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Connections / Adding', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Pre-cleanup: delete prior test_postgres / test_postgres_2 so this spec is idempotent.
  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    for (const n of ['test_postgres', 'test_postgres_2']) {
      const c = await g.dapi.connections.filter(`name = "${n}"`).first();
      if (c) await g.dapi.connections.delete(c);
    }
  });

  await softStep('Step 1: Go to Browse > Databases', async () => {
    await page.goto(`${baseUrl}/db?browse=db`, {waitUntil: 'networkidle', timeout: 30_000});
    await page.waitForTimeout(2500);
    await expect(page.locator('.d4-tree-view-group-label', {hasText: /^Postgres$/}).first()).toBeVisible({timeout: 15_000});
  });

  await softStep('Step 2: Right-click Postgres → New connection...', async () => {
    await page.evaluate(() => {
      const node = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find(el => el.textContent?.trim() === 'Postgres');
      node?.closest('.d4-tree-view-node')?.dispatchEvent(
        new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}),
      );
    });
    await page.waitForTimeout(500);
    await page.locator('[name="div-New-connection..."]').click();
    await page.waitForTimeout(1500);
    await expect(page.locator('.d4-dialog')).toBeVisible({timeout: 5000});
  });

  await softStep('Steps 3-4: Fill Name, Server, Port, Db, Login, Password', async () => {
    // Dart inputs only react to native key events — not .value setters or fill(). Use keyboard.
    const fillField = async (name: string, value: string) => {
      const loc = page.locator(`[name="${name}"]`);
      await loc.click();
      await page.keyboard.press('Control+A');
      await page.keyboard.press('Delete');
      await page.keyboard.type(value);
    };
    await fillField('input-Name', 'test_postgres');
    await fillField('input-Server', 'db.datagrok.ai');
    await fillField('input-Port', '54322');
    await fillField('input-Db', 'northwind');
    await fillField('input-Login', 'datagrok');
    await fillField('input-Password', 'placeholder');
    expect(await page.locator('[name="input-Name"]').inputValue()).toBe('test_postgres');
    expect(await page.locator('[name="input-Server"]').inputValue()).toBe('db.datagrok.ai');
    expect(await page.locator('[name="input-Port"]').inputValue()).toBe('54322');
    expect(await page.locator('[name="input-Db"]').inputValue()).toBe('northwind');
    expect(await page.locator('[name="input-Login"]').inputValue()).toBe('datagrok');
  });

  await softStep('Step 5: Click TEST — dialog stays open (auth error expected with placeholder)', async () => {
    await page.evaluate(() => {
      const btn = Array.from(document.querySelectorAll('.d4-dialog button'))
        .find(b => b.textContent?.trim() === 'TEST') as HTMLButtonElement | undefined;
      btn?.click();
    });
    await page.waitForTimeout(12_000);
    // The TEST run is async; with a placeholder password the dialog stays open and OK is enabled.
    // Assert the dialog is still here so we know the click was wired up and the test attempt ran.
    await expect(page.locator('.d4-dialog')).toBeVisible();
  });

  await softStep('Step 6: Click OK — connection saved', async () => {
    await page.evaluate(() => {
      const ok = Array.from(document.querySelectorAll('.d4-dialog button'))
        .find(b => b.textContent?.trim() === 'OK') as HTMLButtonElement | undefined;
      ok?.click();
    });
    await page.waitForTimeout(3000);
    await expect(page.locator('.d4-dialog')).toHaveCount(0);
    const exists = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const c = await g.dapi.connections.filter('name = "test_postgres"').first();
      return !!c;
    });
    expect(exists).toBe(true);
  });

  await softStep('Step 7: Create test_postgres_2 from the NEW CONNECTION... button', async () => {
    await page.evaluate(() => {
      const btn = Array.from(document.querySelectorAll('button.ui-btn-ok'))
        .find(b => b.textContent?.trim().toLowerCase() === 'new connection...') as HTMLButtonElement | undefined;
      btn?.click();
    });
    await page.waitForTimeout(1500);
    await expect(page.locator('.d4-dialog')).toBeVisible({timeout: 5000});

    // Pick Postgres data source
    await page.evaluate(() => {
      const ds = document.querySelector('[name="input-Data-Source"]') as HTMLSelectElement | null;
      if (!ds) return;
      ds.value = 'Postgres';
      ds.dispatchEvent(new Event('input', {bubbles: true}));
      ds.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.waitForTimeout(1200);

    const fillField = async (name: string, value: string) => {
      const loc = page.locator(`[name="${name}"]`);
      await loc.click();
      await page.keyboard.press('Control+A');
      await page.keyboard.press('Delete');
      await page.keyboard.type(value);
    };
    await fillField('input-Name', 'test_postgres_2');
    await fillField('input-Server', 'db.datagrok.ai');
    await fillField('input-Port', '54322');
    await fillField('input-Db', 'northwind');
    await fillField('input-Login', 'datagrok');
    await fillField('input-Password', 'placeholder');

    await page.evaluate(() => {
      const ok = Array.from(document.querySelectorAll('.d4-dialog button'))
        .find(b => b.textContent?.trim() === 'OK') as HTMLButtonElement | undefined;
      ok?.click();
    });
    await page.waitForTimeout(3000);
    await expect(page.locator('.d4-dialog')).toHaveCount(0);
    const exists = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const c = await g.dapi.connections.filter('name = "test_postgres_2"').first();
      return !!c;
    });
    expect(exists).toBe(true);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
