import {test, expect} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Prerequisite: Adding-spec.ts must have run to create test_postgres

test('Connections / Edit', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await page.goto(`${baseUrl}/connect?browse=connections`, {waitUntil: 'networkidle', timeout: 30000});
  await page.waitForTimeout(3000);
  await page.evaluate(async () => {
    const postgresNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
      .find((el) => el.textContent?.trim() === 'Postgres');
    (postgresNode as HTMLElement | undefined)?.click();
    await new Promise((r) => setTimeout(r, 2000));
  });

  await softStep('Step 1: Reload the tree', async () => {
    await page.evaluate(async () => {
      const sync = document.querySelector('[name="icon-sync"]') as HTMLElement | null;
      sync?.click();
      await new Promise((r) => setTimeout(r, 1500));
    });
    await expect(page.locator('[data-link="/db/null.test_postgres"]').or(
      page.locator('[name="span-test-postgres"]'))).toBeVisible({timeout: 5000});
  });

  await softStep('Steps 2-4: Right-click test_postgres → Edit → rename to new_test_postgres', async () => {
    await page.evaluate(async () => {
      // Locate by data-link or name to avoid ambiguity if multiple connections share friendlyName
      const target = (document.querySelector('[data-link="/db/null.test_postgres"]')
        ?? document.querySelector('[name="span-test-postgres"]')
        ?? Array.from(document.querySelectorAll('.d4-link-label, label'))
            .filter((el) => el.textContent?.trim() === 'test_postgres').pop()) as HTMLElement | undefined;
      if (!target) return;
      target.scrollIntoView({behavior: 'instant' as ScrollBehavior, block: 'center'});
      const r = target.getBoundingClientRect();
      target.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: r.left + 10, clientY: r.top + 5,
      }));
      await new Promise((r) => setTimeout(r, 500));
      const edit = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find((el) => el.textContent?.trim() === 'Edit...') as HTMLElement | undefined;
      edit?.click();
      await new Promise((r) => setTimeout(r, 1500));
    });
    await expect(page.locator('text=Edit Connection')).toBeVisible({timeout: 5000});
    // Wait for the dialog form to render its inputs (rendered async)
    await page.waitForFunction(
      () => (document.querySelector('.d4-dialog')?.querySelectorAll('input')?.length ?? 0) >= 5,
      null, {timeout: 10000});

    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog') as HTMLElement;
      const nameInput = dialog.querySelector('[name="input-Name"]') as HTMLInputElement;
      const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      ns.call(nameInput, 'new_test_postgres');
      nameInput.dispatchEvent(new Event('input', {bubbles: true}));
      nameInput.dispatchEvent(new Event('change', {bubbles: true}));
      const ok = (dialog.querySelector('[name="button-OK"]')
        ?? Array.from(dialog.querySelectorAll('button, .ui-btn'))
            .find((b) => b.textContent?.trim() === 'OK')) as HTMLElement;
      ok.click();
    });
    await page.waitForTimeout(3000);

    // Verify via JS API — DOM may show duplicate friendlyName labels
    const renamed = await page.evaluate(async () => {
      const c = await (window as any).grok.dapi.connections.find('af9bcf40-21a0-11f1-89e2-7b1321b80948');
      return c?.friendlyName;
    });
    expect(renamed).toBe('new_test_postgres');
  });

  await softStep('Steps 5-6: Set wrong credentials, verify Test connection fails', async () => {
    await page.evaluate(async () => {
      const sync = document.querySelector('[name="icon-sync"]') as HTMLElement | null;
      sync?.click();
      await new Promise((r) => setTimeout(r, 1500));
    });

    await page.evaluate(async () => {
      const target = document.querySelector('[data-link="/db/null.test_postgres"]') as HTMLElement | null;
      if (!target) return;
      target.scrollIntoView({behavior: 'instant' as ScrollBehavior, block: 'center'});
      const r = target.getBoundingClientRect();
      target.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: r.left + 10, clientY: r.top + 5,
      }));
      await new Promise((r) => setTimeout(r, 500));
      const edit = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find((el) => el.textContent?.trim() === 'Edit...') as HTMLElement | undefined;
      edit?.click();
    });
    await expect(page.locator('text=Edit Connection')).toBeVisible({timeout: 5000});
    await page.waitForFunction(
      () => (document.querySelector('.d4-dialog')?.querySelectorAll('input')?.length ?? 0) >= 5,
      null, {timeout: 10000});

    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog') as HTMLElement;
      const setInput = (selector: string, value: string) => {
        const inp = dialog.querySelector(selector) as HTMLInputElement;
        const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        ns.call(inp, value);
        inp.dispatchEvent(new Event('input', {bubbles: true}));
        inp.dispatchEvent(new Event('change', {bubbles: true}));
      };
      setInput('[name="input-Login"]', 'wronguser');
      setInput('[name="input-Password"]', 'wrongpassword');
      const ok = (dialog.querySelector('[name="button-OK"]')
        ?? Array.from(dialog.querySelectorAll('button, .ui-btn'))
            .find((b) => b.textContent?.trim() === 'OK')) as HTMLElement;
      ok.click();
    });
    await page.waitForTimeout(2500);

    // Trigger Test connection via JS API on the entity (UI balloon may auto-dismiss)
    const result = await page.evaluate(async () => {
      const conn = await (window as any).grok.dapi.connections.find('af9bcf40-21a0-11f1-89e2-7b1321b80948');
      return conn ? await conn.test() : null;
    });
    expect(result).toMatch(/password authentication failed for user "wronguser"/);
  });

  // Step 7: Restoring with correct credentials requires the real DB password
  // (managed by DevOps). Skipped — see edit-run.md.

  // Cleanup: restore friendlyName for the next run
  await page.evaluate(async () => {
    const conn = await (window as any).grok.dapi.connections.find('af9bcf40-21a0-11f1-89e2-7b1321b80948');
    if (conn) {
      conn.friendlyName = 'test_postgres';
      await (window as any).grok.dapi.connections.save(conn);
    }
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
