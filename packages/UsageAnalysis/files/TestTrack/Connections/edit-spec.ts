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

  await softStep('Steps 2-4: Edit test_postgres, rename to new_test_postgres', async () => {
    await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('.d4-link-label, label'))
        .filter((el) => el.textContent?.trim() === 'test_postgres').pop() as HTMLElement | undefined;
      if (!target) return;
      target.scrollIntoView({behavior: 'instant' as ScrollBehavior, block: 'center'});
      target.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: target.getBoundingClientRect().left + 10,
        clientY: target.getBoundingClientRect().top + 5,
      }));
      await new Promise((r) => setTimeout(r, 400));
      const edit = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find((el) => el.textContent?.trim() === 'Edit...') as HTMLElement | undefined;
      edit?.click();
      await new Promise((r) => setTimeout(r, 1000));
    });
    await expect(page.locator('text=Edit Connection')).toBeVisible({timeout: 5000});

    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog')!;
      const labels = Array.from(dialog.querySelectorAll('label, .d4-label, .ui-label'));
      const nameLabel = labels.find((l) => l.textContent?.trim() === 'Name');
      const nameInput = nameLabel?.closest('.d4-input-base, .d4-flex-row, div')?.querySelector('input') as HTMLInputElement | null;
      if (nameInput) {
        const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        ns.call(nameInput, 'new_test_postgres');
        nameInput.dispatchEvent(new Event('input', {bubbles: true}));
        nameInput.dispatchEvent(new Event('change', {bubbles: true}));
      }
      const ok = Array.from(document.querySelectorAll('button,.ui-btn'))
        .find((b) => b.textContent?.trim() === 'OK') as HTMLElement | undefined;
      ok?.click();
    });
    await page.waitForTimeout(2000);
    await expect(page.locator('text=new_test_postgres').first()).toBeVisible({timeout: 5000});
  });

  await softStep('Steps 5-6: Set wrong credentials, verify test fails', async () => {
    await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('.d4-link-label, label'))
        .find((el) => el.textContent?.trim() === 'new_test_postgres') as HTMLElement | undefined;
      if (!target) return;
      target.scrollIntoView({behavior: 'instant' as ScrollBehavior, block: 'center'});
      target.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: target.getBoundingClientRect().left + 10,
        clientY: target.getBoundingClientRect().top + 5,
      }));
      await new Promise((r) => setTimeout(r, 400));
      const edit = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find((el) => el.textContent?.trim() === 'Edit...') as HTMLElement | undefined;
      edit?.click();
      await new Promise((r) => setTimeout(r, 1000));
      const dialog = document.querySelector('.d4-dialog')!;
      const labels = Array.from(dialog.querySelectorAll('label, .d4-label, .ui-label'));
      const setInput = (labelText: string, value: string) => {
        const label = labels.find((l) => l.textContent?.trim() === labelText);
        const inp = label?.closest('.d4-input-base, .d4-flex-row, div')?.querySelector('input') as HTMLInputElement | null;
        if (!inp) return;
        const ns = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
        ns.call(inp, value);
        inp.dispatchEvent(new Event('input', {bubbles: true}));
        inp.dispatchEvent(new Event('change', {bubbles: true}));
      };
      setInput('Login', 'wronguser');
      setInput('Password', 'wrongpassword');
      const ok = Array.from(document.querySelectorAll('button,.ui-btn'))
        .find((b) => b.textContent?.trim() === 'OK') as HTMLElement | undefined;
      ok?.click();
    });
    await page.waitForTimeout(2000);

    await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('.d4-link-label, label'))
        .find((el) => el.textContent?.trim() === 'new_test_postgres') as HTMLElement | undefined;
      if (!target) return;
      target.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: target.getBoundingClientRect().left + 10,
        clientY: target.getBoundingClientRect().top + 5,
      }));
      await new Promise((r) => setTimeout(r, 400));
      const testConn = Array.from(document.querySelectorAll('.d4-menu-item'))
        .find((el) => el.textContent?.trim() === 'Test connection') as HTMLElement | undefined;
      testConn?.click();
    });
    await page.waitForTimeout(12000);
    await expect(page.locator('text=wronguser').first()).toBeVisible({timeout: 15000});
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
