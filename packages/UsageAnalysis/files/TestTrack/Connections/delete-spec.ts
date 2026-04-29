import {test, expect, Page} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Prerequisite: Adding-spec.ts and Edit-spec.ts must have run
// (new_test_postgres and test_postgres_2 must exist)

const deleteConnection = async (page: Page, name: string) => {
  await page.evaluate(async (connName) => {
    const target = Array.from(document.querySelectorAll('.d4-link-label, label'))
      .find((el) => el.textContent?.trim() === connName) as HTMLElement | undefined;
    if (!target) return;
    target.scrollIntoView({behavior: 'instant' as ScrollBehavior, block: 'center'});
    await new Promise((r) => setTimeout(r, 200));
    target.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
      clientX: target.getBoundingClientRect().left + 10,
      clientY: target.getBoundingClientRect().top + 5,
    }));
    await new Promise((r) => setTimeout(r, 400));
    const deleteItem = Array.from(document.querySelectorAll('.d4-menu-item'))
      .find((el) => el.textContent?.trim() === 'Delete...') as HTMLElement | undefined;
    deleteItem?.click();
    await new Promise((r) => setTimeout(r, 1000));
    const deleteBtn = Array.from(document.querySelectorAll('button,.ui-btn'))
      .find((b) => b.textContent?.trim() === 'DELETE') as HTMLElement | undefined;
    deleteBtn?.click();
    await new Promise((r) => setTimeout(r, 2000));
  }, name);

  await page.evaluate(() => {
    (document.querySelector('.d4-refresh, .fa-sync, .fa-sync-alt') as HTMLElement | null)?.click();
  });
  await page.waitForTimeout(2000);
};

test('Connections / Delete', async ({page}) => {
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

  await softStep('Steps 1-2: Delete new_test_postgres', async () => {
    await deleteConnection(page, 'new_test_postgres');
    const present = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-link-label, label'))
        .some((el) => el.textContent?.trim() === 'new_test_postgres'),
    );
    expect(present).toBe(false);
  });

  await softStep('Steps 3-4: Delete test_postgres_2', async () => {
    await deleteConnection(page, 'test_postgres_2');
    const present = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-link-label, label'))
        .some((el) => el.textContent?.trim() === 'test_postgres_2'),
    );
    expect(present).toBe(false);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
