import {test, expect} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Prerequisite: Edit-spec.ts must have run (new_test_postgres must exist)

test('Connections / Browser', async ({page}) => {
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

  await softStep('Steps 1-3: Search for new_test connection', async () => {
    const searchInput = page.locator('input[placeholder="Search connections by name or by #tags"]');
    await searchInput.fill('new_test');
    await page.waitForTimeout(1000);
    await expect(page.locator('text=new_test_postgres').first()).toBeVisible({timeout: 5000});
    await expect(page.locator('text=1 / 1')).toBeVisible({timeout: 3000});
  });

  await softStep('Steps 4-5: Click connection, verify Context Pane tabs', async () => {
    await page.evaluate(async () => {
      const link = Array.from(document.querySelectorAll('.d4-link-label, label'))
        .find((el) => el.textContent?.trim() === 'new_test_postgres');
      (link as HTMLElement | undefined)?.click();
      await new Promise((r) => setTimeout(r, 2000));
    });

    await expect(page.locator('text=Details').first()).toBeVisible({timeout: 5000});

    await page.evaluate(async () => {
      const header = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'))
        .find((el) => el.textContent?.includes('Activity'));
      (header as HTMLElement | undefined)?.click();
      await new Promise((r) => setTimeout(r, 500));
    });
    const activityContent = await page.evaluate(() => {
      const header = Array.from(document.querySelectorAll('.d4-accordion-pane-header, .grok-accordion-pane-header'))
        .find((el) => el.textContent?.includes('Activity'));
      return header?.nextElementSibling?.textContent || '';
    });
    expect(activityContent).toContain('created');
  });

  await softStep('Step 6: Click dropdown icon near connection name', async () => {
    const menuItems = await page.evaluate(async () => {
      const arrow = document.querySelector('.grok-context-arrow-down') as HTMLElement | null;
      if (!arrow) return [] as string[];
      arrow.click();
      await new Promise((r) => setTimeout(r, 500));
      const items = Array.from(document.querySelectorAll('.d4-menu-item'))
        .map((el) => el.textContent?.trim()).filter((t): t is string => !!t);
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      return items;
    });
    expect(menuItems).toContain('Edit...');
    expect(menuItems).toContain('Delete...');
    expect(menuItems).toContain('Browse queries');
    expect(menuItems).toContain('Test connection');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
