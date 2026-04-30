import {test, expect} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Connections / Schema', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await page.goto(`${baseUrl}/connect?browse=connections`, {waitUntil: 'networkidle', timeout: 30000});
  await page.waitForTimeout(3000);

  await softStep('Steps 1-2: Navigate to Databases > Postgres > Northwind', async () => {
    await page.evaluate(async () => {
      const postgresNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((el) => el.textContent?.trim() === 'Postgres');
      const row = postgresNode?.closest('.d4-tree-view-node');
      const tri = row?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri') as HTMLElement | null;
      if (tri) tri.click();
      else (postgresNode as HTMLElement | undefined)?.click();
      await new Promise((r) => setTimeout(r, 3000));
    });
    await expect(page.locator('text=Northwind').first()).toBeVisible({timeout: 10000});
  });

  await softStep('Step 3: Expand Northwind Schemas (Browse schema alternative)', async () => {
    await page.evaluate(async () => {
      const northwindNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((el) => el.textContent?.trim() === 'Northwind');
      const row = northwindNode?.closest('.d4-tree-view-node');
      const tri = row?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri') as HTMLElement | null;
      if (tri) tri.click();
      else (northwindNode as HTMLElement | undefined)?.click();
      await new Promise((r) => setTimeout(r, 3000));
    });
    const hasSchemas = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .some((el) => el.textContent?.trim() === 'Schemas'),
    );
    expect(hasSchemas).toBe(true);
  });

  await softStep('Step 4: Tables have DB interaction context menus', async () => {
    await page.evaluate(async () => {
      const schemasNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((el) => el.textContent?.trim() === 'Schemas');
      const sRow = schemasNode?.closest('.d4-tree-view-node');
      (sRow?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri') as HTMLElement | null)?.click();
      await new Promise((r) => setTimeout(r, 2000));
      const publicNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((el) => el.textContent?.trim() === 'public');
      const pRow = publicNode?.closest('.d4-tree-view-node');
      (pRow?.querySelector('.d4-tree-view-group-tri, .d4-tree-view-tri') as HTMLElement | null)?.click();
      await new Promise((r) => setTimeout(r, 3000));
    });

    const hasTables = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .some((el) => el.textContent?.trim() === 'customers'),
    );
    expect(hasTables).toBe(true);

    const menuItems = await page.evaluate(async () => {
      const target = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((el) => el.textContent?.trim() === 'customers');
      if (!target) return [] as string[];
      (target.closest('.d4-tree-view-node') as HTMLElement)?.dispatchEvent(
        new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}),
      );
      await new Promise((r) => setTimeout(r, 400));
      return Array.from(document.querySelectorAll('.d4-menu-item'))
        .map((e) => e.textContent?.trim())
        .filter((t): t is string => !!t);
    });
    expect(menuItems).toContain('Get All');
    expect(menuItems).toContain('Get Top 100');
    expect(menuItems).toContain('New SQL Query...');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
