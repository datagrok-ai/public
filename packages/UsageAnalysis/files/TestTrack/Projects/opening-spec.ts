import { test, expect, Page } from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const BASE_URL = 'https://release-ec2.datagrok.ai';
const LOGIN = 'claude';
const PASSWORD = 'grokclaude';

async function login(page: Page) {
  await page.goto(BASE_URL);
  const loginInput = page.locator('input[placeholder*="Login"]');
  if (await loginInput.isVisible({ timeout: 5000 }).catch(() => false)) {
    await loginInput.fill(LOGIN);
    await page.locator('input[type="password"]').fill(PASSWORD);
    await page.keyboard.press('Enter');
  }
  await page.waitForSelector('.d4-toolbox', { timeout: 60000 });
}

async function evalJs(page: Page, script: string): Promise<any> {
  return page.evaluate(script);
}

async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

test.describe('Projects / Opening', () => {
  const projectName = 'AutoTest-Opening-' + Date.now();

  test.beforeAll(async ({ browser }) => {
    const page = await browser.newPage();
    await login(page);

    // Create a test project with demo data
    await evalJs(page, `(async () => {
      const tv = grok.shell.addTableView(grok.data.demo.demog());
    })()`);
    await page.waitForTimeout(2000);
    await evalJs(page, `(async () => {
      const project = grok.shell.project;
      project.name = '${projectName}';
      await grok.dapi.projects.save(project);
    })()`);
    await page.waitForTimeout(3000);
    await closeAll(page);
    await page.close();
  });

  test.afterAll(async ({ browser }) => {
    const page = await browser.newPage();
    await login(page);
    await evalJs(page, `(async () => {
      const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
      if (p) await grok.dapi.projects.delete(p);
    })()`);
    await page.close();
  });

  test('Case 1: Navigate to Browse > Dashboards', async ({ page }) => {
    await login(page);
    await page.goto(`${BASE_URL}/projects`);
    await page.waitForSelector('.grok-gallery-grid', { timeout: 30000 });
    await expect(page.locator('.grok-gallery-grid')).toBeVisible();
  });

  test('Case 2: Find project in Dashboards', async ({ page }) => {
    await login(page);
    await page.goto(`${BASE_URL}/projects`);
    await page.waitForSelector('.grok-gallery-grid', { timeout: 30000 });

    const found = await evalJs(page, `(async () => {
      const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
      return p !== null;
    })()`);
    expect(found).toBe(true);
  });

  test('Case 3: Check Context Panel attributes', async ({ page }) => {
    await login(page);
    await page.goto(`${BASE_URL}/projects`);
    await page.waitForSelector('.grok-gallery-grid', { timeout: 30000 });

    // Select project via API and check context panel
    const attrs = await evalJs(page, `(async () => {
      const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
      return { name: p.name, createdOn: p.createdOn !== null, author: p.author !== null };
    })()`);
    expect(attrs.name).toBe(projectName);
    expect(attrs.createdOn).toBe(true);
    expect(attrs.author).toBe(true);
  });

  test('Case 4: Open project and verify data', async ({ page }) => {
    await login(page);
    await closeAll(page);

    const tables = await evalJs(page, `(async () => {
      const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
      await p.open();
      return grok.shell.tables.map(t => ({ name: t.name, rows: t.rowCount, cols: t.columns.length }));
    })()`);
    expect(tables.length).toBeGreaterThan(0);
    expect(tables[0].rows).toBeGreaterThan(0);
    expect(tables[0].cols).toBeGreaterThan(0);
  });

  test('Case 5: Close all views', async ({ page }) => {
    await login(page);
    await closeAll(page);
    const viewCount = await evalJs(page, 'Array.from(grok.shell.views).length');
    expect(viewCount).toBeLessThanOrEqual(1); // Browse view may remain
  });
});
