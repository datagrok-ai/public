import { test, expect, Page } from '@playwright/test';

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

test.describe('Projects / Browser', () => {

  test.beforeEach(async ({ page }) => {
    await login(page);
    await closeAll(page);
  });

  test.afterEach(async ({ page }) => {
    await closeAll(page);
  });

  test('Case 1: Navigate to Browse > Dashboards', async ({ page }) => {
    await page.goto(`${BASE_URL}/projects`);
    await page.waitForSelector('.grok-gallery-grid', { timeout: 30000 });
    const gallery = page.locator('.grok-gallery-grid');
    await expect(gallery).toBeVisible();
  });

  test('Case 2: Search for projects in Dashboards', async ({ page }) => {
    await page.goto(`${BASE_URL}/projects`);
    await page.waitForSelector('.grok-gallery-grid', { timeout: 30000 });

    const searchInput = page.locator('input[placeholder*="Search"]');
    await searchInput.fill('Test Project');
    await page.waitForTimeout(2000);

    const cards = page.locator('.grok-gallery-grid .d4-item-card');
    const count = await cards.count();
    expect(count).toBeGreaterThan(0);
  });

  test('Case 5: Check Context Panel for selected project', async ({ page }) => {
    await page.goto(`${BASE_URL}/projects`);
    await page.waitForSelector('.grok-gallery-grid', { timeout: 30000 });

    // Click first project card
    const firstCard = page.locator('.grok-gallery-grid .d4-item-card').first();
    await firstCard.click();
    await page.waitForTimeout(1000);

    // Check Context Panel sections
    const contextPanel = page.locator('.grok-prop-panel');
    await expect(contextPanel).toBeVisible();
    const panelText = await contextPanel.textContent();
    expect(panelText).toContain('Sharing');
  });

  test('Case 8: Review Context Panel tabs', async ({ page }) => {
    await page.goto(`${BASE_URL}/projects`);
    await page.waitForSelector('.grok-gallery-grid', { timeout: 30000 });

    const firstCard = page.locator('.grok-gallery-grid .d4-item-card').first();
    await firstCard.click();
    await page.waitForTimeout(1000);

    const contextPanel = page.locator('.grok-prop-panel');
    const panelText = await contextPanel.textContent();
    for (const section of ['Details', 'Sharing', 'Activity']) {
      expect(panelText).toContain(section);
    }
  });

  test('Case 9: Open a project from Dashboards', async ({ page }) => {
    const projectName = 'AutoTest-Browser-' + Date.now();

    // Create a project to open
    await evalJs(page, `(async () => {
      const tv = grok.shell.addTableView(grok.data.demo.demog());
    })()`);
    await page.waitForTimeout(2000);

    // Save it
    await evalJs(page, `(async () => {
      const project = grok.shell.project;
      project.name = '${projectName}';
      await grok.dapi.projects.save(project);
    })()`);
    await page.waitForTimeout(3000);

    await closeAll(page);

    // Open via API
    const tables = await evalJs(page, `(async () => {
      const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
      if (!p) return [];
      await p.open();
      return grok.shell.tables.map(t => t.name);
    })()`);
    expect(tables.length).toBeGreaterThan(0);

    // Cleanup
    await evalJs(page, `(async () => {
      const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
      if (p) await grok.dapi.projects.delete(p);
    })()`);
  });
});
