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

test.describe('Projects / Deleting', () => {
  const projectNames = [
    'AutoTest-Delete-1-' + Date.now(),
    'AutoTest-Delete-2-' + Date.now(),
  ];

  test.beforeAll(async ({ browser }) => {
    const page = await browser.newPage();
    await login(page);

    for (const name of projectNames) {
      await evalJs(page, `(async () => {
        const tv = grok.shell.addTableView(grok.data.demo.demog());
      })()`);
      await page.waitForTimeout(2000);
      await evalJs(page, `(async () => {
        const project = grok.shell.project;
        project.name = '${name}';
        await grok.dapi.projects.save(project);
      })()`);
      await page.waitForTimeout(3000);
      await closeAll(page);
    }
    await page.close();
  });

  // Cleanup in case tests fail
  test.afterAll(async ({ browser }) => {
    const page = await browser.newPage();
    await login(page);
    for (const name of projectNames) {
      await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${name}"').first();
        if (p) await grok.dapi.projects.delete(p);
      })()`);
    }
    await page.close();
  });

  test('Case 1: Find test projects', async ({ page }) => {
    await login(page);

    for (const name of projectNames) {
      const exists = await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${name}"').first();
        return p !== null;
      })()`);
      expect(exists).toBe(true);
    }
  });

  test('Case 2-3: Delete projects via API', async ({ page }) => {
    await login(page);

    for (const name of projectNames) {
      await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${name}"').first();
        if (p) await grok.dapi.projects.delete(p);
      })()`);
    }
    await page.waitForTimeout(2000);
  });

  test('Case 4: Verify projects are deleted', async ({ page }) => {
    await login(page);

    for (const name of projectNames) {
      const exists = await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${name}"').first();
        return p !== null;
      })()`);
      expect(exists).toBe(false);
    }
  });
});
