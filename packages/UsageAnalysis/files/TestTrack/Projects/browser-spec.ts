import {test, expect, Page} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function evalJs(page: Page, script: string): Promise<any> {
  return page.evaluate(script);
}

async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

test('Projects / Browser', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await closeAll(page);

  await softStep('Case 1: Navigate to Browse > Dashboards', async () => {
    await page.goto(`${baseUrl}/projects`);
    await page.waitForSelector('.grok-gallery-grid', {timeout: 30000});
    await expect(page.locator('.grok-gallery-grid')).toBeVisible();
  });

  await softStep('Case 2: Search for projects in Dashboards', async () => {
    await page.goto(`${baseUrl}/projects`);
    await page.waitForSelector('.grok-gallery-grid', {timeout: 30000});

    const searchInput = page.locator('input[placeholder*="Search"]');
    await searchInput.fill('Test Project');
    await page.waitForTimeout(2000);

    const cards = page.locator('.grok-gallery-grid .d4-item-card');
    const count = await cards.count();
    expect(count).toBeGreaterThan(0);
  });

  await softStep('Case 5: Check Context Panel for selected project', async () => {
    await page.goto(`${baseUrl}/projects`);
    await page.waitForSelector('.grok-gallery-grid', {timeout: 30000});

    const firstCard = page.locator('.grok-gallery-grid .d4-item-card').first();
    await firstCard.click();
    await page.waitForTimeout(1000);

    const contextPanel = page.locator('.grok-prop-panel');
    await expect(contextPanel).toBeVisible();
    const panelText = await contextPanel.textContent();
    expect(panelText).toContain('Sharing');
  });

  await softStep('Case 8: Review Context Panel tabs', async () => {
    await page.goto(`${baseUrl}/projects`);
    await page.waitForSelector('.grok-gallery-grid', {timeout: 30000});

    const firstCard = page.locator('.grok-gallery-grid .d4-item-card').first();
    await firstCard.click();
    await page.waitForTimeout(1000);

    const contextPanel = page.locator('.grok-prop-panel');
    const panelText = await contextPanel.textContent() ?? '';
    for (const section of ['Details', 'Sharing', 'Activity'])
      expect(panelText).toContain(section);
  });

  await softStep('Case 9: Open a project from Dashboards', async () => {
    const projectName = 'AutoTest-Browser-' + Date.now();

    await evalJs(page, `(async () => {
      grok.shell.addTableView(grok.data.demo.demog());
    })()`);
    await page.waitForTimeout(2000);

    await evalJs(page, `(async () => {
      const project = grok.shell.project;
      project.name = '${projectName}';
      await grok.dapi.projects.save(project);
    })()`);
    await page.waitForTimeout(3000);

    await closeAll(page);

    const tables = await evalJs(page, `(async () => {
      const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
      if (!p) return [];
      await p.open();
      return grok.shell.tables.map(t => t.name);
    })()`);
    expect(tables.length).toBeGreaterThan(0);

    await evalJs(page, `(async () => {
      const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
      if (p) await grok.dapi.projects.delete(p);
    })()`);
  });

  await closeAll(page);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
