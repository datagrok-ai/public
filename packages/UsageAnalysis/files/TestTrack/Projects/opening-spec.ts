import {test, expect, Page} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function evalJs(page: Page, script: string): Promise<any> {
  return page.evaluate(script);
}

async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

test('Projects / Opening', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const projectName = 'AutoTest-Opening-' + Date.now();

  await loginToDatagrok(page);

  try {
    await softStep('Setup: Create test project with demo data', async () => {
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
    });

    await softStep('Case 1: Navigate to Browse > Dashboards', async () => {
      await page.goto(`${baseUrl}/projects`);
      await page.waitForSelector('.grok-gallery-grid', {timeout: 30000});
      await expect(page.locator('.grok-gallery-grid')).toBeVisible();
    });

    await softStep('Case 2: Find project in Dashboards', async () => {
      await page.goto(`${baseUrl}/projects`);
      await page.waitForSelector('.grok-gallery-grid', {timeout: 30000});

      const found = await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
        return p !== null;
      })()`);
      expect(found).toBe(true);
    });

    await softStep('Case 3: Check Context Panel attributes', async () => {
      await page.goto(`${baseUrl}/projects`);
      await page.waitForSelector('.grok-gallery-grid', {timeout: 30000});

      const attrs = await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
        return { name: p.name, createdOn: p.createdOn !== null, author: p.author !== null };
      })()`);
      expect(attrs.name).toBe(projectName);
      expect(attrs.createdOn).toBe(true);
      expect(attrs.author).toBe(true);
    });

    await softStep('Case 4: Open project and verify data', async () => {
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

    await softStep('Case 5: Close all views', async () => {
      await closeAll(page);
      const viewCount = await evalJs(page, 'Array.from(grok.shell.views).length');
      expect(viewCount).toBeLessThanOrEqual(1);
    });
  } finally {
    await evalJs(page, `(async () => {
      const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
      if (p) await grok.dapi.projects.delete(p);
    })()`).catch(() => {});
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
