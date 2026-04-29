import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function evalJs(page: Page, script: string): Promise<any> {
  return page.evaluate(script);
}

async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

test('Projects / Deleting', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const projectNames = [
    'AutoTest-Delete-1-' + Date.now(),
    'AutoTest-Delete-2-' + Date.now(),
  ];

  await loginToDatagrok(page);

  try {
    await softStep('Setup: Create test projects', async () => {
      for (const name of projectNames) {
        await evalJs(page, `(async () => {
          grok.shell.addTableView(grok.data.demo.demog());
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
    });

    await softStep('Case 1: Find test projects', async () => {
      for (const name of projectNames) {
        const exists = await evalJs(page, `(async () => {
          const p = await grok.dapi.projects.filter('name = "${name}"').first();
          return p !== null;
        })()`);
        expect(exists).toBe(true);
      }
    });

    await softStep('Case 2-3: Delete projects via API', async () => {
      for (const name of projectNames) {
        await evalJs(page, `(async () => {
          const p = await grok.dapi.projects.filter('name = "${name}"').first();
          if (p) await grok.dapi.projects.delete(p);
        })()`);
      }
      await page.waitForTimeout(2000);
    });

    await softStep('Case 4: Verify projects are deleted', async () => {
      for (const name of projectNames) {
        const exists = await evalJs(page, `(async () => {
          const p = await grok.dapi.projects.filter('name = "${name}"').first();
          return p !== null;
        })()`);
        expect(exists).toBe(false);
      }
    });
  } finally {
    for (const name of projectNames) {
      await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${name}"').first();
        if (p) await grok.dapi.projects.delete(p);
      })()`).catch(() => {});
    }
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
