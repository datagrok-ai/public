import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function evalJs(page: Page, script: string): Promise<any> {
  return page.evaluate(script);
}

async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

async function saveProject(page: Page, name: string) {
  await page.click('button:has-text("SAVE"), .ui-btn:has-text("SAVE")');
  const dialog = page.locator('.d4-dialog');
  await dialog.waitFor({timeout: 10000});
  const nameInput = dialog.locator('input[type="text"]');
  await nameInput.fill(name);
  await dialog.locator('.ui-btn.ui-btn-ok').click();
  await page.waitForSelector('text=Share ' + name, {timeout: 30000});
  await page.locator('.d4-dialog .ui-btn:has-text("CANCEL")').click();
}

async function deleteProject(page: Page, name: string) {
  await evalJs(page, `(async () => {
    const p = await grok.dapi.projects.filter('name = "${name}"').first();
    if (p) await grok.dapi.projects.delete(p);
  })()`);
}

test('Projects / Uploading', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await closeAll(page);

  await softStep('Case 1: Save project from two local tables', async () => {
    const projectName = 'AutoTest-Upload-Case1-' + Date.now();
    try {
      await evalJs(page, `(async () => {
        grok.shell.addTableView(grok.data.demo.demog());
        const df2 = await grok.data.getDemoTable('demog.csv');
        df2.name = 'demog2';
        grok.shell.addTableView(df2);
      })()`);

      await page.waitForTimeout(2000);
      const tables = await evalJs(page, 'grok.shell.tables.map(t => t.name)');
      expect(tables.length).toBe(2);

      await saveProject(page, projectName);

      const exists = await evalJs(page,
        `(async () => {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          return p !== null;
        })()`,
      );
      expect(exists).toBe(true);
    } finally {
      await deleteProject(page, projectName).catch(() => {});
      await closeAll(page);
    }
  });

  await softStep('Case 3: Save project from Browse > Files', async () => {
    const projectName = 'AutoTest-Upload-Case3-' + Date.now();
    try {
      await evalJs(page, `(async () => {
        const df1 = await grok.data.files.openTable('System:DemoFiles/demog.csv');
        grok.shell.addTableView(df1);
        const df2 = await grok.data.files.openTable('System:DemoFiles/geo/earthquakes.csv');
        grok.shell.addTableView(df2);
      })()`);

      await page.waitForTimeout(2000);
      const tables = await evalJs(page, 'grok.shell.tables.map(t => t.name)');
      expect(tables.length).toBe(2);

      await saveProject(page, projectName);

      const exists = await evalJs(page,
        `(async () => {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          return p !== null;
        })()`,
      );
      expect(exists).toBe(true);
    } finally {
      await deleteProject(page, projectName).catch(() => {});
      await closeAll(page);
    }
  });

  await softStep('Case 4: Save project from query results', async () => {
    const projectName = 'AutoTest-Upload-Case4-' + Date.now();
    try {
      await evalJs(page, `(async () => {
        const queries = await grok.dapi.queries.list({limit: 500});
        const q1 = queries.find(q => q.name === 'PostgresProducts');
        const q2 = queries.find(q => q.name === 'PostgresCustomers');
        const r1 = await q1.executeTable();
        grok.shell.addTableView(r1);
        const r2 = await q2.executeTable();
        grok.shell.addTableView(r2);
      })()`);

      await page.waitForTimeout(5000);
      const tables = await evalJs(page, 'grok.shell.tables.map(t => t.name)');
      expect(tables.length).toBe(2);

      await saveProject(page, projectName);

      const exists = await evalJs(page,
        `(async () => {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          return p !== null;
        })()`,
      );
      expect(exists).toBe(true);
    } finally {
      await deleteProject(page, projectName).catch(() => {});
      await closeAll(page);
    }
  });

  await softStep('Case 9: Save project with joined tables', async () => {
    const projectName = 'AutoTest-Upload-Case9-' + Date.now();
    try {
      await evalJs(page, `(async () => {
        const df1 = await grok.data.files.openTable('System:DemoFiles/demog.csv');
        grok.shell.addTableView(df1);
        const df2 = await grok.data.files.openTable('System:DemoFiles/demog.csv');
        df2.name = 'demog_copy';
        grok.shell.addTableView(df2);

        const innerJoin = grok.data.joinTables(df1, df2,
          ['USUBJID'], ['USUBJID'],
          ['AGE', 'SEX'], ['RACE', 'DIS_POP'], 'inner');
        grok.shell.addTableView(innerJoin);

        const leftJoin = grok.data.joinTables(df1, df2,
          ['USUBJID'], ['USUBJID'],
          ['HEIGHT', 'WEIGHT'], ['SEVERITY'], 'left');
        grok.shell.addTableView(leftJoin);
      })()`);

      await page.waitForTimeout(3000);
      const tables = await evalJs(page, 'grok.shell.tables.map(t => t.name)');
      expect(tables.length).toBe(4);

      await saveProject(page, projectName);

      const exists = await evalJs(page,
        `(async () => {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          return p !== null;
        })()`,
      );
      expect(exists).toBe(true);
    } finally {
      await deleteProject(page, projectName).catch(() => {});
      await closeAll(page);
    }
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
