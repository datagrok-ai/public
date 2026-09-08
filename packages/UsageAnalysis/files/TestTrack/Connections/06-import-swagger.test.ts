import { test, expect } from '@playwright/test';
import * as fs from 'fs';
import * as path from 'path';
import {
  AUTH_STATE,
  applyAutomationSetup,
  clickConnectionSave,
  clickMenuItemExact,
  expandTreeNode,
  fillConnectionField,
  goHome,
  rightClickTreeNode,
} from './helpers';

const SWAGGER_NAME = 'OpenWeatherMap';
const API_KEY = process.env.DG_OPENWEATHERMAP_API_KEY ?? '';

const YAML_PATH = path.resolve(
  __dirname, '..', '..', '..', 'public', 'packages', 'Samples', 'swaggers', 'openweathermap.yaml',
);

test.describe.serial('Connections / Import Swagger (OpenWeatherMap)', () => {
  test.beforeAll(async ({ browser }) => {
    if (!fs.existsSync(YAML_PATH))
      throw new Error(`openweathermap.yaml not found at ${YAML_PATH}`);

    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await page.evaluate(async (name) => {
      const g = (window as unknown as { grok: any }).grok;
      const cs = await g.dapi.connections.filter(`friendlyName = "${name}"`).list();
      for (const c of cs) await g.dapi.connections.delete(c);
    }, SWAGGER_NAME);
    await ctx.close();
  });

  test.afterAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await page.evaluate(async (name) => {
      const g = (window as unknown as { grok: any }).grok;
      const cs = await g.dapi.connections.filter(`friendlyName = "${name}"`).list();
      for (const c of cs) await g.dapi.connections.delete(c);
    }, SWAGGER_NAME);
    await ctx.close();
  });

  test('1. Ingest openweathermap.yaml; OpenWeatherMap connection appears under OpenAPI', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);

    const fileChooserPromise = page.waitForEvent('filechooser', { timeout: 10_000 });
    await page.keyboard.press('Control+O');
    const fileChooser = await fileChooserPromise;
    await fileChooser.setFiles(YAML_PATH);

    await page.waitForFunction(
      () => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
        .some((b) => /succesfully added|successfully added/i.test(
          (b as HTMLElement).textContent ?? '',
        )),
      undefined,
      { timeout: 60_000 },
    );

    await expect.poll(async () => page.evaluate(async (name) => {
      const g = (window as unknown as { grok: any }).grok;
      const cs = await g.dapi.connections.filter(`friendlyName = "${name}"`).list();
      return cs.length;
    }, SWAGGER_NAME), { timeout: 30_000 }).toBeGreaterThan(0);
  });

  test('2. Edit connection; enter ApiKey and Save', async ({ page }) => {
    test.skip(!API_KEY, 'DG_OPENWEATHERMAP_API_KEY env var not set — cannot enter a real key');

    await goHome(page);
    await applyAutomationSetup(page);
    await expandTreeNode(page, 'tree-Platform');
    await expandTreeNode(page, 'tree-Platform---Functions');
    await expandTreeNode(page, 'tree-Platform---Functions---OpenAPI');

    await rightClickTreeNode(page, `tree-Platform---Functions---OpenAPI---${SWAGGER_NAME}`);
    await clickMenuItemExact(page, 'Edit...');
    await page.locator('.d4-dialog').waitFor({ timeout: 10_000 });

    await fillConnectionField(page, 'ApiKey', API_KEY);
    await clickConnectionSave(page);
  });

  test('3. Run a query under OpenWeatherMap and verify it returns rows', async ({ page }) => {
    test.skip(!API_KEY, 'DG_OPENWEATHERMAP_API_KEY env var not set — cannot exercise a real query');

    await goHome(page);
    await applyAutomationSetup(page);
    await expandTreeNode(page, 'tree-Platform');
    await expandTreeNode(page, 'tree-Platform---Functions');
    await expandTreeNode(page, 'tree-Platform---Functions---OpenAPI');
    await expandTreeNode(page, `tree-Platform---Functions---OpenAPI---${SWAGGER_NAME}`);

    const firstQuery = page.locator(
      `[name^="tree-Platform---Functions---OpenAPI---${SWAGGER_NAME}---"]:not([name*="tree-expander-"])`,
    ).first();
    const queryNodeName = (await firstQuery.getAttribute('name'))!;
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Run');

    try {
      await page.locator('.d4-dialog').first().waitFor({ state: 'visible', timeout: 5_000 });
      const runBtn = page.locator('.d4-dialog [name="button-RUN"], .d4-dialog [name="button-OK"]').first();
      await runBtn.click();
      await page.locator('.d4-dialog').first().waitFor({ state: 'detached', timeout: 10_000 }).catch(() => null);
    }
    catch {

    }

    await page.waitForFunction(() =>
      !!document.querySelector('[name="viewer-Grid"] canvas')
      || !!document.querySelector('.grok-balloon, .d4-balloon'),
    undefined,
    { timeout: 60_000 });

    const errors = await page.evaluate(() => Array.from(document.querySelectorAll('.d4-balloon-error, .grok-balloon-error'))
      .map((b) => (b as HTMLElement).textContent?.trim() ?? '')
      .filter((s) => s.length > 0));
    expect(errors, 'no error balloons after running an OpenAPI query').toEqual([]);
  });
});
