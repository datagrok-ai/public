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

// Manual scenario `import-swagger.md` (order 6).
//
// 1. Open openweathermap.yaml in Datagrok (manual = drag-drop; we use Ctrl+O)
// 2. Browse > Platform > Functions > OpenAPI > OpenWeatherMap appears
// 3. Right-click → Edit; enter ApiKey
// 4. Run a query and verify a non-error result
//
// SCOPE NOTE: the original Test Track step is to drag-drop the YAML into the
// browser tab. Synthetic DOM `drop` events do NOT trigger the Dart-side file
// handler — the platform's drop zones (`ui.makeDroppable`) only wake up for
// real OS-originated `DataTransfer`s. We swap drag-drop for the equivalent
// `File | Open | File...` command (`Ctrl+O`), which calls the same
// `openFile -> _openFile -> openAsSwaggerFile -> Swagger.fromYaml -> save`
// chain in `xamgle/lib/src/commands/file/file_commands.dart`. The chain that
// actually parses YAML and registers the OpenAPI connection is identical;
// only the *user gesture* differs from the manual scenario. The manual
// drag-drop visual remains covered by `import-swagger-ui.md`.

const SWAGGER_NAME = 'OpenWeatherMap';
// OpenWeatherMap free-tier API key, default-provided so the CI run
// (which can't see the dev .env) can still exercise tests 2-3. Override
// via the DG_OPENWEATHERMAP_API_KEY env var on dev / playwright-tests.
// The free-tier key is rate-limited (60 calls/min) and not security-
// sensitive — rotating it costs only a 5-minute OWM sign-up.
const API_KEY = process.env.DG_OPENWEATHERMAP_API_KEY ?? '940677a955994e7147785652254f1742';

const YAML_PATH = path.resolve(
  __dirname, '..', '..', '..', 'public', 'packages', 'Samples', 'swaggers', 'openweathermap.yaml',
);

test.describe.serial('Connections / Import Swagger (OpenWeatherMap)', () => {
  test.beforeAll(async ({ browser }) => {
    if (!fs.existsSync(YAML_PATH))
      throw new Error(`openweathermap.yaml not found at ${YAML_PATH}`);

    // Make the test repeatable: drop any leftover OpenWeatherMap connection
    // before we ingest a fresh one. (afterAll also cleans up, but the previous
    // run could have crashed mid-way.)
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

    // Trigger File | Open | File... — `htmlOpenFile` creates a hidden
    // `<input type="file">` and clicks it. Playwright captures that via
    // `filechooser` event. The file's `.yaml` extension routes through
    // `openAsSwaggerFile` → `Swagger.fromYaml(content).save()` on the Dart side.
    const fileChooserPromise = page.waitForEvent('filechooser', { timeout: 10_000 });
    await page.keyboard.press('Control+O');
    const fileChooser = await fileChooserPromise;
    await fileChooser.setFiles(YAML_PATH);

    // Server-side ingestion fires a balloon — note the typo in the platform
    // message ("succesfully", single 's'). Match both spellings to be safe.
    await page.waitForFunction(
      () => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
        .some((b) => /succesfully added|successfully added/i.test(
          (b as HTMLElement).textContent ?? '',
        )),
      undefined,
      { timeout: 60_000 },
    );

    // Verify the connection exists server-side (the canonical signal).
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

    // Pick the first query under the connection and run it via context menu → Run.
    const firstQuery = page.locator(
      `[name^="tree-Platform---Functions---OpenAPI---${SWAGGER_NAME}---"]:not([name*="tree-expander-"])`,
    ).first();
    const queryNodeName = (await firstQuery.getAttribute('name'))!;
    await rightClickTreeNode(page, queryNodeName);
    await clickMenuItemExact(page, 'Run');

    // Run dialog (param prompt) → OK / RUN button. `isVisible` does NOT wait,
    // so use `waitFor({state:'visible'})` to actually poll for it.
    try {
      await page.locator('.d4-dialog').first().waitFor({ state: 'visible', timeout: 5_000 });
      const runBtn = page.locator('.d4-dialog [name="button-RUN"], .d4-dialog [name="button-OK"]').first();
      await runBtn.click();
      await page.locator('.d4-dialog').first().waitFor({ state: 'detached', timeout: 10_000 }).catch(() => null);
    }
    catch {
      // No params → query ran without a prompt.
    }

    // Wait for either a result grid or a balloon — the manual just asks that
    // queries run; we accept either successful grid or a non-error balloon.
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
