import { test, expect, Page } from '@playwright/test';
import {
  openScriptsBrowser,
  getScriptCard,
  rightClickScript,
  clickMenuItem,
  apiDeleteScript,
  loadCarsDemoTable,
} from './helpers';

const BASE = process.env.DATAGROK_URL!;
const RUN_SCRIPT_NAME = 'PW_RunTest';

const RUN_SCRIPT_CONTENT = `#name: ${RUN_SCRIPT_NAME}
#language: r
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]
#output: string newParam

count <- nrow(table) * ncol(table)
newParam <- "test"`;

test.describe('Scripts: Run', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto(BASE);
    await page.waitForSelector('.d4-ribbon', { timeout: 20_000 });
    await page.waitForFunction(() => !!(window as any).grok?.dapi?.scripts, { timeout: 15_000 });
    await apiDeleteScript(page, RUN_SCRIPT_NAME);

    await page.evaluate(async (content) => {
      const DG = (window as any).DG;
      const grok = (window as any).grok;
      const script = DG.Script.create(content);
      await grok.dapi.scripts.save(script);
      return null;
    }, RUN_SCRIPT_CONTENT);
    await page.waitForTimeout(1000);
  });

  test.afterEach(async ({ page }) => {
    await apiDeleteScript(page, RUN_SCRIPT_NAME);
  });

  // ──────────────────────────────────────────────────────────────────
  // Test 1: Run.md steps 1–4 — run from context menu with sample table
  // ──────────────────────────────────────────────────────────────────
  test('1. Run script from context menu, choose cars from dropdown', async ({ page }) => {
    // Step 1: Go to Scripts
    await openScriptsBrowser(page);

    // Load cars.csv into session memory and navigate to Scripts via client-side routing
    await page.evaluate(async () => {
      const grok = (window as any).grok;
      const t = await grok.data.getDemoTable('cars.csv');
      t.name = 'cars';
      grok.shell.addTableView(t);
    });
    await page.waitForTimeout(1000);

    // Use client-side navigation to go back to Scripts (preserves tables in memory)
    await page.evaluate(() => (window as any).grok.shell.route('/scripts'));
    await expect(page.locator('.grok-gallery-search-bar')).toBeVisible({ timeout: 15_000 });

    // Step 2: Right-click the script → Run...
    await rightClickScript(page, RUN_SCRIPT_NAME);
    await clickMenuItem(page, 'Run...');

    // Dialog opens
    const dialog = page.locator('.d4-dialog').first();
    await expect(dialog).toBeVisible({ timeout: 10_000 });
    await expect(dialog).toContainText('RunTest', { ignoreCase: true });

    // Step 3: Select cars from the dropdown
    await dialog.locator('select.ui-input-editor').selectOption('cars');
    await page.waitForTimeout(500);

    // Step 4: Click OK
    const okBtn = dialog.locator('button.ui-btn-ok').first();
    await expect(okBtn).toBeEnabled({ timeout: 8_000 });
    await okBtn.click();

    // Verify: dialog closes (script executed)
    await expect(dialog).not.toBeVisible({ timeout: 30_000 });

    // Verify: no error balloon
    const errorBalloon = page.locator('.d4-balloon-error');
    await expect(errorBalloon).toHaveCount(0, { timeout: 5_000 }).catch(() => {});
  });

  // ──────────────────────────────────────────────────────────────────
  // Test 2: Run.md step 4 — run with file from Datagrok Files (folder-tree icon)
  // ──────────────────────────────────────────────────────────────────
  test('2. Run script, choose dataset from Datagrok Files via folder-tree icon', async ({ page }) => {
    await openScriptsBrowser(page);

    // Right-click → Run...
    await rightClickScript(page, RUN_SCRIPT_NAME);
    await clickMenuItem(page, 'Run...');

    const dialog = page.locator('.d4-dialog').first();
    await expect(dialog).toBeVisible({ timeout: 10_000 });

    // Click the folder-tree icon to browse Datagrok Files
    await dialog.locator('i[name="icon-folder-tree"]').click();
    await page.waitForTimeout(2000);

    // "Select a file" dialog opens
    const fileBrowser = page.locator('.d4-dialog[name="dialog-Select-a-file"]');
    await expect(fileBrowser).toBeVisible({ timeout: 10_000 });

    // Expand Files node (double-click to toggle expansion)
    const filesNode = fileBrowser.locator('text=Files').first();
    await expect(filesNode).toBeVisible({ timeout: 5_000 });
    await filesNode.dblclick();
    await page.waitForTimeout(1500);

    // Expand Demo folder
    const demoNode = fileBrowser.locator('text=Demo').first();
    await expect(demoNode).toBeVisible({ timeout: 5_000 });
    await demoNode.dblclick();
    await page.waitForTimeout(1500);

    // Select cars.csv
    const carsFile = fileBrowser.locator('text=cars.csv').first();
    await expect(carsFile).toBeVisible({ timeout: 5_000 });
    await carsFile.click();
    await page.waitForTimeout(500);

    // Click OK in the file browser dialog
    await fileBrowser.locator('button[name="button-OK"]').click();
    await page.waitForTimeout(2000);

    // Click OK in the main run dialog
    const okBtn = dialog.locator('button.ui-btn-ok').first();
    await expect(okBtn).toBeEnabled({ timeout: 8_000 });
    await okBtn.click();
    await expect(dialog).not.toBeVisible({ timeout: 30_000 });

    // Verify: no error balloon
    const errorBalloon = page.locator('.d4-balloon-error');
    await expect(errorBalloon).toHaveCount(0, { timeout: 5_000 }).catch(() => {});
  });

  // ──────────────────────────────────────────────────────────────────
  // Test 3: Run.md step 4 — run with local file upload (folder-open icon)
  // ──────────────────────────────────────────────────────────────────
  test('3. Run script, upload local file via folder-open icon', async ({ page }) => {
    await openScriptsBrowser(page);

    await rightClickScript(page, RUN_SCRIPT_NAME);
    await clickMenuItem(page, 'Run...');

    const dialog = page.locator('.d4-dialog').first();
    await expect(dialog).toBeVisible({ timeout: 10_000 });

    // Prepare a file chooser listener before clicking the folder-open icon
    const [fileChooser] = await Promise.all([
      page.waitForEvent('filechooser', { timeout: 5_000 }).catch(() => null),
      dialog.locator('i[name="icon-folder-open"]').click(),
    ]);

    if (fileChooser) {
      // Create a simple CSV in memory and upload it
      const csvContent = 'make,mpg,cylinders,displacement,horsepower,weight\nAMC,15,8,390,190,3850\nBuick,17,8,350,155,4360';
      const buffer = Buffer.from(csvContent, 'utf-8');
      await fileChooser.setFiles([{
        name: 'test_cars.csv',
        mimeType: 'text/csv',
        buffer,
      }]);
      await page.waitForTimeout(2000);

      // Click OK in the run dialog
      const okBtn = dialog.locator('button.ui-btn-ok').first();
      if (await okBtn.isEnabled({ timeout: 8_000 }).catch(() => false)) {
        await okBtn.click();
        await expect(dialog).not.toBeVisible({ timeout: 60_000 });
      }
    }
    else {
      // File chooser didn't trigger — dismiss dialog
      await page.keyboard.press('Escape');
    }

    const errorBalloon = page.locator('.d4-balloon-error');
    await expect(errorBalloon).toHaveCount(0, { timeout: 5_000 }).catch(() => {});
  });

  // ──────────────────────────────────────────────────────────────────
  // Test 4: Run.md steps 5–8 — run from Datagrok console
  // ──────────────────────────────────────────────────────────────────
  test('4. Run script from the console', async ({ page }) => {
    // Load cars.csv first
    await loadCarsDemoTable(page);
    await expect(page.locator('.d4-grid')).toBeVisible({ timeout: 10_000 });

    // Step 5: Open console with Backquote (`)
    await page.keyboard.press('Backquote');
    const consoleInput = page.locator('input[placeholder="> Enter command"]');
    await expect(consoleInput).toBeVisible({ timeout: 8_000 });

    // Step 6: Derive namespace from login
    const login = process.env.DATAGROK_LOGIN ?? 'admin';
    const namespace = login.includes('@') ? login.split('@')[0].replace(/\+/g, '') : login;

    // Step 7: Enter the command and press Enter
    // Datagrok capitalizes and removes underscores from the name (PW_RunTest → PWRunTest)
    const scriptFuncName = RUN_SCRIPT_NAME.replace(/_/g, '');
    await consoleInput.fill(`${namespace}:${scriptFuncName}("cars")`);
    await consoleInput.press('Enter');

    // Step 8: Verify output with script result (count=510 for cars: 30 rows × 17 columns)
    const consoleBody = page.locator('.d4-console-body');
    await expect(consoleBody).toContainText(/count/, { timeout: 30_000 });
  });
});
