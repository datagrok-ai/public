import { test, expect, Page, BrowserContext } from '@playwright/test';
import * as path from 'path';
import {
  openScriptsBrowser,
  getScriptCard,
  rightClickScript,
  clickMenuItem,
  apiDeleteScript,
  loadCarsDemoTable,
} from './helpers';

const BASE = process.env.DATAGROK_URL!;
const AUTH_STATE = path.resolve(__dirname, '..', '.auth.json');
const RUN_SCRIPT_NAME = 'PW_RunTest';

const RUN_SCRIPT_CONTENT = `#name: ${RUN_SCRIPT_NAME}
#language: r
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]
#output: string newParam

count <- nrow(table) * ncol(table)
newParam <- "test"`;

async function resetToScripts(page: Page) {
  await page.evaluate(() => {
    const g = (window as any).grok;
    if (g?.shell?.closeAll) g.shell.closeAll();
    document.querySelectorAll('.d4-dialog').forEach((d: any) => { try { d.remove(); } catch(_) {} });
    document.querySelectorAll('.d4-toast, .d4-balloon, .d4-menu').forEach(e => e.remove());
  });
  await page.waitForTimeout(300);
  const scriptsLabel = page.locator('.d4-tree-view-item-label', { hasText: /^Scripts$/i }).first();
  if (await scriptsLabel.isVisible({ timeout: 1_000 }).catch(() => false)) {
    await scriptsLabel.click();
    await expect(page.locator('.grok-gallery-search-bar')).toBeVisible({ timeout: 10_000 });
  } else {
    await openScriptsBrowser(page);
  }
}

test.describe.serial('Scripts: Run', () => {
  let sharedContext: BrowserContext;
  let page: Page;

  test.beforeAll(async ({ browser }) => {

    test.setTimeout(180_000);
    sharedContext = await browser.newContext({ storageState: AUTH_STATE });
    page = await sharedContext.newPage();
    await openScriptsBrowser(page);
  });

  test.afterAll(async () => {
    await sharedContext?.close();
  });

  test.beforeEach(async () => {
    await page.evaluate(() => {
      const g = (window as any).grok;
      if (g?.shell?.closeAll) g.shell.closeAll();
    });
    await page.waitForTimeout(300);
    await page.waitForFunction(() => !!(window as any).grok?.dapi?.scripts, { timeout: 10_000 });
    await apiDeleteScript(page, RUN_SCRIPT_NAME);

    await page.evaluate(async (content) => {
      const DG = (window as any).DG;
      const grok = (window as any).grok;
      const script = DG.Script.create(content);
      await grok.dapi.scripts.save(script);
      return null;
    }, RUN_SCRIPT_CONTENT);
    await page.waitForTimeout(300);
  });

  test.afterEach(async () => {
    await apiDeleteScript(page, RUN_SCRIPT_NAME);
  });

  test('1. Run script from context menu, choose cars from dropdown', async () => {

    test.setTimeout(180_000);

    await resetToScripts(page);

    await page.evaluate(async () => {
      const grok = (window as any).grok;
      const t = await grok.data.getDemoTable('cars.csv');
      t.name = 'cars';
      grok.shell.addTableView(t);
    });
    await page.waitForTimeout(300);

    await page.evaluate(() => (window as any).grok.shell.route('/scripts'));
    await expect(page.locator('.grok-gallery-search-bar')).toBeVisible({ timeout: 15_000 });

    await rightClickScript(page, RUN_SCRIPT_NAME);
    await clickMenuItem(page, 'Run...');

    const dialog = page.locator('.d4-dialog').first();
    await expect(dialog).toBeVisible({ timeout: 10_000 });
    await expect(dialog).toContainText('RunTest', { ignoreCase: true });

    await dialog.locator('select.ui-input-editor').selectOption('cars');

    const okBtn = dialog.locator('button.ui-btn-ok').first();
    await expect(okBtn).toBeEnabled({ timeout: 8_000 });
    await okBtn.click();

    await expect(dialog).not.toBeVisible({ timeout: 120_000 });

    const errorBalloon = page.locator('.d4-balloon-error');
    await expect(errorBalloon).toHaveCount(0, { timeout: 5_000 }).catch(() => {});
  });

  test('2. Run script, choose dataset from Datagrok Files via folder-tree icon', async () => {
    test.setTimeout(180_000);
    await resetToScripts(page);

    await rightClickScript(page, RUN_SCRIPT_NAME);
    await clickMenuItem(page, 'Run...');

    const dialog = page.locator('.d4-dialog').first();
    await expect(dialog).toBeVisible({ timeout: 10_000 });

    await dialog.locator('i[name="icon-folder-tree"]').click();

    const fileBrowser = page.locator('.d4-dialog', { has: page.locator('.dlg-select-a-file') }).first();
    await expect(fileBrowser).toBeVisible({ timeout: 10_000 });

    const expandTreeGroup = async (name: string) => {
      const label = fileBrowser
        .locator(`.d4-tree-view-group-label:text-is("${name}")`)
        .first();

      await label.evaluate((el: Element) => {
        el.scrollIntoView({ block: 'center' });
      });
      await expect(label).toBeVisible({ timeout: 10_000 });
      const handle = await label.elementHandle();
      const info = await handle!.evaluate((el: Element) => {
        const node = el.parentElement!;
        const group = node.parentElement!;
        const tri = node.querySelector(':scope > .d4-tree-view-tri') as HTMLElement | null;
        return {
          expanded: !!tri?.classList.contains('d4-tree-view-tri-expanded'),
          hasTri: !!tri,
          groupClasses: group.className,
        };
      });
      if (!info.expanded && info.hasTri) {

        const node = label.locator('xpath=..');
        await node.locator('.d4-tree-view-tri').first().click();

        await fileBrowser.page().waitForFunction((lbl) => {
          const labelEl = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
            .find(l => (l.textContent ?? '').trim() === lbl);
          if (!labelEl) return false;
          const nodeEl = labelEl.parentElement!;
          const tri = nodeEl.querySelector(':scope > .d4-tree-view-tri');
          const group = nodeEl.parentElement!;
          const host = group.querySelector(':scope > .d4-tree-view-group-host');
          return !!tri?.classList.contains('d4-tree-view-tri-expanded')
            && !!host && host.children.length > 0;
        }, name, { timeout: 15_000 }).catch(() => {});
        await fileBrowser.page().waitForTimeout(500);
      }
    };

    await expandTreeGroup('Files');
    await expandTreeGroup('Demo');

    const carsFile = fileBrowser.locator('.d4-tree-view-item-label:text-is("cars.csv")').first();
    await expect(carsFile).toBeVisible({ timeout: 10_000 });
    await carsFile.click();

    await fileBrowser.locator('button.ui-btn-ok, button[name="button-OK"]').first().click();
    await page.waitForTimeout(500);

    const okBtn = dialog.locator('button.ui-btn-ok').first();
    await expect(okBtn).toBeEnabled({ timeout: 8_000 });
    await okBtn.click();
    await expect(dialog).not.toBeVisible({ timeout: 60_000 });

    const errorBalloon = page.locator('.d4-balloon-error');
    await expect(errorBalloon).toHaveCount(0, { timeout: 5_000 }).catch(() => {});
  });

  test('3. Run script, upload local file via folder-open icon', async () => {
    await resetToScripts(page);

    await rightClickScript(page, RUN_SCRIPT_NAME);
    await clickMenuItem(page, 'Run...');

    const dialog = page.locator('.d4-dialog').first();
    await expect(dialog).toBeVisible({ timeout: 10_000 });

    const [fileChooser] = await Promise.all([
      page.waitForEvent('filechooser', { timeout: 5_000 }).catch(() => null),
      dialog.locator('i[name="icon-folder-open"]').click(),
    ]);

    if (fileChooser) {
      const csvContent = 'make,mpg,cylinders,displacement,horsepower,weight\nAMC,15,8,390,190,3850\nBuick,17,8,350,155,4360';
      const buffer = Buffer.from(csvContent, 'utf-8');
      await fileChooser.setFiles([{
        name: 'test_cars.csv',
        mimeType: 'text/csv',
        buffer,
      }]);
      await page.waitForTimeout(500);

      const okBtn = dialog.locator('button.ui-btn-ok').first();
      if (await okBtn.isEnabled({ timeout: 8_000 }).catch(() => false)) {
        await okBtn.click();
        await expect(dialog).not.toBeVisible({ timeout: 60_000 });
      }
    }
    else {
      await page.keyboard.press('Escape');
    }

    const errorBalloon = page.locator('.d4-balloon-error');
    await expect(errorBalloon).toHaveCount(0, { timeout: 5_000 }).catch(() => {});
  });

  test('4. Run script from the console', async () => {

    await loadCarsDemoTable(page);
    await expect(page.locator('.d4-grid')).toBeVisible({ timeout: 10_000 });

    await page.keyboard.press('Backquote');
    const consoleInput = page.locator('input[placeholder="> Enter command"]');
    await expect(consoleInput).toBeVisible({ timeout: 8_000 });

    const login = process.env.DATAGROK_LOGIN ?? 'admin';
    const namespace = login.includes('@') ? login.split('@')[0].replace(/\+/g, '') : login;

    const scriptFuncName = RUN_SCRIPT_NAME.replace(/_/g, '');
    await consoleInput.fill(`${namespace}:${scriptFuncName}("cars")`);
    await consoleInput.press('Enter');

    const consoleBody = page.locator('.d4-console-body');
    await expect(consoleBody).toContainText(/count/, { timeout: 30_000 });
  });
});
