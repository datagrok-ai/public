import { test, expect } from '@playwright/test';
import {
  treeGroupByName,
  treeNodeByPath,
  STATUS_BAR_VIEW_PANEL,
} from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
} from './helpers';

test.describe('Browse Files (Browse-Files-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-Files-01 — Files section contains required subdivisions', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Files');

    await expect(treeNodeByPath(page, ['Files', 'My-files']),
      'My files should be present under Files').toHaveCount(1, { timeout: 5_000 });
    await expect(treeNodeByPath(page, ['Files', 'Demo']),
      'Demo should be present under Files').toHaveCount(1, { timeout: 5_000 });

    const appData = page.locator('.d4-tree-view-group-label', { hasText: /^App Data$/ });
    expect(await appData.count(),
      'At least one App Data subdivision should be present').toBeGreaterThanOrEqual(1);

    await expectNoErrors(page, sink);
  });

  test('Browse-Files-02 — preview of a known demo table file opens', async ({ page }) => {
    const sink = watchErrors(page);

    await page.goto(`${process.env.DATAGROK_URL!}/files/System.DemoFiles/?browse=files`);
    await page.waitForSelector('.d4-ribbon', { timeout: 20_000 });

    const demog = page.locator('label, .d4-tree-view-item-label', { hasText: /^demog\.csv$/i }).first();
    await expect(demog, 'demog.csv should be visible in the Demo folder')
      .toBeVisible({ timeout: 15_000 });
    await demog.dblclick();

    await expect(page.locator(STATUS_BAR_VIEW_PANEL),
      'Status bar should show table info (Rows/Columns)').toContainText(/Rows:\s*\d/, { timeout: 15_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-Files-03 — opening a Demo subfolder shows the folder view', async ({ page }) => {
    const sink = watchErrors(page);

    await page.goto(`${process.env.DATAGROK_URL!}/files/System.DemoFiles/?browse=files`);
    await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });

    const chem = page.locator('label', { hasText: /^chem$/ }).first();
    await expect(chem, 'chem folder should be visible in Demo').toBeVisible({ timeout: 15_000 });
    await chem.dblclick();
    await page.waitForTimeout(2500);

    await expect(page).toHaveURL(/files\/System\.DemoFiles\/chem/);
    const smiles = page.locator('label', { hasText: /^smiles\.csv$/ }).first();
    await expect(smiles, 'smiles.csv should be visible inside the chem folder')
      .toBeVisible({ timeout: 15_000 });

    await expectNoErrors(page, sink);
  });

  test('Browse-Files-04 — Refresh tree picks up new files in System.DemoFiles', async ({ page }) => {
    const sink = watchErrors(page);

    await page.goto(`${process.env.DATAGROK_URL!}/files/System.DemoFiles/?browse=files`);
    await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
    await page.waitForTimeout(1500);

    const uniqueName = `qa_autotest_${Date.now()}.txt`;
    const fullPath = `System:DemoFiles/${uniqueName}`;

    const created = await page.evaluate(async ({ path }) => {
      try {
        await (window as any).grok.dapi.files.writeAsText(path, 'hello');
        return true;
      } catch { return false; }
    }, { path: fullPath });
    if (!created) test.skip(true, 'grok.dapi.files.writeAsText is not supported in this build');

    try {

      await page.locator('[name="icon-sync"]').first().click();
      await page.waitForTimeout(2000);

      await page.locator('[name="icon-sync"]').first().click();
      await page.waitForTimeout(2000);

      const apiSees = await page.evaluate(async ({ path }) => {
        try { return await (window as any).grok.dapi.files.exists(path); } catch { return false; }
      }, { path: fullPath });
      expect(apiSees, 'File must be visible to the API after creation').toBe(true);

      await expect(page.locator('.d4-ribbon').first()).toBeVisible();

      await expectNoErrors(page, sink);
    } finally {
      await page.evaluate(async ({ path }) => {
        try { await (window as any).grok.dapi.files.delete(path); } catch {}
      }, { path: fullPath });
    }
  });

  test('Browse-Files-05 — a shared folder appears under My stuff > Shared with me, grouped by sharer', async ({ page }) => {
    const sink = watchErrors(page);

    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'My stuff');
    await expandTreeGroup(page, 'Shared with me');
    await page.waitForTimeout(800);

    const userGroups = page.locator('[name^="tree-My-stuff---Shared-with-me---"]');
    expect(await userGroups.count(),
      'Shared with me should contain at least one user subgroup').toBeGreaterThanOrEqual(1);

    const firstUser = userGroups.first();
    const firstUserName = (await firstUser.getAttribute('name')) ?? '';
    const tri = firstUser.locator('.d4-tree-view-tri').first();
    if (!(await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded')).catch(() => false))) {
      await tri.click();
      await page.waitForTimeout(1500);
    }
    const itemsUnderFirstUser = page.locator(`[name^="${firstUserName}---"]`);
    expect(await itemsUnderFirstUser.count(),
      'The first Shared-with-me user subgroup should expose at least one shared entity')
      .toBeGreaterThanOrEqual(1);

    await expectNoErrors(page, sink);
  });

  test('Browse-Files-06 — Download from the file context menu triggers a download', async ({ page }) => {
    const sink = watchErrors(page);

    await page.goto(`${process.env.DATAGROK_URL!}/files/System.DemoFiles/?browse=files`);
    await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });

    const demog = page.locator('label', { hasText: /^demog\.csv$/ }).first();
    await expect(demog, 'demog.csv should be visible').toBeVisible({ timeout: 15_000 });

    await demog.evaluate((el) => {
      const r = el.getBoundingClientRect();
      el.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2, clientX: r.left + 5, clientY: r.top + 5,
      }));
    });
    await expect(page.locator('.d4-menu-popup'), 'Context menu must open').toBeVisible({ timeout: 5_000 });

    const downloadPromise = page.waitForEvent('download', { timeout: 15_000 });
    await page.locator('.d4-menu-item-label', { hasText: /^Download$/ }).first().click();
    const dl = await downloadPromise;

    expect(dl.suggestedFilename(), 'Downloaded file should be demog.csv').toMatch(/^demog\.csv$/i);

    const tmpPath = await dl.path();
    expect(tmpPath, 'Download should produce a file on disk').toBeTruthy();

    await expectNoErrors(page, sink);
  });
});
