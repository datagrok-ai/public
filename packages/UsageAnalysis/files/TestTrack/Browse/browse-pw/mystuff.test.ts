import { test, expect } from '@playwright/test';
import {
  CONTEXT_MENU,
  CONTEXT_MENU_ADD_FAVORITES,
  contextMenuItem,
  treeGroupByName,
  treeItemByName,
  treeNodeByPath,
  viewTabHandle,
} from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
} from './helpers';

const BASE: string = process.env.DATAGROK_URL!;

const REQUIRED_MYSTUFF = ['Recent', 'Favorites', 'Shared with me', 'My files'];

test.describe.configure({ mode: 'serial' });

test.describe('Browse My stuff (Browse-MyStuff-*)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelOpen(page);
  });

  test('Browse-MyStuff-01 — My stuff contains all required subgroups', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'My stuff');

    for (const name of REQUIRED_MYSTUFF) {
      const locator = treeNodeByPath(page, ['My-stuff', name]);
      await expect(locator, `Subgroup "${name}" should be present in My stuff`)
        .toHaveCount(1, { timeout: 5_000 });
    }

    await expectNoErrors(page, sink);
  });

  test('Browse-MyStuff-02 — Recent lists a recently opened entity', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    const tutorials = treeNodeByPath(page, ['Apps', 'Tutorials']);
    await tutorials.waitFor({ state: 'visible', timeout: 10_000 });
    await tutorials.click();
    await page.waitForTimeout(2500);

    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'My stuff');
    await expandTreeGroup(page, 'Recent');
    await page.waitForTimeout(1000);

    const recentTut = page.locator('.d4-tree-view-item-label, .d4-tree-view-group-label',
      { hasText: /^Tutorials$/ });
    expect(await recentTut.count(),
      'Tutorials should appear in My stuff > Recent after being opened').toBeGreaterThanOrEqual(1);

    await expectNoErrors(page, sink);
  });

  test('Browse-MyStuff-03 — Favorites lists an entity added via context menu', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'Apps');
    const target = treeNodeByPath(page, ['Apps', 'Tutorials']);
    await target.waitFor({ state: 'visible', timeout: 10_000 });
    const targetLabel = target.locator('.d4-tree-view-item-label, .d4-tree-view-group-label').first();

    await targetLabel.click({ button: 'right' });
    await expect(page.locator(CONTEXT_MENU)).toBeVisible({ timeout: 5_000 });
    if (await contextMenuItem(page, 'Remove from favorites').isVisible().catch(() => false)) {
      await contextMenuItem(page, 'Remove from favorites').click();
      await page.waitForTimeout(800);
    } else {
      await page.keyboard.press('Escape');
    }

    await targetLabel.click({ button: 'right' });
    await contextMenuItem(page, CONTEXT_MENU_ADD_FAVORITES).click();
    await page.waitForTimeout(1500);

    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'My stuff');
    await expandTreeGroup(page, 'Favorites');
    await page.waitForTimeout(800);
    const favEntry = page.locator('.d4-tree-view-item-label, .d4-tree-view-group-label',
      { hasText: /^Tutorials$/ });
    expect(await favEntry.count(),
      'Tutorials should appear under Favorites').toBeGreaterThanOrEqual(1);

    await expectNoErrors(page, sink);

    await ensureBrowsePanelOpen(page);
    await expandTreeGroup(page, 'Apps');
    await targetLabel.click({ button: 'right' });
    if (await contextMenuItem(page, 'Remove from favorites').isVisible().catch(() => false)) {
      await contextMenuItem(page, 'Remove from favorites').click();
    } else {
      await page.keyboard.press('Escape');
    }
  });

test('Browse-MyStuff-05 — Add to Favorites is reachable from My stuff > My files', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'My stuff');
    const myFiles = treeNodeByPath(page, ['My-stuff', 'My-files']);
    await myFiles.waitFor({ state: 'visible', timeout: 10_000 });

    const label = myFiles.locator('.d4-tree-view-item-label, .d4-tree-view-group-label').first();
    await label.click({ button: 'right' });
    const menu = page.locator(CONTEXT_MENU);
    await expect(menu).toBeVisible({ timeout: 5_000 });

    const add = contextMenuItem(page, CONTEXT_MENU_ADD_FAVORITES);
    const remove = contextMenuItem(page, 'Remove from favorites');
    const addVisible = await add.isVisible().catch(() => false);
    const removeVisible = await remove.isVisible().catch(() => false);
    expect(addVisible || removeVisible,
      'My stuff > My files context menu must expose Add or Remove from favorites (ref: GROK-19848)')
      .toBe(true);

    await page.keyboard.press('Escape');
    await expectNoErrors(page, sink);
  });

  test('Browse-MyStuff-06 — a script saved via UI lands under My stuff > My scripts', async ({ page }) => {
    const sink = watchErrors(page);

    const unique = `qa_autotest_${Date.now()}`;
    let createdId: string | null = null;

    async function expandByName(name: string): Promise<void> {
      const node = page.locator(`[name="${name}"]`);
      await node.waitFor({ state: 'visible', timeout: 10_000 });
      const tri = node.locator('.d4-tree-view-tri').first();
      const expanded = await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded'))
        .catch(() => false);
      if (!expanded) {
        await tri.click();
        await page.waitForTimeout(1500);
      }
    }

    try {

      await page.goto(`${BASE}/scripts`);
      await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
      await page.locator('button[name="button-New"]').click();
      await page.locator('.d4-menu-item-label, .d4-list-item label, .d4-list-item',
        { hasText: /^JavaScript Script\.\.\.$/ }).first().click();
      await page.waitForSelector('.CodeMirror', { timeout: 15_000 });

      await page.evaluate((name) => {
        const cm = document.querySelector('.CodeMirror');
        const cur = cm.CodeMirror.getValue();
        cm.CodeMirror.setValue(cur.replace(/\/\/name:.*/, `//name: ${name}`));
      }, unique);
      await page.waitForTimeout(800);

      const saveBtn = page.locator('button[name="button-Save"]');
      await expect(saveBtn, 'Save button must become enabled after the name edit')
        .not.toHaveClass(/disabled/, { timeout: 5_000 });
      await saveBtn.click();
      await page.waitForTimeout(2500);

      await page.goto(BASE);
      await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
      await ensureBrowsePanelOpen(page);
      await page.locator('[name="icon-sync"]').first().click(); 
      await page.waitForTimeout(2000);
      await expandTreeGroup(page, 'My stuff');

      const ts = unique.replace(/^qa_autotest_/, '');
      const variants = [
        `tree-My-stuff---qa-autotest-${ts}`,
        `tree-My-stuff---qa_autotest_${ts}`,
      ];
      const matchSelector = variants.map((n) => `[name="${n}"]`).join(', ');

      await expect.poll(async () => {
        const found = await page.locator(matchSelector).count();
        if (found > 0) return found;
        const loadMore = page.locator(
          '[name^="tree-My-stuff---"] >> text=/Load more|Show more/i'
        ).first();
        if (await loadMore.isVisible().catch(() => false)) {
          await loadMore.click().catch(() => undefined);
        }
        return 0;
      }, {
        message: `A newly-saved script "${unique}" should appear under My stuff after reload`,
        timeout: 20_000,
        intervals: [1000, 1500, 2000, 2500, 3000],
      }).toBeGreaterThanOrEqual(1);

      createdId = await page.evaluate(async (n) => {
        try {
          const list = await (window as any).grok.dapi.scripts.filter(`name LIKE "%${n}%"`).list();
          return list?.[0]?.id ?? null;
        } catch { return null; }
      }, ts);

      await expectNoErrors(page, sink);
    } finally {

      if (createdId) {
        await page.evaluate(async (id) => {
          try {
            const s = await (window as any).grok.dapi.scripts.find(id);
            if (s) await (window as any).grok.dapi.scripts.delete(s);
          } catch {}
        }, createdId);
      } else {

        await page.evaluate(async (n) => {
          try {
            const list = await (window as any).grok.dapi.scripts.filter(`name LIKE "%${n}%"`).list();
            for (const s of list ?? []) await (window as any).grok.dapi.scripts.delete(s);
          } catch {}
        }, unique.replace(/^qa_autotest_/, ''));
      }
    }
  });
});
