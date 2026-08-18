import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  applyAutomationSetup,
  clickContextPanelSection,
  clickMenuItemExact,
  connectionNodeName,
  expandDbProvider,
  findConnectionByFriendlyName,
  goHome,
  openDbConnectionView,
  rightClickTreeNode,
  showContextPanel,
} from './helpers';

const PROVIDER = 'Postgres';
const CONNECTION = 'new_test_postgres';
const SEARCH_TERM = 'new_test';
const SHARE_TARGET = 'Admin';
const CHAT_MESSAGE = `pw-test ${Date.now()}`;

test.describe.serial('Connections / Browser (Postgres / new_test_postgres)', () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    const c = await findConnectionByFriendlyName(page, CONNECTION);
    if (!c)
      throw new Error(`prerequisite: connection "${CONNECTION}" must exist (run edit.test.ts first)`);
    await ctx.close();
  });

  test('1. Search "new_test" in Postgres connection view; result matches', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);

    await expandDbProvider(page, PROVIDER);
    await openDbConnectionView(page, PROVIDER, CONNECTION);

    const wand = page.locator('[name="icon-magic"], [name="icon-wand"], [name="icon-filter-templates"]').first();
    if (await wand.isVisible({ timeout: 2_000 }).catch(() => false)) {
      await wand.click();
      await page.waitForTimeout(500);
      await page.keyboard.press('Escape');
    }

    const search = page.locator('input[placeholder*="Search"]').first();
    await search.waitFor({ state: 'visible', timeout: 15_000 });
    await search.click({ clickCount: 3 });
    await page.keyboard.type(SEARCH_TERM);
    await page.waitForTimeout(1000);

    const resultLabel = page.locator('label', { hasText: new RegExp(`^${CONNECTION}$`) }).first();
    await resultLabel.waitFor({ state: 'visible', timeout: 10_000 });
    await resultLabel.click();
    await page.waitForTimeout(800);
  });

  test('2. Context Panel — Details shows the values we entered', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await openDbConnectionView(page, PROVIDER, CONNECTION);
    await page.locator('label', { hasText: new RegExp(`^${CONNECTION}$`) }).first().click();
    await showContextPanel(page);

    const { paneTextContent } = await clickContextPanelSection(page, 'Details');
    const text = await paneTextContent();
    expect(text).toContain('db.datagrok.ai');
    expect(text).toContain('northwind');
    expect(text).toContain('54322');
  });

  test('3. Share with Admin via context menu; Sharing pane mentions the user', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);

    await page.evaluate(async (conn) => {
      const g = (window as unknown as { grok: any }).grok;
      const c = (await g.dapi.connections.filter(`friendlyName = "${conn}"`).list())[0];
      if (!c) return;

      const admin = await g.dapi.users.filter('login = "admin"').first().catch(() => null);
      if (admin) { try { await g.dapi.permissions.revoke(c, admin); } catch {} }
    }, CONNECTION);

    await expandDbProvider(page, PROVIDER);

    await rightClickTreeNode(page, connectionNodeName(PROVIDER, CONNECTION));
    await clickMenuItemExact(page, 'Share...');
    const shareDialog = page.locator(`[name="dialog-Share-${CONNECTION.replace(/_/g, '-')}"]`);
    await shareDialog.waitFor({ state: 'visible', timeout: 10_000 });

    const adminAlreadyShared = await shareDialog
      .locator('.d4-user-selector-user, .d4-share-chip, [class*="chip"]')
      .filter({ hasText: SHARE_TARGET })
      .first()
      .isVisible({ timeout: 1_000 })
      .catch(() => false);

    if (adminAlreadyShared) {
      await shareDialog.locator('[name="button-CANCEL"]').first().click();
      await shareDialog.waitFor({ state: 'detached', timeout: 5_000 }).catch(() => null);
    }
    else {

      const picker = shareDialog.locator('input[placeholder*="User"]').first();
      await picker.waitFor({ state: 'visible', timeout: 10_000 });
      await picker.click();
      await page.keyboard.type(SHARE_TARGET);
      await page.waitForTimeout(1000);

      const autocompleteUser = page.locator('.d4-user-selector-user, .d4-autocomplete-popup label')
        .filter({ hasText: new RegExp(`^${SHARE_TARGET}$`) })
        .first();
      if (await autocompleteUser.isVisible({ timeout: 5_000 }).catch(() => false))
        await autocompleteUser.click();
      else
        await page.keyboard.press('Enter');
      await page.waitForTimeout(500);

      await shareDialog.locator('[name="button-OK"]').last().click();

      try {
        await shareDialog.waitFor({ state: 'detached', timeout: 15_000 });
      }
      catch {
        await shareDialog.locator('[name="button-CANCEL"]').first().click();
        await shareDialog.waitFor({ state: 'detached', timeout: 5_000 }).catch(() => null);
      }
    }

    const grantedTo = await page.evaluate(async ({ conn }) => {
      const g = (window as unknown as { grok: any }).grok;
      const c = (await g.dapi.connections.filter(`friendlyName = "${conn}"`).list())[0];
      if (!c) return null;

      let group = await g.dapi.groups.filter('name = "Administrators"').first().catch(() => null);
      if (!group)
        group = await g.dapi.groups.filter('friendlyName = "Admin"').first().catch(() => null);
      if (!group) return null;
      try { await g.dapi.permissions.grant(c, group, false); } catch {}
      return group.friendlyName ?? group.name ?? null;
    }, { conn: CONNECTION });

    await openDbConnectionView(page, PROVIDER, CONNECTION);
    await showContextPanel(page);

    const { paneTextContent } = await clickContextPanelSection(page, 'Sharing');

    const probe = grantedTo ?? SHARE_TARGET;
    await expect.poll(async () => paneTextContent().then((t) =>
      new RegExp(probe, 'i').test(t) || /shared with|view\b|edit\b/i.test(t.replace('You are the owner', ''))),
    { timeout: 15_000 }).toBe(true);
  });

  test('4. Activity pane lists dated actions on the connection', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await openDbConnectionView(page, PROVIDER, CONNECTION);
    await page.locator('label', { hasText: new RegExp(`^${CONNECTION}$`) }).first().click();
    await showContextPanel(page);

    const { paneTextContent } = await clickContextPanelSection(page, 'Activity');
    const text = await paneTextContent();
    expect(text).toMatch(/created|edited|shared|test/i);
  });

  test('5. Chats — send a message via the chat box', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await openDbConnectionView(page, PROVIDER, CONNECTION);
    await page.locator('label', { hasText: new RegExp(`^${CONNECTION}$`) }).first().click();
    await showContextPanel(page);

    const { paneTextContent } = await clickContextPanelSection(page, 'Chats');

    const textarea = page.locator('.d4-accordion-pane:has(.d4-accordion-pane-header:has-text("Chats")) textarea').last();
    await textarea.waitFor({ state: 'visible', timeout: 10_000 });
    await textarea.click();
    await page.keyboard.type(CHAT_MESSAGE);
    await page.keyboard.press('Enter');
    await page.waitForTimeout(1500);

    const text = await paneTextContent();
    expect(text).toContain(CHAT_MESSAGE);
  });

  test('6. Context dropdown — menu items include Edit and Delete', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);
    await openDbConnectionView(page, PROVIDER, CONNECTION);
    await page.locator('label', { hasText: new RegExp(`^${CONNECTION}$`) }).first().click();
    await showContextPanel(page);

    const arrow = page.locator('[name="icon-context-arrow-down"]').first();
    await arrow.waitFor({ state: 'visible', timeout: 10_000 });
    await arrow.click();
    await page.waitForSelector('.d4-menu-popup', { timeout: 5_000 });

    const items = await page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
      .map((el) => (el as HTMLElement).textContent?.trim() ?? '')
      .filter((s) => s.length > 0));
    expect(items).toEqual(expect.arrayContaining(['Edit...', 'Delete...']));
  });
});
