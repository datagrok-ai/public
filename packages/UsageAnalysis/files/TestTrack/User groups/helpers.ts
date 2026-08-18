import { Page, Locator, Browser, expect } from '@playwright/test';
import {
  GALLERY_COUNTS,
  GALLERY_GRID,
  CONTEXT_MENU,
  NEW_USER_BUTTON,
  DIALOG_TITLE,
  dialogInput,
  dialogButton,
  galleryCardByName,
  gallerySearch,
  ribbonButtonByText,
} from './selectors';
import { RIBBON } from '../browse/selectors';

export {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  expandTreeGroup,
  watchErrors,
  expectNoErrors,
} from '../browse/helpers';

import { BASE, ensureContextPanelOpen } from '../browse/helpers';

export type PlatformView = 'Users' | 'Groups' | 'Roles';

const VIEW_URL: Record<PlatformView, RegExp> = {
  Users: /\/users\b/,
  Groups: /\/groups\b/,
  Roles: /\/roles\b/,
};

export async function openPlatformView(page: Page, name: PlatformView): Promise<void> {

  await page.goto(`${BASE}/${name.toLowerCase()}`, { waitUntil: 'domcontentloaded', timeout: 60_000 });
  await page.waitForSelector(RIBBON, { timeout: 60_000 });
  await page.waitForURL(VIEW_URL[name], { timeout: 15_000 });
  await ensureContextPanelOpen(page);
  await page.locator(GALLERY_GRID).first().waitFor({ state: 'visible', timeout: 20_000 });
  await page.waitForTimeout(400);
}

export async function readGalleryCount(page: Page): Promise<{ shown: number; total: number }> {
  const counter = page.locator(GALLERY_COUNTS).first();
  let m: RegExpMatchArray | null = null;
  for (let i = 0; i < 20; i++) {
    const text = (await counter.textContent().catch(() => ''))?.trim() ?? '';
    m = text.match(/(\d+)\s*\/\s*(\d+)/);
    if (m) break;
    await page.waitForTimeout(250);
  }
  if (!m) return { shown: NaN, total: NaN };
  return { shown: Number(m[1]), total: Number(m[2]) };
}

export async function searchGallery(page: Page, kind: 'users' | 'groups' | 'roles', text: string): Promise<void> {
  const input = gallerySearch(page, kind);
  await input.click();
  await input.fill(text);
  await page.waitForTimeout(1200);
}

export async function clearGallerySearch(page: Page, kind: 'users' | 'groups' | 'roles'): Promise<void> {
  const reset = page.locator('.d4-search-bar .d4-reset-input, .d4-search-bar [name="icon-times"]').first();
  if (await reset.isVisible().catch(() => false)) await reset.click();
  else await gallerySearch(page, kind).fill('');

  await page.waitForURL((url) => !url.search.includes('q='), { timeout: 10_000 }).catch(() => {});
  await page.waitForTimeout(500);
}

export async function searchAndWaitCard(page: Page, kind: 'users' | 'groups' | 'roles', name: string): Promise<void> {
  for (let i = 0; i < 8; i++) {
    await searchGallery(page, kind, name);
    const card = galleryCardByName(page, name);
    if (await card.isVisible().catch(() => false)) {
      await page.waitForTimeout(800); 
      if (await card.isVisible().catch(() => false)) return;
    }

    await clearGallerySearch(page, kind);
    const refresh = page.locator('.d4-search-bar [name="icon-sync"]').first();
    if (await refresh.isVisible().catch(() => false)) await refresh.click();
    await page.waitForTimeout(1500);
  }
  throw new Error(`Entity "${name}" not found in the ${kind} gallery after refresh retries`);
}

export async function openCardContextMenu(page: Page, name: string): Promise<Locator> {
  const card = galleryCardByName(page, name);
  const menu = page.locator(CONTEXT_MENU);
  for (let attempt = 0; attempt < 4; attempt++) {

    if (!(await card.isVisible().catch(() => false)))
      await card.waitFor({ state: 'visible', timeout: 10_000 }).catch(() => {});
    await card.scrollIntoViewIfNeeded().catch(() => {}); 

    if (await menu.first().isVisible().catch(() => false)) {
      await page.keyboard.press('Escape');
      await page.waitForTimeout(300);
    }
    await card.click({ button: 'right' }).catch(() => {});
    if (await menu.first().isVisible({ timeout: 2_500 }).catch(() => false)) return card;
    await page.waitForTimeout(500);
  }
  await expect(menu.first(), `Context menu for "${name}" should open`).toBeVisible({ timeout: 3_000 });
  return card;
}

export async function closeMenu(page: Page): Promise<void> {
  await page.keyboard.press('Escape');
  await page.waitForTimeout(300);
}

const PERM_DIALOG = '.d4-dialog';

export async function addMembershipBySearch(page: Page, name: string): Promise<void> {
  const input = page.locator(`${PERM_DIALOG} .d4-user-selector-input`);
  await input.click();
  await input.fill('');
  await input.pressSequentially(name, { delay: 60 });
  const addRow = page.locator(`${PERM_DIALOG} .membership-add-row`, { hasText: name }).first();
  await addRow.waitFor({ state: 'visible', timeout: 10_000 });
  await addRow.locator('button[name="button-Add"]').first().click();

  await expect(page.locator(`${PERM_DIALOG} .membership-row`, { hasText: name }).first())
    .toBeVisible({ timeout: 5_000 });
}

export async function removeMembershipRow(page: Page, name: string): Promise<void> {
  const row = page.locator(`${PERM_DIALOG} .membership-row`, { hasText: name }).first();
  await row.waitFor({ state: 'visible', timeout: 5_000 });
  await row.locator('button[name="button-Remove"]').first().click();
  await page.waitForTimeout(300);
}

export async function saveDialog(page: Page): Promise<void> {
  await page.locator(`${PERM_DIALOG} .ui-btn`, { hasText: /^SAVE$/i }).first().click();
  await page.locator(PERM_DIALOG).first().waitFor({ state: 'hidden', timeout: 8_000 }).catch(() => {});
}

export async function cancelDialog(page: Page): Promise<void> {
  await page.locator(`${PERM_DIALOG} .ui-btn`, { hasText: /^CANCEL$/i }).first().click();
  await page.locator(PERM_DIALOG).first().waitFor({ state: 'hidden', timeout: 5_000 }).catch(() => {});
}

export async function openUserMembershipDialog(
  page: Page, user: string, menuName: 'Groups...' | 'Roles...',
): Promise<void> {
  await openCardContextMenu(page, user);
  await page.locator(`${CONTEXT_MENU} .d4-menu-item[d4-name="${menuName}"]`).first().click();
  await page.locator(`${PERM_DIALOG} .d4-user-selector-input`).waitFor({ state: 'visible', timeout: 10_000 });
}

export async function getUserStatus(page: Page, login: string): Promise<string> {
  return page.evaluate(async (lg) => {
    const g = (window as any).grok;
    const u = await g.dapi.users.filter(`login="${lg}"`).first();
    const auth = g.dapi.token || '';
    const j = await (await fetch(`${location.origin}/api/users/${u.id}`, { headers: { Authorization: auth } })).json();
    return j.status as string;
  }, login);
}

export async function unblockUser(page: Page, login: string): Promise<void> {
  await page.evaluate(async (lg) => {
    const g = (window as any).grok;
    const u = await g.dapi.users.filter(`login="${lg}"`).first();
    const auth = g.dapi.token || '';
    const userJson = await (await fetch(`${location.origin}/api/users/${u.id}`, { headers: { Authorization: auth } })).json();
    await fetch(`${location.origin}/api/users/unblock`, {
      method: 'POST', headers: { 'content-type': 'application/json', Authorization: auth }, body: JSON.stringify(userJson),
    });
  }, login);
}

export async function openNewUserOption(page: Page, option: string): Promise<void> {
  await page.locator(NEW_USER_BUTTON).click();
  await page.waitForTimeout(400);
  await page.locator('.d4-menu-popup .d4-menu-item-label', { hasText: new RegExp(`^${option.replace(/[.]/g, '\\.')}$`) })
    .first().click();
  await page.waitForTimeout(600);
}

export async function selectCard(page: Page, name: string): Promise<Locator> {
  const card = galleryCardByName(page, name);

  for (let attempt = 0; attempt < 4; attempt++) {
    if (!(await card.isVisible().catch(() => false)))
      await card.waitFor({ state: 'visible', timeout: 10_000 }).catch(() => {});
    await card.scrollIntoViewIfNeeded().catch(() => {});
    try {
      await card.click({ timeout: 5_000 });
      await page.waitForTimeout(1200);
      return card;
    } catch { await page.waitForTimeout(500); }
  }
  await card.click();
  await page.waitForTimeout(1200);
  return card;
}

export async function createGroup(page: Page, name: string, description?: string): Promise<void> {
  await ribbonButtonByText(page, 'NEW GROUP').click();
  await page.locator(DIALOG_TITLE).filter({ hasText: /Create New Group/i }).waitFor({ state: 'visible', timeout: 10_000 });
  await dialogInput(page, 'Name').fill(name);
  if (description) await dialogInput(page, 'Description').fill(description);
  await dialogButton(page, 'OK').click();
  await page.waitForTimeout(1500);
}

export async function createRole(page: Page, name: string, description?: string): Promise<void> {
  await ribbonButtonByText(page, 'New Role').click();
  await page.locator(DIALOG_TITLE).filter({ hasText: /Create New Role/i }).waitFor({ state: 'visible', timeout: 10_000 });
  await dialogInput(page, 'Name').fill(name);
  if (description) await dialogInput(page, 'Description').fill(description);
  await dialogButton(page, 'OK').click();
  await page.waitForTimeout(1500);
}

export async function deleteEntityViaContextMenu(page: Page, name: string): Promise<void> {
  await openCardContextMenu(page, name);
  await page.locator(`${CONTEXT_MENU} .d4-menu-item[d4-name="Delete"]`).first().click();

  const yes = page.locator('.d4-dialog button[name="button-YES"], .d4-dialog .ui-btn').filter({ hasText: /^(YES|OK)$/i }).first();
  if (await yes.isVisible({ timeout: 3_000 }).catch(() => false)) await yes.click();
  await page.waitForTimeout(1200);
}

export async function apiDeleteGroupsByPrefix(page: Page, prefix: string): Promise<void> {
  await page.evaluate(async (p) => {
    const g = (window as any).grok;
    const auth = g.dapi.token || '';
    const list = await g.dapi.groups.filter(`name like "${p}%"`).list();
    for (const grp of list) {
      try {
        await fetch(`${location.origin}/api/entities/${grp.id}`, { method: 'DELETE', headers: { Authorization: auth } });
      } catch {  }
    }
  }, prefix);
}

export async function sweepGroupsByPrefix(browser: Browser, prefix: string): Promise<void> {
  const ctx = await browser.newContext({ storageState: 'e2e/.auth.json' });
  try {
    const page = await ctx.newPage();
    await page.goto(BASE, { waitUntil: 'domcontentloaded', timeout: 60_000 });
    await page.waitForSelector(RIBBON, { timeout: 60_000 });
    await apiDeleteGroupsByPrefix(page, prefix);
  } finally {
    await ctx.close();
  }
}

export async function createGroupAsAdmin(browser: Browser, name: string): Promise<string> {
  const ctx = await browser.newContext({ storageState: 'e2e/.auth.json' });
  try {
    const page = await ctx.newPage();
    await page.goto(BASE, { waitUntil: 'domcontentloaded', timeout: 60_000 });
    await page.waitForSelector(RIBBON, { timeout: 60_000 });
    return await page.evaluate(async (n) => {
      const g = (window as any).grok;
      const DG = (window as any).DG;
      const saved = await g.dapi.groups.save(DG.Group.create(n));
      return (saved.friendlyName || saved.name) as string;
    }, name);
  } finally {
    await ctx.close();
  }
}

export async function openManageFromPane(page: Page, paneName: 'Members' | 'Assigned to'): Promise<void> {
  const pane = page.locator(`.grok-prop-panel .d4-accordion-pane[name="pane-${paneName.replace(/ /g, '-')}"]`).first();
  await pane.waitFor({ state: 'visible', timeout: 10_000 });
  const manage = pane.getByText('MANAGE', { exact: true }).first();
  if (!(await manage.isVisible().catch(() => false))) {
    await pane.locator('.d4-accordion-pane-header').first().click();
    await page.waitForTimeout(600);
  }
  await manage.click();
  await page.locator(`${PERM_DIALOG} .d4-user-selector-input`).waitFor({ state: 'visible', timeout: 10_000 });
}

export async function setMemberRowToggle(page: Page, name: string, on: boolean): Promise<void> {
  const row = page.locator(`${PERM_DIALOG} .membership-row`, { hasText: name }).first();
  const cb = row.locator('input[type="checkbox"]').first();
  if ((await cb.isChecked()) !== on) await cb.click();
}
