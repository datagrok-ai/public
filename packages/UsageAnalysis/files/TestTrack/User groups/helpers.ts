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

// Single import surface for the test files — re-export the Browse helpers we reuse.
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

/** URL fragment each Platform view routes to. */
const VIEW_URL: Record<PlatformView, RegExp> = {
  Users: /\/users\b/,
  Groups: /\/groups\b/,
  Roles: /\/roles\b/,
};

/**
 * Open a Platform leaf view (Users / Groups / Roles) through the Browse tree — the real
 * user path. Falls back to direct navigation only if the tree node can't be reached.
 * Leaves the gallery grid rendered and waits for the matching URL.
 */
export async function openPlatformView(page: Page, name: PlatformView): Promise<void> {
  // Open the view directly by route — boots the app and opens the gallery in a single navigation
  // (faster than Home → expand tree → click, and routing is a supported UI entry point). The
  // Browse-tree open path itself is covered by browse/platform.test.ts.
  await page.goto(`${BASE}/${name.toLowerCase()}`, { waitUntil: 'domcontentloaded', timeout: 60_000 });
  await page.waitForSelector(RIBBON, { timeout: 60_000 });
  await page.waitForURL(VIEW_URL[name], { timeout: 15_000 });
  await ensureContextPanelOpen(page);
  await page.locator(GALLERY_GRID).first().waitFor({ state: 'visible', timeout: 20_000 });
  await page.waitForTimeout(400);
}

/** Parse the "shown / total" gallery counter into numbers, polling until it settles. */
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

/** Type into the gallery search and wait for the count to settle. */
export async function searchGallery(page: Page, kind: 'users' | 'groups' | 'roles', text: string): Promise<void> {
  const input = gallerySearch(page, kind);
  await input.click();
  await input.fill(text);
  await page.waitForTimeout(1200);
}

/** Clear the gallery search via the reset (X) icon, then wait for the filter to drop. */
export async function clearGallerySearch(page: Page, kind: 'users' | 'groups' | 'roles'): Promise<void> {
  const reset = page.locator('.d4-search-bar .d4-reset-input, .d4-search-bar [name="icon-times"]').first();
  if (await reset.isVisible().catch(() => false)) await reset.click();
  else await gallerySearch(page, kind).fill('');
  // The query is reflected in the URL — wait until it's gone rather than guessing a delay.
  await page.waitForURL((url) => !url.search.includes('q='), { timeout: 10_000 }).catch(() => {});
  await page.waitForTimeout(500);
}

/**
 * Right-click an entity card by name and wait for the context menu to open. Re-queries the card on
 * each attempt and tolerates the gallery re-rendering (cards detach after search / save / delete).
 */
/**
 * Search for an entity by name and wait until a STABLE matching card is shown. A freshly created
 * group/role can lag the server search index (the card briefly appears from cache, then a server
 * search returns 0 and removes it), so this refreshes the gallery and retries until the card sticks.
 */
export async function searchAndWaitCard(page: Page, kind: 'users' | 'groups' | 'roles', name: string): Promise<void> {
  for (let i = 0; i < 8; i++) {
    await searchGallery(page, kind, name);
    const card = galleryCardByName(page, name);
    if (await card.isVisible().catch(() => false)) {
      await page.waitForTimeout(800); // let the server search settle, then confirm it stayed
      if (await card.isVisible().catch(() => false)) return;
    }
    // Not there (or flickered out) — clear, refresh the list, and retry.
    await clearGallerySearch(page, kind);
    const refresh = page.locator('.d4-search-bar [name="icon-sync"]').first();
    if (await refresh.isVisible().catch(() => false)) await refresh.click();
    await page.waitForTimeout(1500);
  }
  throw new Error(`Entity "${name}" not found in the ${kind} gallery after refresh retries`);
}

/** Right-click an entity card by name and wait for the context menu to open (retries). */
export async function openCardContextMenu(page: Page, name: string): Promise<Locator> {
  const card = galleryCardByName(page, name);
  const menu = page.locator(CONTEXT_MENU);
  for (let attempt = 0; attempt < 4; attempt++) {
    // Re-wait visible each attempt: a freshly searched / just-mutated gallery re-renders and
    // detaches cards, so scrollIntoView/click can hit a stale node. The card locator re-resolves
    // on every action, and the detach-prone calls are guarded.
    if (!(await card.isVisible().catch(() => false)))
      await card.waitFor({ state: 'visible', timeout: 10_000 }).catch(() => {});
    await card.scrollIntoViewIfNeeded().catch(() => {}); // tolerate "not attached" during re-render
    // Dismiss any stray popup that could intercept the right-click.
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

/** Close any open context menu / popup. */
export async function closeMenu(page: Page): Promise<void> {
  await page.keyboard.press('Escape');
  await page.waitForTimeout(300);
}

// ---------- Membership / assignment editor (MembershipEditor widget) ----------
// The same `.grok-permissions` widget backs: user Groups.../Roles..., group Members (MANAGE),
// and role Assigned to (MANAGE). Typing in `.d4-user-selector-input` runs an async search whose
// matches appear as `.membership-add-row` with a `button[name="button-Add"]`; existing rows are
// `.membership-row` with `button[name="button-Remove"]` and a per-row Admin/"Can assign" checkbox.

const PERM_DIALOG = '.d4-dialog';

/** Within an open membership/assignment dialog, search for `name` and click its Add (+) button. */
export async function addMembershipBySearch(page: Page, name: string): Promise<void> {
  const input = page.locator(`${PERM_DIALOG} .d4-user-selector-input`);
  await input.click();
  await input.fill('');
  await input.pressSequentially(name, { delay: 60 });
  const addRow = page.locator(`${PERM_DIALOG} .membership-add-row`, { hasText: name }).first();
  await addRow.waitFor({ state: 'visible', timeout: 10_000 });
  await addRow.locator('button[name="button-Add"]').first().click();
  // After adding, the item moves to the existing-members table.
  await expect(page.locator(`${PERM_DIALOG} .membership-row`, { hasText: name }).first())
    .toBeVisible({ timeout: 5_000 });
}

/** Within an open membership/assignment dialog, remove an existing member row by name. */
export async function removeMembershipRow(page: Page, name: string): Promise<void> {
  const row = page.locator(`${PERM_DIALOG} .membership-row`, { hasText: name }).first();
  await row.waitFor({ state: 'visible', timeout: 5_000 });
  await row.locator('button[name="button-Remove"]').first().click();
  await page.waitForTimeout(300);
}

/** Save an open dialog (SAVE) and wait for it to close. */
export async function saveDialog(page: Page): Promise<void> {
  await page.locator(`${PERM_DIALOG} .ui-btn`, { hasText: /^SAVE$/i }).first().click();
  await page.locator(PERM_DIALOG).first().waitFor({ state: 'hidden', timeout: 8_000 }).catch(() => {});
}

/** Cancel an open dialog and wait for it to close. */
export async function cancelDialog(page: Page): Promise<void> {
  await page.locator(`${PERM_DIALOG} .ui-btn`, { hasText: /^CANCEL$/i }).first().click();
  await page.locator(PERM_DIALOG).first().waitFor({ state: 'hidden', timeout: 5_000 }).catch(() => {});
}

/** Right-click a user card and open one of its membership editors (Groups... / Roles...). */
export async function openUserMembershipDialog(
  page: Page, user: string, menuName: 'Groups...' | 'Roles...',
): Promise<void> {
  await openCardContextMenu(page, user);
  await page.locator(`${CONTEXT_MENU} .d4-menu-item[d4-name="${menuName}"]`).first().click();
  await page.locator(`${PERM_DIALOG} .d4-user-selector-input`).waitFor({ state: 'visible', timeout: 10_000 });
}

/** Read a user's server-side status ('active' | 'blocked') by login. API read, for verification only. */
export async function getUserStatus(page: Page, login: string): Promise<string> {
  return page.evaluate(async (lg) => {
    const g = (window as any).grok;
    const u = await g.dapi.users.filter(`login="${lg}"`).first();
    const auth = g.dapi.token || '';
    const j = await (await fetch(`${location.origin}/api/users/${u.id}`, { headers: { Authorization: auth } })).json();
    return j.status as string;
  }, login);
}

/**
 * Unblock a user by login via the server API. The Users gallery has NO UI unblock action
 * (the context menu does not flip to "Unblock"), so this is the only way to revert a Block —
 * used strictly for test cleanup.
 */
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

/** Open the Users "NEW" dropdown and click one of its options. */
export async function openNewUserOption(page: Page, option: string): Promise<void> {
  await page.locator(NEW_USER_BUTTON).click();
  await page.waitForTimeout(400);
  await page.locator('.d4-menu-popup .d4-menu-item-label', { hasText: new RegExp(`^${option.replace(/[.]/g, '\\.')}$`) })
    .first().click();
  await page.waitForTimeout(600);
}

/** Select an entity card (single click) so the Context Panel follows it. */
export async function selectCard(page: Page, name: string): Promise<Locator> {
  const card = galleryCardByName(page, name);
  // The gallery re-renders after a mutation, detaching the card between resolve and click;
  // retry, tolerating the transient "not attached" detach.
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

/** Create a group via the ribbon "NEW GROUP..." dialog. Returns nothing; navigate back to verify. */
export async function createGroup(page: Page, name: string, description?: string): Promise<void> {
  await ribbonButtonByText(page, 'NEW GROUP').click();
  await page.locator(DIALOG_TITLE).filter({ hasText: /Create New Group/i }).waitFor({ state: 'visible', timeout: 10_000 });
  await dialogInput(page, 'Name').fill(name);
  if (description) await dialogInput(page, 'Description').fill(description);
  await dialogButton(page, 'OK').click();
  await page.waitForTimeout(1500);
}

/** Create a role via the ribbon "New Role..." dialog. */
export async function createRole(page: Page, name: string, description?: string): Promise<void> {
  await ribbonButtonByText(page, 'New Role').click();
  await page.locator(DIALOG_TITLE).filter({ hasText: /Create New Role/i }).waitFor({ state: 'visible', timeout: 10_000 });
  await dialogInput(page, 'Name').fill(name);
  if (description) await dialogInput(page, 'Description').fill(description);
  await dialogButton(page, 'OK').click();
  await page.waitForTimeout(1500);
}

/** Delete an entity (group / role) from its card context menu, confirming "Are you sure?". */
export async function deleteEntityViaContextMenu(page: Page, name: string): Promise<void> {
  await openCardContextMenu(page, name);
  await page.locator(`${CONTEXT_MENU} .d4-menu-item[d4-name="Delete"]`).first().click();
  // Optional confirmation dialog.
  const yes = page.locator('.d4-dialog button[name="button-YES"], .d4-dialog .ui-btn').filter({ hasText: /^(YES|OK)$/i }).first();
  if (await yes.isVisible({ timeout: 3_000 }).catch(() => false)) await yes.click();
  await page.waitForTimeout(1200);
}

// ---------- API cleanup for groups / roles ----------
// Roles are groups in the API (`grok.dapi.groups`, `personal: false`), so one call covers both.
// The UI lifecycle tests create `qa_autotest_g_*` / `qa_autotest_r_*` entities; deleting them only
// through the UI leaks them whenever a run fails mid-test, so every run also sweeps by name prefix.

/**
 * Delete every group/role whose name starts with `prefix`, via the API (uses the test's page).
 * Uses the generic entity-delete endpoint (`DELETE /api/entities/{id}`) rather than the typed
 * `groups.delete`: the latter throws a foreign-key error (23503) for a group that has a pending
 * membership request, while entity-delete removes the group and its requests together.
 */
export async function apiDeleteGroupsByPrefix(page: Page, prefix: string): Promise<void> {
  await page.evaluate(async (p) => {
    const g = (window as any).grok;
    const auth = g.dapi.token || '';
    const list = await g.dapi.groups.filter(`name like "${p}%"`).list();
    for (const grp of list) {
      try {
        await fetch(`${location.origin}/api/entities/${grp.id}`, { method: 'DELETE', headers: { Authorization: auth } });
      } catch { /* ignore individual failures */ }
    }
  }, prefix);
}

/**
 * Sweep all groups/roles matching `prefix` in a fresh authenticated context — for `beforeAll` /
 * `afterAll` hooks, which have no `page` fixture. Reuses the stored `e2e/.auth.json` session.
 */
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

/**
 * Create a group via the API in a fresh ADMIN context (`e2e/.auth.json`) and return its server
 * `friendlyName`. Used to seed a clean target group for the non-admin membership-request test:
 * a freshly created group's only admin is its creator (a real user), so the server's
 * notify-admins step has no null user to choke on. Groups are world-visible on dev, so the
 * non-admin sees it in the gallery.
 */
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

/**
 * Open the membership / assignment editor (MANAGE) from a Context-Panel pane of the selected entity:
 * 'Members' for a group, 'Assigned to' for a role. Expands the pane first if needed.
 */
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

/** Toggle the per-row Admin / "Can assign" checkbox for a member in an open membership dialog. */
export async function setMemberRowToggle(page: Page, name: string, on: boolean): Promise<void> {
  const row = page.locator(`${PERM_DIALOG} .membership-row`, { hasText: name }).first();
  const cb = row.locator('input[type="checkbox"]').first();
  if ((await cb.isChecked()) !== on) await cb.click();
}
