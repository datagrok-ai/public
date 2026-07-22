import { Page, expect } from '@playwright/test';

export const BASE = process.env.DATAGROK_URL!;

/**
 * Go to the home page and wait for the ribbon. Self-heals the session: if the login form is shown
 * (the shared `.auth.json` token can be invalidated by the relogin scenario's logout, a redeploy, or
 * expiry), log in through the UI.
 */
export async function gotoHome(page: Page): Promise<void> {
  await page.goto(BASE);
  const ribbon = page.locator('.d4-ribbon').first();
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  await Promise.race([
    ribbon.waitFor({ timeout: 90_000 }).catch(() => {}),
    loginInput.waitFor({ state: 'visible', timeout: 90_000 }).catch(() => {}),
  ]);
  if (await loginInput.isVisible().catch(() => false)) {
    await loginInput.fill(process.env.DATAGROK_LOGIN!);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).fill(process.env.DATAGROK_PASSWORD!);
    await page.keyboard.press('Enter');
  }
  await ribbon.waitFor({ timeout: 90_000 });
}

/** Standard automation setup: selenium class, tabs mode, close all views. */
export async function setupEnv(page: Page): Promise<void> {
  await page.evaluate(() => {
    const g = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
  });
  await page.waitForTimeout(400);
}

/** Unique suffix for this worker/run so parallel/leftover state never collides. */
export function uniqueSuffix(): string {
  return `${process.pid}_${Math.floor(Math.random() * 1e6)}`;
}

/**
 * Create a sticky-meta schema via the JS API (setup helper).
 * `matchBy` defaults to a molecule semtype matcher.
 */
export async function apiCreateSchema(
  page: Page,
  name: string,
  properties: { name: string; type: string }[],
  matchBy = 'semtype=molecule',
): Promise<void> {
  await page.evaluate(async ({ name, properties, matchBy }) => {
    const g = (window as any).grok;
    await g.dapi.stickyMeta.createSchema(name, [{ name: name + '_type', matchBy }], properties);
  }, { name, properties, matchBy });
}

/**
 * Delete a sticky-meta schema by name via the JS API (cleanup helper).
 * The JS Schema wrapper has no `id` getter, so the id is read from the Dart handle.
 * `deleteSchema` removes only the schema, NOT the inline entity types it created — those are
 * deleted explicitly via the generic entity-delete endpoint so no orphan types leak.
 */
export async function apiDeleteSchema(page: Page, name: string): Promise<void> {
  await page.evaluate(async (name) => {
    const g = (window as any).grok;
    const getId = (window as any).grok_Entity_Get_Id;
    const schema = (await g.dapi.stickyMeta.getSchemas()).find((s: any) => s.name === name);
    if (!schema) return;
    const typeIds = (schema.entityTypes || []).map((t: any) => getId(t.dart));
    await g.dapi.stickyMeta.deleteSchema(getId(schema.dart));
    const auth = g.dapi.token || '';
    for (const tid of typeIds)
      await fetch(`${location.origin}/api/entities/types/${tid}`, { method: 'DELETE', headers: { Authorization: auth } });
  }, name);
}

/**
 * Delete every leftover `PW_SM_*` schema via the API (defensive pre-cleanup).
 * Guarantees the test's own schema is the only molecule-matching one, so the cell-edit dialog
 * and Sticky-meta pane show a single section. Also deletes each schema's inline entity types
 * (deleteSchema leaves them behind), preventing orphan-type buildup across runs.
 */
export async function apiDeleteAllTestSchemas(page: Page): Promise<void> {
  await page.evaluate(async () => {
    const g = (window as any).grok;
    const getId = (window as any).grok_Entity_Get_Id;
    const auth = g.dapi.token || '';
    for (const s of (await g.dapi.stickyMeta.getSchemas()).filter((x: any) => /^PW_SM_/.test(x.name))) {
      const typeIds = (s.entityTypes || []).map((t: any) => getId(t.dart));
      await g.dapi.stickyMeta.deleteSchema(getId(s.dart));
      for (const tid of typeIds)
        await fetch(`${location.origin}/api/entities/types/${tid}`, { method: 'DELETE', headers: { Authorization: auth } });
    }
  });
}

/** Whether a schema with the given name currently exists (API read). */
export async function apiSchemaExists(page: Page, name: string): Promise<boolean> {
  return page.evaluate(async (name) => {
    const g = (window as any).grok;
    return (await g.dapi.stickyMeta.getSchemas()).some((s: any) => s.name === name);
  }, name);
}

/** Open SPGI.csv via the JS API and wait for the molecule grid to render. */
export async function openSpgi(page: Page): Promise<void> {
  await page.evaluate(async () => {
    const g = (window as any).grok;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 4000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 4000));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({ timeout: 30_000 });
}

/**
 * Navigate to a Sticky-Meta sub-view (`Types` or `Schemas`) through the Browse tree.
 * `goto()` between the sibling routes is ignored by the SPA, so we click the tree nodes.
 */
export async function gotoStickyMeta(page: Page, tab: 'Types' | 'Schemas'): Promise<void> {
  const click = () => page.evaluate(async (tab) => {
    const findLabel = (t: string) =>
      Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label'))
        .find((el) => el.textContent?.trim() === t) as HTMLElement | undefined;
    const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
    (window as any).grok.shell.windows.showBrowse = true;
    await wait(400);
    if (!findLabel('Sticky Meta')) { findLabel('Platform')?.click(); await wait(800); }
    if (!findLabel(tab)) { findLabel('Sticky Meta')?.click(); await wait(800); }
    findLabel(tab)?.click();
    await wait(1500);
    return location.href;
  }, tab);

  const route = tab === 'Types' ? /\/meta\/types(\?|$)/ : /\/meta\/schemas(\?|$)/;
  const button = tab === 'Types' ? '[name="button-New-Entity-Type..."]' : '[name="button-New-Schema..."]';
  for (let attempt = 0; attempt < 4; attempt++) {
    const url = await click();
    if (route.test(url)) break;
    await page.waitForTimeout(1500);
  }
  await page.locator(button).waitFor({ timeout: 30_000 });
  await expect(page).toHaveURL(route);
}

/** Type a value into the list's search box (the lists are paginated). */
export async function searchList(page: Page, query: string): Promise<void> {
  await page.evaluate((q) => {
    const search = document.querySelector('input[placeholder*="by name"]') as HTMLInputElement | null;
    if (!search) return;
    search.focus();
    const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
    setter.call(search, q);
    search.dispatchEvent(new Event('input', { bubbles: true }));
  }, query);
  await page.waitForTimeout(1500);
}

/** Whether a list card with the exact label exists (search first for paginated lists). */
export async function listHasCard(page: Page, label: string): Promise<boolean> {
  return page.evaluate((label) =>
    Array.from(document.querySelectorAll('label')).some((l) => l.textContent?.trim() === label),
  label);
}

/** Right-click a list card by label and click a context-menu item (e.g. Edit, Delete). */
export async function cardContextMenu(page: Page, label: string, item: string): Promise<void> {
  await page.evaluate(async ({ label, item }) => {
    const card = Array.from(document.querySelectorAll('label'))
      .find((l) => l.textContent?.trim() === label) as HTMLElement | undefined;
    if (!card) return;
    const rect = card.getBoundingClientRect();
    card.dispatchEvent(new MouseEvent('contextmenu',
      { bubbles: true, cancelable: true, button: 2, clientX: rect.left + 10, clientY: rect.top + 5 }));
    await new Promise((r) => setTimeout(r, 500));
    const menu = document.querySelector('.d4-menu-popup, .d4-menu');
    const mi = Array.from(menu?.querySelectorAll('.d4-menu-item') || [])
      .find((el) => (el.textContent?.trim() || '').startsWith(item)) as HTMLElement | undefined;
    mi?.click();
  }, { label, item });
  await page.waitForTimeout(800);
}

/** Click a dialog button (OK/CANCEL/DELETE/…) via a real mouse click at its box. */
export async function clickDialogButton(page: Page, dialogName: string, buttonName: string): Promise<void> {
  const box = await page.locator(`.d4-dialog[name="${dialogName}"] [name="${buttonName}"]`).boundingBox();
  if (box) await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
}
