import { Page, expect } from '@playwright/test';

// Auth state file produced by the default `e2e/global-setup.ts` (dev login).
export const AUTH_STATE = 'e2e/.auth.json';

// Default Postgres credentials for the test_postgres / external-provider scenarios.
// `db.datagrok.ai` is reachable from the dev environment; the password is provided
// via the `DG_PG_PASSWORD` env var (.env / shell), never hard-coded here.
export const PG_SERVER = process.env.DG_PG_SERVER ?? 'db.datagrok.ai';
export const PG_PORT = process.env.DG_PG_PORT ?? '54322';
export const PG_DB = process.env.DG_PG_DB ?? 'northwind';
export const PG_LOGIN = process.env.DG_PG_LOGIN ?? 'datagrok';
export const PG_PASSWORD = process.env.DG_PG_PASSWORD ?? '';

// External-provider Postgres (port 54327, db `test`, user `superuser`).
export const PG_EXT_SERVER = process.env.DG_PG_EXT_SERVER ?? 'db.datagrok.ai';
export const PG_EXT_PORT = process.env.DG_PG_EXT_PORT ?? '54327';
export const PG_EXT_DB = process.env.DG_PG_EXT_DB ?? 'test';
export const PG_EXT_LOGIN = process.env.DG_PG_EXT_LOGIN ?? 'superuser';
export const PG_EXT_PASSWORD = process.env.DG_PG_EXT_PASSWORD ?? '';

// SPARQL test endpoint (public).
export const SPARQL_ENDPOINT =
  process.env.DG_SPARQL_ENDPOINT ?? 'http://data.ontotext.com/repositories/data-last';

// ---------------------------------------------------------------------------
// Page setup
// ---------------------------------------------------------------------------

/**
 * Navigate to home, wait for the Browse sidebar, suppress hover tooltips, and
 * make sure the Browse panel is visible. Mirrors `e2e/queries/helpers.ts:goHome`
 * — same readiness signals, same tooltip-disable trick, same Browse-tab activation.
 */
export async function goHome(page: Page): Promise<void> {
  await page.goto('/', { waitUntil: 'domcontentloaded', timeout: 30_000 });
  await page.locator('[name="Browse"]').first().waitFor({ state: 'visible', timeout: 60_000 });
  // Wait for #grok-preloader inside #rootDiv to detach. On a cold CI Datlas the
  // Browse sidebar becomes visible BEFORE the preloader is dismissed; while
  // present, the preloader covers the entire rootDiv subtree and intercepts
  // every click — including those inside dialogs — which manifests as
  // "subtree intercepts pointer events" timeouts on the first interactive step.
  await page.waitForFunction(
    () => document.querySelector('#grok-preloader, .grok-preloader') == null,
    undefined, { timeout: 90_000 },
  );
  await page.waitForTimeout(500);
  // Datagrok hover-tooltips intercept Playwright's stability checks on tree nodes.
  // The platform's tooltip JS still runs (handlers fire, content builds) — only
  // visual rendering and pointer interception are disabled for the page lifetime.
  // Additionally, make `#grok-preloader` transparent to pointer events: on cold
  // CI Datlas the preloader can re-appear AFTER goHome (e.g. while a freshly
  // opened dialog fetches its schema) and time-out the next click. The platform
  // JS still drives it; only the pointer-events intercept is neutralised.
  await page.addStyleTag({ content: `
    .d4-tooltip { display: none !important; }
    #grok-preloader, .grok-preloader { pointer-events: none !important; }
  ` });
  // Activate Browse so the Databases tree is actually visible — otherwise tree
  // helpers operate on hidden DOM and UI-mode reviewers see no navigation.
  const databasesRoot = treeNodeLocator(page, 'tree-Databases');
  if (!(await databasesRoot.isVisible().catch(() => false))) {
    await page.locator('[name="Browse"]').first().click();
    await databasesRoot.waitFor({ state: 'visible', timeout: 10_000 });
  }
}

/**
 * Resolve a tree-view node by its `name=` attribute, excluding the lazy-loading
 * "Show more" footer (`.d4-tree-view-list-more`) which the platform tags with
 * the same name as its parent group and would otherwise trigger strict-mode
 * locator violations.
 */
function treeNodeLocator(page: Page, nodeName: string) {
  return page.locator(`[name="${nodeName}"]:not(.d4-tree-view-list-more)`).first();
}

/** Expand any Browse tree node by its `name=` attribute (idempotent). */
export async function expandTreeNode(page: Page, nodeName: string): Promise<void> {
  const node = treeNodeLocator(page, nodeName);
  await node.waitFor({ state: 'visible', timeout: 15_000 });
  await node.scrollIntoViewIfNeeded();
  const tri = node.locator('.d4-tree-view-tri').first();
  const expanded = await tri.evaluate((el) => el.classList.contains('d4-tree-view-tri-expanded'))
    .catch(() => false);
  if (!expanded)
    await tri.click();
  await page.waitForTimeout(800);
}

/** Expand a connection provider (e.g. "Postgres", "MS SQL") inside Browse > Databases. */
export async function expandDbProvider(page: Page, provider: string): Promise<void> {
  // The top-level `Databases` group can be collapsed by default — expand it first
  // so the provider node is actually visible before we click its expander.
  await expandTreeNode(page, 'tree-Databases');
  await expandTreeNode(page, providerNodeName(provider));
}

/** Expand an individual DB connection (e.g. Postgres > test_postgres) under its provider. */
export async function expandDbConnection(
  page: Page, provider: string, connection: string,
): Promise<void> {
  await expandTreeNode(page, connectionNodeName(provider, connection));
}

/** Build the `name=` attribute for a provider tree node. */
export function providerNodeName(provider: string): string {
  return `tree-Databases---${provider.replace(/ /g, '-')}`;
}

/** Build the `name=` attribute for a connection tree node under its provider. */
export function connectionNodeName(provider: string, connection: string): string {
  return `${providerNodeName(provider)}---${connection.replace(/_/g, '-').replace(/ /g, '-')}`;
}

/**
 * Open the "Show all data sources" footer on the Databases group — exposes
 * providers that don't currently have any saved connection (e.g. `Sparql`).
 *
 * The footer carries `.d4-tree-view-list-more` and the text "Show more...";
 * a real click expands it inline.
 */
export async function showAllDatabaseProviders(page: Page): Promise<void> {
  await expandTreeNode(page, 'tree-Databases');
  // The "Show more / ..." footer at the bottom of the Databases group may be
  // rendered as a `.d4-tree-view-list-more` row, OR as a plain element whose
  // textContent is "..." or "Show more". Try class-based first; fall back to
  // visible-text scan inside the Databases subtree.
  let clicked = false;
  const moreByClass = page.locator('[name="tree-Databases"] .d4-tree-view-list-more').first();
  if (await moreByClass.isVisible({ timeout: 1_000 }).catch(() => false)) {
    await moreByClass.scrollIntoViewIfNeeded();
    await moreByClass.click();
    clicked = true;
  }
  else {
    clicked = await page.evaluate(() => {
      const root = document.querySelector('[name="tree-Databases"]');
      if (!root) return false;
      const candidates = Array.from(root.querySelectorAll('div, span, label')) as HTMLElement[];
      const more = candidates.find((el) => {
        const t = el.textContent?.trim() ?? '';
        return (t === '...' || t === 'Show more' || t === 'more') && el.offsetParent !== null;
      });
      if (!more) return false;
      more.scrollIntoView({ block: 'center' });
      more.click();
      return true;
    });
  }
  if (clicked) await page.waitForTimeout(1500);
}

/**
 * Expand a Browse-tree expander whose `name=` follows the `tree-expander-...`
 * convention. Used for nested groups whose tree-NODE selector is shadowed by a
 * wrapper div and so isn't reachable via `expandTreeNode`.
 *
 * SCOPE NOTE: the inline Schemas / Catalogs group on a DB connection is rendered
 * inside a `div-{Provider}-{Conn}-Schemas` / `...Catalogs` wrapper, and its inner
 * `.d4-tree-view-node` reuses the parent connection's `name=` — so we can't
 * target its expander via the usual `tree-Databases---Provider---Conn---Group`
 * path. The dedicated `tree-expander-...` element is unambiguous; clicking it
 * via JS avoids hover-interception flakes on collapsed groups (the platform
 * listens at the DOM level either way).
 */
export async function clickTreeExpander(page: Page, expanderName: string): Promise<void> {
  await page.waitForSelector(`[name="${expanderName}"]`, { state: 'attached', timeout: 15_000 });
  await page.evaluate((name) => {
    const el = document.querySelector(`[name="${name}"]`) as HTMLElement | null;
    if (el && !el.classList.contains('d4-tree-view-tri-expanded')) el.click();
  }, expanderName);
  await page.waitForTimeout(1200);
}

/**
 * Expand the inline `Catalogs` (or `Schemas`) group directly under a DB
 * connection. The wrapper is rendered with `[name="div-{Provider}-{Conn}-{Group}"]`
 * and its inner tree-view-node reuses the parent connection's `name=`, so we
 * cannot drive its expander via the usual `tree-Databases---Provider---Conn---Group`
 * selector. Click the inner triangle directly via JS — the platform listens at
 * the DOM level either way (same pattern as queries/helpers.ts:expandDbSchemas).
 */
export async function expandDbGroupWrapper(
  page: Page, provider: string, connServerName: string, group: 'Catalogs' | 'Schemas',
): Promise<string> {
  const wrapper = `[name="div-${provider.replace(/ /g, '-')}-${connServerName}-${group}"]`;
  // The wrapper lazy-loads after the connection is expanded — server fetch of
  // catalog/schema metadata can take 5–20s on a cold connection.
  await page.waitForSelector(wrapper, { state: 'attached', timeout: 45_000 });
  await page.evaluate((sel) => {
    const w = document.querySelector(sel);
    const innerNode = w?.querySelector(':scope > .d4-tree-view-node');
    (innerNode as HTMLElement | null)?.scrollIntoView({ block: 'center' });
    const tri = innerNode?.querySelector(':scope > .d4-tree-view-tri');
    if (tri && !tri.classList.contains('d4-tree-view-tri-expanded')) (tri as HTMLElement).click();
  }, wrapper);
  await page.waitForTimeout(1500);
  return wrapper;
}

/**
 * Expand a tree-view node by its visible label, scoped to a parent `name=`.
 * Catalog / schema groups under a DB connection lack stable `name=` attributes
 * on dev, so we resolve them by label inside the parent `tree-Databases---...`
 * subtree and click the expander triangle directly.
 */
export async function expandChildByLabel(
  page: Page, parentName: string, label: string,
): Promise<void> {
  await page.waitForSelector(`[name="${parentName}"]`, { state: 'attached', timeout: 15_000 });
  await page.evaluate(({ p, l }) => {
    const parent = document.querySelector(`[name="${p}"]`);
    if (!parent) throw new Error(`parent tree node not found: ${p}`);
    // Tree groups render `.d4-tree-view-group-label`; items render `.d4-tree-view-item-label`.
    const candidates = Array.from(parent.querySelectorAll(
      '.d4-tree-view-group-label, .d4-tree-view-item-label',
    )) as HTMLElement[];
    const match = candidates.find((el) => el.textContent?.trim() === l);
    if (!match) throw new Error(`tree node "${l}" not found under ${p}`);
    const node = match.closest('.d4-tree-view-node') as HTMLElement | null;
    if (!node) throw new Error(`tree-view-node ancestor missing for "${l}"`);
    node.scrollIntoView({ block: 'center' });
    const tri = node.querySelector(':scope > .d4-tree-view-tri') as HTMLElement | null;
    if (tri && !tri.classList.contains('d4-tree-view-tri-expanded')) tri.click();
  }, { p: parentName, l: label });
  await page.waitForTimeout(1500);
}

/**
 * Right-click a tree-view node located by its visible label inside a parent
 * `name=` subtree. Same purpose as `expandChildByLabel` — works around the
 * platform reusing `name=` attributes on Catalogs / Schemas nested groups.
 */
export async function rightClickChildByLabel(
  page: Page, parentName: string, label: string,
): Promise<void> {
  await page.evaluate(({ p, l }) => {
    const parent = document.querySelector(`[name="${p}"]`);
    if (!parent) throw new Error(`parent tree node not found: ${p}`);
    const candidates = Array.from(parent.querySelectorAll(
      '.d4-tree-view-group-label, .d4-tree-view-item-label',
    )) as HTMLElement[];
    const match = candidates.find((el) => el.textContent?.trim() === l);
    if (!match) throw new Error(`tree node "${l}" not found under ${p}`);
    const node = match.closest('.d4-tree-view-node') as HTMLElement | null;
    if (!node) throw new Error(`tree-view-node ancestor missing for "${l}"`);
    node.scrollIntoView({ block: 'center' });
    node.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
    }));
  }, { p: parentName, l: label });
  await page.waitForSelector('.d4-menu-popup', { timeout: 5_000 });
}

/** Click the connection label in the Browse tree to open its preview/queries view. */
export async function openDbConnectionView(
  page: Page, provider: string, connection: string,
): Promise<void> {
  const label = page
    .locator(
      `[name="${connectionNodeName(provider, connection)}"]:not(.d4-tree-view-list-more)`
      + ' .d4-tree-view-group-label',
    )
    .first();
  await label.waitFor({ state: 'visible', timeout: 15_000 });
  await label.scrollIntoViewIfNeeded();
  await label.click();
}

// ---------------------------------------------------------------------------
// Right-click + menus
// ---------------------------------------------------------------------------

/** Right-click a Browse tree node by its `name=` attribute and wait for the popup. */
export async function rightClickTreeNode(page: Page, nodeName: string): Promise<void> {
  const node = treeNodeLocator(page, nodeName);
  await node.waitFor({ state: 'visible', timeout: 15_000 });
  await node.scrollIntoViewIfNeeded();
  // Real right-click first — UI-mode shows the action visually. Datagrok renders
  // rich `.d4-tooltip` overlays on hover that occasionally intercept the
  // pointer-event check; fall back to a dispatched `contextmenu` event in that
  // case (the platform listens at the DOM level, the menu still pops up).
  try {
    await node.click({ button: 'right', position: { x: 80, y: 8 }, timeout: 5_000 });
  }
  catch {
    await node.evaluate((el) => el.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
    })));
  }
  await page.waitForSelector('.d4-menu-popup', { timeout: 5_000 });
}

/** Click a context-menu item by exact label text (scoped to the open popup). */
export async function clickMenuItemExact(page: Page, text: string): Promise<void> {
  const escaped = text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  await page.locator('.d4-menu-popup .d4-menu-item-label')
    .filter({ hasText: new RegExp(`^${escaped}$`) })
    .first()
    .click();
}

/** Read the visible context-menu items in a single shot (used by assertions). */
export async function readMenuItems(page: Page): Promise<string[]> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
    .map((el) => (el as HTMLElement).textContent?.trim() ?? '')
    .filter((s) => s.length > 0));
}

/** Close any currently-open context menu by clicking outside / pressing Escape. */
export async function closeMenuPopup(page: Page): Promise<void> {
  await page.keyboard.press('Escape');
  await page.waitForTimeout(150);
}

// ---------------------------------------------------------------------------
// Add / edit / delete connection dialog
// ---------------------------------------------------------------------------

/** Right-click a provider in Browse > Databases and choose "Add connection..." (a.k.a. "New connection..."). */
export async function openAddConnectionDialog(page: Page, provider: string): Promise<void> {
  await expandDbProvider(page, provider);
  await rightClickTreeNode(page, providerNodeName(provider));
  // Both labels exist in the codebase across versions ("Add connection..." vs
  // "New connection..."). Click whichever is present.
  const items = await readMenuItems(page);
  const label = items.find((t) => /^Add connection\.\.\.$|^New connection\.\.\.$/.test(t));
  if (!label) throw new Error(`Add/New connection menu item not found for ${provider}: got ${JSON.stringify(items)}`);
  await clickMenuItemExact(page, label);
  await page.locator('.d4-dialog').waitFor({ timeout: 10_000 });
}

/**
 * Type a value into a connection-dialog input by its `Caption` (e.g. "Name", "Server").
 *
 * Uses real keyboard typing into the Datagrok input host
 * (`[name="input-host-{Caption}"] input`). Triple-click selects the existing
 * value so the type call replaces it cleanly.
 */
export async function fillConnectionField(
  page: Page, host: string, value: string,
): Promise<void> {
  const input = page.locator(`.d4-dialog [name="input-host-${host}"] input`).first();
  await input.waitFor({ state: 'visible', timeout: 10_000 });
  await input.click({ clickCount: 3 });
  await page.keyboard.type(value);
  // Verify the input accepted the keystrokes — fail fast if a Dart input listener
  // dropped the event.
  await expect(input).toHaveValue(value, { timeout: 5_000 });
}

/** Read a connection-dialog input value (e.g. for assertions on prefilled fields). */
export async function readConnectionField(page: Page, host: string): Promise<string> {
  const input = page.locator(`.d4-dialog [name="input-host-${host}"] input`).first();
  return input.inputValue();
}

/**
 * Pick a value in a `<select>` rendered inside a connection-dialog input host.
 * Falls back to the first `<select>` in the dialog when the host scope is missing
 * — some single-field dialogs (e.g. "Select primary schema") don't tag the host.
 */
export async function selectConnectionField(
  page: Page, host: string, value: string,
): Promise<void> {
  const scoped = page.locator(`.d4-dialog [name="input-host-${host}"] select`).first();
  const fallback = page.locator('.d4-dialog select').first();
  const target = (await scoped.count()) ? scoped : fallback;
  await target.waitFor({ state: 'visible', timeout: 10_000 });
  await target.selectOption(value);
}

/** Click the TEST button in the connection dialog and wait for any balloon to appear. */
export async function clickConnectionTest(page: Page, timeout = 60_000): Promise<void> {
  // Clear stale balloons from prior actions so the next observed one is the
  // result of THIS test click (not, e.g., a leftover Activity toast from the
  // last navigation step). `.grok-balloon` is the class shared across versions.
  await page.evaluate(() => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
    .forEach((b) => (b as HTMLElement).remove()));
  await page.locator('.d4-dialog [name="button-TEST"]').click();
  // The platform fires a green or red balloon depending on the driver result.
  await page.waitForFunction(
    () => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
      .some((b) => (b as HTMLElement).textContent && (b as HTMLElement).textContent!.trim().length > 0),
    undefined,
    { timeout },
  );
}

/** Click the OK button in the connection dialog and wait for the dialog to detach. */
export async function clickConnectionOk(page: Page): Promise<void> {
  await page.locator('.d4-dialog [name="button-OK"]').click();
  await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 15_000 });
}

/** Click the SAVE button inside the connection dialog (used by Edit... in some versions). */
export async function clickConnectionSave(page: Page): Promise<void> {
  // Edit dialogs sometimes label the primary button "SAVE" instead of "OK".
  const ok = page.locator('.d4-dialog [name="button-OK"]');
  const save = page.locator('.d4-dialog [name="button-SAVE"]');
  if (await save.isVisible().catch(() => false))
    await save.click();
  else
    await ok.click();
  await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 15_000 });
}

/** Click CANCEL in the connection dialog (best-effort cleanup helper). */
export async function clickConnectionCancel(page: Page): Promise<void> {
  const cancel = page.locator('.d4-dialog [name="button-CANCEL"]').first();
  if (await cancel.isVisible({ timeout: 1_000 }).catch(() => false)) {
    await cancel.click();
    await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 5_000 }).catch(() => null);
  }
}

// ---------------------------------------------------------------------------
// Delete confirmation dialog
// ---------------------------------------------------------------------------

/** Right-click a tree node, click `Delete...`, confirm DELETE in the dialog. */
export async function deleteTreeNodeViaContext(page: Page, nodeName: string): Promise<void> {
  await rightClickTreeNode(page, nodeName);
  await clickMenuItemExact(page, 'Delete...');
  // Confirmation dialog: "Are you sure?" — buttons DELETE / CANCEL (or YES / NO).
  await page.locator('.d4-dialog').waitFor({ timeout: 10_000 });
  const del = page.locator('.d4-dialog [name="button-DELETE"]');
  const yes = page.locator('.d4-dialog [name="button-YES"]');
  if (await del.isVisible().catch(() => false))
    await del.click();
  else
    await yes.click();
  await page.locator('.d4-dialog').waitFor({ state: 'detached', timeout: 10_000 });
}

// ---------------------------------------------------------------------------
// Test connection (from context menu)
// ---------------------------------------------------------------------------

/** Right-click a connection node and pick "Test connection"; returns when a balloon appears. */
export async function testConnectionViaContext(
  page: Page, nodeName: string, timeout = 60_000,
): Promise<void> {
  // Clear any leftover balloons so the next call sees only the new result.
  await page.evaluate(() => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
    .forEach((b) => (b as HTMLElement).remove()));
  await rightClickTreeNode(page, nodeName);
  await clickMenuItemExact(page, 'Test connection');
  await page.waitForFunction(
    () => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
      .some((b) => (b as HTMLElement).textContent && (b as HTMLElement).textContent!.trim().length > 0),
    undefined,
    { timeout },
  );
}

/** Read the concatenated text of all currently-visible balloons. */
export async function getAllBalloonsText(page: Page): Promise<string> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.grok-balloon, .d4-balloon'))
    .map((b) => (b as HTMLElement).textContent ?? '')
    .join(' | '));
}

// ---------------------------------------------------------------------------
// Server-side helpers (cleanup / verification — JS API is the canonical path)
// ---------------------------------------------------------------------------

/**
 * Find a connection by friendly name. Friendly name is the user-visible label
 * (e.g. `test_postgres`) — distinct from the server-side stored name (`TestPostgres`).
 */
export async function findConnectionByFriendlyName(
  page: Page, friendlyName: string,
): Promise<{ id: string; name: string; friendlyName: string; dataSource: string } | null> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const cs = await g.dapi.connections.filter(`friendlyName = "${fn}"`).list();
    if (!cs.length) return null;
    const c = cs[0];
    return { id: c.id, name: c.name, friendlyName: c.friendlyName, dataSource: c.dataSource };
  }, friendlyName);
}

/** Delete every connection with a given friendly name (safe cleanup before/after). */
export async function deleteConnectionByFriendlyName(
  page: Page, friendlyName: string,
): Promise<number> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const cs = await g.dapi.connections.filter(`friendlyName = "${fn}"`).list();
    for (const c of cs) await g.dapi.connections.delete(c);
    return cs.length;
  }, friendlyName);
}

/**
 * Refresh the Browse tree so newly-created or renamed connections appear /
 * deleted ones disappear. Equivalent to clicking the Refresh icon in the Browse
 * toolbar, but using the platform refresh function — Browse panel re-fetches
 * on `goHome` reload anyway, this is the lighter option mid-test.
 */
export async function refreshBrowseTree(page: Page): Promise<void> {
  await page.evaluate(() => {
    const g = (window as unknown as { grok: any }).grok;
    // The Browse panel reacts to a tree-refresh event — `simpleMode = true` keeps
    // the `Browse` view active, so we simply trigger a dataSource refresh by
    // re-opening the panel via the sidebar.
    g.shell.windows.simpleMode = true;
  });
  await page.locator('[name="Browse"]').first().click();
  await page.waitForTimeout(500);
}

// ---------------------------------------------------------------------------
// Automation setup (selenium class, simple mode, close all)
// ---------------------------------------------------------------------------

/** Apply the Datagrok automation env (selenium body class, filter icons, simple mode, close-all). */
export async function applyAutomationSetup(page: Page): Promise<void> {
  await page.evaluate(() => {
    const g = (window as unknown as { grok: any }).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
  });
}

// ---------------------------------------------------------------------------
// Context Panel sections
// ---------------------------------------------------------------------------

/**
 * Make sure the Context Panel (right-side property panel) is visible.
 * The platform default is `showProperties = true`; the panel is opened by
 * selection events, not a dedicated UI toggle — idempotent state-set is the
 * canonical approach here (see grok-browser SKILL).
 */
export async function showContextPanel(page: Page): Promise<void> {
  await page.evaluate(() => {
    const g = (window as unknown as { grok: any }).grok;
    g.shell.windows.showProperties = true;
  });
}

/**
 * Click an accordion-pane header on the Context Panel by its visible text.
 * Returns the pane locator so callers can read its body text for assertions.
 */
export async function clickContextPanelSection(
  page: Page, section: string,
): Promise<{ paneTextContent: () => Promise<string> }> {
  const escaped = section.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  // Header may include a trailing count badge tucked right against the title
  // ("Activity6" / "Activity 6" / "Activity(12)") — match starts-with rather
  // than enforcing a separator after the section name.
  const header = page.locator('.d4-accordion-pane-header')
    .filter({ hasText: new RegExp(`^${escaped}`) })
    .first();
  await header.waitFor({ state: 'visible', timeout: 10_000 });
  await header.scrollIntoViewIfNeeded();
  await header.click();
  await page.waitForTimeout(400);
  return {
    paneTextContent: async () => page.evaluate((t) => {
      const headers = Array.from(document.querySelectorAll('.d4-accordion-pane-header'));
      const h = headers.find((el) => (el.textContent ?? '').trim().startsWith(t)) as HTMLElement | undefined;
      return h?.parentElement?.textContent ?? '';
    }, section),
  };
}
