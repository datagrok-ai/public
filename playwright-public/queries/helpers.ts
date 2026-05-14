import { Page, expect } from '@playwright/test';

// Connection name constant — easy to swap between environments.
//
// CI: the ephemeral Datlas auto-provisions the platform's own metadata DB
// as a Postgres connection named "Datagrok" (System:Datagrok, see
// `ServiceConnectionsMigration.createDatagrokConnection`). It always
// exists, contains the `public` schema with platform tables
// (`users`, `groups`, `entities`, …), and is reachable from grok_connect
// inside the test docker network — so queries-suite specs target it
// instead of Northwind / NorthwindTest, which exist only on dev/public.
// MS_SQL_CONNECTION isn't used in the CI flow (mssql-query-lifecycle is
// skipped — no MS SQL service in the CI compose); kept for parity with
// the dev playwright-tests/ copy.
export const POSTGRES_CONNECTION = 'Datagrok';
export const MS_SQL_CONNECTION = 'Northwind';

// Auth state file produced by queries/global-setup.ts (public-env login).
// Used by beforeAll/afterAll hooks that build their own browser context.
export const AUTH_STATE = 'e2e/.auth.public.json';

export async function goHome(page: Page): Promise<void> {
  // `'/'` resolves against `use.baseURL` from the active config — works for both
  // dev and public profiles without hardcoding a URL here. `domcontentloaded`
  // (instead of the default `load`) avoids waiting on every late-loaded asset on
  // slower envs like `public`; the ribbon-selector wait below confirms readiness.
  await page.goto('/', { waitUntil: 'domcontentloaded', timeout: 30_000 });
  // Wait for the sidebar Browse tab — it's the most reliable readiness signal.
  // (`.d4-ribbon` may be empty on home view, and `view-handle: Browse` only exists
  //  when Browse is opened as a view rather than the side panel.)
  await page.locator('[name="Browse"]').first().waitFor({ state: 'visible', timeout: 60_000 });
  // Wait for #grok-preloader to detach. On a cold CI Datlas the Browse sidebar
  // becomes visible BEFORE the preloader is dismissed; while present, the
  // preloader covers rootDiv and intercepts every click (including inside
  // dialogs), surfacing as "subtree intercepts pointer events" timeouts.
  await page.waitForFunction(
    () => document.querySelector('#grok-preloader, .grok-preloader') == null,
    undefined, { timeout: 90_000 },
  );
  await page.waitForTimeout(500);
  // Suppress Datagrok hover-tooltips for the duration of this page. They overlay tree
  // nodes and other clickable elements during Playwright's wait-for-stable check,
  // causing intermittent intercept failures. The platform's tooltip JS still fires
  // (handlers run, content is built) — only visual rendering and pointer interception
  // are disabled. `addStyleTag` lives in the page DOM, so it is re-injected after each
  // `goto('/')` automatically because every test calls `goHome`.
  // Additionally make `#grok-preloader` transparent to pointer events: it can
  // re-appear AFTER goHome on cold CI Datlas (e.g. while opening a dialog that
  // does an async server fetch) and time-out the next click. Platform JS still
  // drives it; only the pointer-events intercept is neutralised.
  await page.addStyleTag({ content: `
    .d4-tooltip { display: none !important; }
    #grok-preloader, .grok-preloader { pointer-events: none !important; }
  ` });
  // Activate the Browse sidebar so the Databases tree is actually visible. Otherwise
  // tree helpers operate on hidden DOM — the platform still handles the events, but
  // UI-mode timelines and human reviewers see no visual navigation.
  const databasesRoot = treeNodeLocator(page, 'tree-Databases');
  if (!(await databasesRoot.isVisible().catch(() => false))) {
    await page.locator('[name="Browse"]').first().click();
    await databasesRoot.waitFor({ state: 'visible', timeout: 10_000 });
  }
}

/**
 * Resolve a tree-view node by `name=` while excluding the lazy-loading "Show more"
 * footer (`.d4-tree-view-list-more`), which the platform tags with the same name as
 * its parent group and would otherwise trigger strict-mode locator violations.
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
  await page.waitForTimeout(1000);
}

/** Expand a connection provider (e.g. "Postgres", "MS SQL") inside Browse > Databases. */
export async function expandDbProvider(page: Page, provider: string): Promise<void> {
  // The top-level `Databases` group is collapsed by default on some envs (e.g. public).
  // Expand it first so the provider node is actually visible before we click it.
  await expandTreeNode(page, 'tree-Databases');
  await expandTreeNode(page, `tree-Databases---${provider.replace(/ /g, '-')}`);
}

/** Expand an individual DB connection (e.g. Postgres > Northwind) to show its saved queries. */
export async function expandDbConnection(page: Page, provider: string, connection: string): Promise<void> {
  await expandTreeNode(page, `tree-Databases---${provider.replace(/ /g, '-')}---${connection.replace(/ /g, '-')}`);
}

/** Click the connection label in the Browse tree to open its preview/queries view. */
export async function openDbConnectionView(page: Page, provider: string, connection: string): Promise<void> {
  const nodeName = `tree-Databases---${provider.replace(/ /g, '-')}---${connection.replace(/ /g, '-')}`;
  const label = page
    .locator(`[name="${nodeName}"]:not(.d4-tree-view-list-more) .d4-tree-view-group-label`)
    .first();
  await label.waitFor({ state: 'visible', timeout: 15_000 });
  await label.scrollIntoViewIfNeeded();
  await label.click();
}

/**
 * Ensure the Context Panel (right-side property panel) is visible.
 * The platform default is `showProperties = true` (see client_settings.dart),
 * and the panel is opened by selection events, not a dedicated UI toggle —
 * so per grok-browser SKILL ("Setup/cleanup → JS API allowed"), idempotent
 * state-set is the canonical approach here.
 */
export async function showContextPanel(page: Page): Promise<void> {
  await page.evaluate(() => {
    const g = (window as unknown as { grok: { shell: { windows: { showProperties: boolean } } } }).grok;
    g.shell.windows.showProperties = true;
  });
}

/** Switch to the Transformations tab inside the query editor view. */
export async function openTransformationsTab(page: Page): Promise<void> {
  await page.locator('[name="Transformations"]').first().click();
  await page.waitForSelector('[data-source="tab-content-Transformations"]', { state: 'attached', timeout: 10_000 });
  // Wait until the first pipeline step reports it has run — otherwise adding a new step fails
  // with the message "Previous step was not completed".
  await page.waitForFunction(() => {
    const firstStep = document.querySelector('.grok-action-editor');
    // When complete the step is not `.running` and has neither an error span nor a loader.
    return firstStep !== null && !firstStep.classList.contains('running');
  }, undefined, { timeout: 15_000 });
}

/** In the Transformations tab's Actions gallery, click an action by name and return when the dialog appears. */
export async function clickTransformationAction(page: Page, actionName: string): Promise<void> {
  // The action row is a bare `<label>`; Playwright's real mouse click at its coords doesn't reach
  // the handler (overlapping spans steal the hit). Dispatch a synthetic click on the label itself.
  await page.evaluate((name) => {
    const label = Array.from(document.querySelectorAll('.grok-actions-browser label'))
      .find((el) => el.textContent?.trim() === name) as HTMLElement | undefined;
    label?.click();
  }, actionName);
  await page.waitForSelector('.d4-dialog', { timeout: 10_000 });
}

/** Add an "Add New Column" transformation step with the given expression, then confirm the dialog. */
export async function addNewColumnTransformation(page: Page, expression: string): Promise<void> {
  await clickTransformationAction(page, 'Add New Column');
  // The formula editor in `_showAddNewColumnDialog`
  // (core/client/xamgle/lib/src/commands/edit/edit_add_new_column.dart:61)
  // is a plain `TextAreaInput` — a `<textarea class="ui-input-editor">`
  // with placeholder "Formula. Press '$' to select a column, or drag a
  // column here". Datagrok overlays the textarea with a formula-display
  // div, so the textarea itself is hidden — `state: 'visible'` fails (CI
  // build #24: "locator resolved to hidden <textarea ...>"). Wait for
  // attachment only, then set value programmatically and dispatch the
  // events `TextAreaInput.onChanged` listens for.
  const dialog = page.locator('.d4-dialog').first();
  await dialog.waitFor({ state: 'visible', timeout: 10_000 });
  const editor = dialog.locator('textarea').first();
  await editor.waitFor({ state: 'attached', timeout: 30_000 });
  await editor.evaluate((el, val) => {
    const ta = el as HTMLTextAreaElement;
    ta.focus();
    ta.value = val;
    ta.dispatchEvent(new Event('input', { bubbles: true }));
    ta.dispatchEvent(new Event('change', { bubbles: true }));
  }, expression);
  // The OK button used to be tagged `[name="button-Add-New-Column---OK"]` on
  // dev, but newer Datagrok builds drop the per-dialog prefix and the
  // footer button is simply `[name="button-OK"]` inside the Modal. Try
  // both, then fall back to the role/text-based locator so the dialog is
  // confirmed regardless of attribute naming.
  const okSelectors = [
    '[name="button-Add-New-Column---OK"]',
    '[name="button-OK"]',
  ];
  let okClicked = false;
  for (const sel of okSelectors) {
    const btn = dialog.locator(sel).first();
    if (await btn.isVisible({ timeout: 1_500 }).catch(() => false)) {
      await btn.click();
      okClicked = true;
      break;
    }
  }
  if (!okClicked)
    await dialog.getByRole('button', { name: /^OK$/ }).first().click();
  await expect(page.locator('.d4-dialog')).toHaveCount(0, { timeout: 10_000 });
}

/** Return the names of transformation pipeline steps currently listed in the Transformations tab. */
export async function transformationStepNames(page: Page): Promise<string[]> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.grok-action-editor'))
    .map((el) => el.querySelector('#name')?.textContent?.trim() ?? '')
    .filter((s) => s.length > 0));
}

/**
 * Delete the transformation step at the given zero-based index by clicking its delete icon.
 * The icon is rendered hover-only, so we hover the step first to reveal it, then `.click()`.
 * `force: true` is the safety net — if the platform skips opacity transitions on this row
 * (e.g. when fixed-tooltip styles override), Playwright's actionability check would otherwise
 * stall waiting for full visibility.
 */
export async function deleteTransformationStep(page: Page, index: number): Promise<void> {
  const step = page.locator('.grok-action-editor').nth(index);
  await step.scrollIntoViewIfNeeded();
  await step.hover();
  await step.locator('[name="icon-delete"]').click({ force: true });
  await page.waitForTimeout(500);
}

/** Find a project by its user-facing name (stored as `friendlyName`). */
export async function findProjectByFriendlyName(page: Page, friendlyName: string): Promise<{
  id: string; name: string; friendlyName: string;
} | null> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const ps = await g.dapi.projects.filter(`friendlyName = "${fn}"`).list();
    if (!ps.length) return null;
    return { id: ps[0].id, name: ps[0].name, friendlyName: ps[0].friendlyName };
  }, friendlyName);
}

/** Delete every project with the given friendlyName (safe cleanup). */
export async function deleteProjectByFriendlyName(page: Page, friendlyName: string): Promise<number> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const ps = await g.dapi.projects.filter(`friendlyName = "${fn}"`).list();
    for (const p of ps) await g.dapi.projects.delete(p);
    return ps.length;
  }, friendlyName);
}

/** Look up a DB connection's server-side name (e.g. "PostgresNorthwind") from its user-facing name. */
export async function getConnectionServerName(page: Page, provider: string, friendlyName: string): Promise<string> {
  return page.evaluate(async ({ p, f }) => {
    const g = (window as unknown as { grok: any }).grok;
    const c = (await g.dapi.connections.filter(`friendlyName = "${f}" and dataSource = "${p}"`).list())[0];
    return c?.name as string;
  }, { p: provider, f: friendlyName });
}

/**
 * Expand the "Schemas" group directly under a DB connection. The group uses an exceptional
 * `[name="div-{Provider}-{ConnServerName}-Schemas"]` wrapper, and its inner tree-view-node
 * reuses the same `name=` as the parent connection — so we can't target its expander via the
 * usual `tree-Databases---Provider---Conn` path.
 */
export async function expandDbSchemas(page: Page, provider: string, connServerName: string): Promise<void> {
  const groupSelector = `[name="div-${provider}-${connServerName}-Schemas"]`;
  await page.waitForSelector(groupSelector, { state: 'attached', timeout: 15_000 });
  await page.evaluate((sel) => {
    const group = document.querySelector(sel);
    const innerNode = group?.querySelector(':scope > .d4-tree-view-node');
    (innerNode as HTMLElement | null)?.scrollIntoView({ block: 'center' });
    const tri = innerNode?.querySelector(':scope > .d4-tree-view-tri');
    if (tri && !tri.classList.contains('d4-tree-view-tri-expanded')) (tri as HTMLElement).click();
  }, groupSelector);
  await page.waitForTimeout(1500);
}

/** Click a Browse-tree node in a way that registers it as the current object (updates the Context Panel). */
export async function selectTreeNodeAsCurrentObject(page: Page, nodeName: string): Promise<void> {
  const node = page.locator(`[name="${nodeName}"]:not(.d4-tree-view-list-more)`).first();
  await node.waitFor({ state: 'visible', timeout: 15_000 });
  await node.scrollIntoViewIfNeeded();
  // Aim past the expander triangle (clicks on `.d4-tree-view-tri` only toggle expand/collapse).
  // Same offset pattern used by `rightClickTreeNode` so behavior stays consistent.
  await node.click({ position: { x: 80, y: 8 } });
  await page.waitForTimeout(250);
}

/** Read the display name of the platform's current object (usually matches the clicked tree node). */
export async function getCurrentObjectName(page: Page): Promise<string | null> {
  return page.evaluate(() => {
    const g = (window as unknown as { grok: any }).grok;
    const o = g.shell.o;
    if (!o) return null;
    return o.name ?? String(o);
  });
}

/** Return visible error balloons (empty array = no errors). */
export async function getVisibleErrorBalloons(page: Page): Promise<string[]> {
  return page.evaluate(() => Array.from(document.querySelectorAll('.d4-balloon-error, .grok-balloon-error, [class*="balloon"][class*="error"]'))
    .filter((b) => (b as HTMLElement).getBoundingClientRect().width > 0)
    .map((b) => (b as HTMLElement).textContent?.trim() ?? '')
    .filter((s) => s.length > 0));
}

/** List the `name=` attributes of every column under a given DB table tree node. */
export async function listDbTableColumnNodeNames(page: Page, tableNodeName: string): Promise<string[]> {
  return page.evaluate((t) => Array.from(document.querySelectorAll(`[name^="${t}---"]`))
    .map((n) => n.getAttribute('name')!)
    .filter((s) => !!s), tableNodeName);
}

/** Click the `Browse` sidebar tab to bring the tree back into view after a TableView takes over. */
export async function focusBrowseSidebar(page: Page): Promise<void> {
  await page.locator('[name="Browse"]').first().click();
  await page.waitForTimeout(500);
}

/** Convert a query friendlyName to the suffix used in its tree `name=` attribute (underscores become hyphens). */
export function queryTreeNodeSuffix(friendlyName: string): string {
  return friendlyName.replace(/_/g, '-');
}

/** Right-click a Browse tree node by its `name=` attribute. */
export async function rightClickTreeNode(page: Page, nodeName: string): Promise<void> {
  const node = treeNodeLocator(page, nodeName);
  await node.waitFor({ state: 'visible', timeout: 15_000 });
  await node.scrollIntoViewIfNeeded();
  // Try a real right-click first so UI-mode shows the action visually.
  // Datagrok renders rich `.d4-tooltip` overlays on hover that often intercept the
  // pointer-event check; fall back to a dispatched `contextmenu` event in that case
  // (the menu still pops up correctly because the platform listens at the DOM level).
  try {
    // Aim past the expander triangle so the click lands on the label area
    // (clicks on `.d4-tree-view-tri` only toggle expand/collapse).
    await node.click({ button: 'right', position: { x: 80, y: 8 }, timeout: 5_000 });
  }
  catch {
    await node.evaluate((el) => el.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
    })));
  }
  // Wait for the context-menu popup itself — generic `.d4-menu-item-label` matches
  // the (hidden) top-menu placeholder before the popup renders.
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

/** Overwrite the query Name field in the query editor. */
export async function setQueryName(page: Page, name: string): Promise<void> {
  const input = page.locator('[name="input-Name"]');
  await input.waitFor({ timeout: 10_000 });
  // Triple-click selects the entire current value — works reliably against the Dart-managed
  // widget, whereas `fill()` and `pressSequentially` with delays lose keystrokes.
  await input.click({ clickCount: 3 });
  await page.keyboard.type(name);
  await expect(input).toHaveValue(name);
}

/**
 * UI-mode SQL entry: focus the CodeMirror editor and type the value via real
 * keyboard events. Used in "first" Adding tests to exercise the actual user
 * input path (Dart change handlers, autocomplete, syntax highlight). For
 * subsequent tests prefer `setQuerySql` — it's faster and not subject to
 * autocomplete/auto-pairs interference.
 */
export async function typeQuerySql(page: Page, sql: string): Promise<void> {
  await page.waitForSelector('.CodeMirror', { state: 'visible', timeout: 20_000 });
  // Ensure CM has finished mounting before typing — focus events are dropped otherwise.
  await page.waitForFunction(() => {
    const el = document.querySelector('.CodeMirror') as unknown as { CodeMirror?: unknown } | null;
    return !!el?.CodeMirror;
  }, undefined, { timeout: 10_000 });
  // The contenteditable surface inside CM is `.CodeMirror-code`; click there to focus.
  const surface = page.locator('.CodeMirror .CodeMirror-code').first();
  await surface.click();
  // Clear any existing content (Edit-like flows would have a previous SQL body).
  await page.keyboard.press('Control+A');
  await page.keyboard.press('Delete');
  await page.keyboard.type(sql, { delay: 15 });
  // Verify the editor accepted the keystrokes — fail fast if autocomplete/auto-pairs interfered.
  await expect.poll(async () => page.evaluate(() => {
    const el = document.querySelector('.CodeMirror') as unknown as { CodeMirror?: { getValue: () => string } } | null;
    return el?.CodeMirror?.getValue() ?? '';
  }), { timeout: 5_000 }).toBe(sql);
}

/** Replace the SQL body in the query editor's CodeMirror instance. */
export async function setQuerySql(page: Page, sql: string): Promise<void> {
  // CodeMirror mounts asynchronously after the editor view opens — attachment is enough to call its API.
  await page.waitForSelector('.CodeMirror', { state: 'attached', timeout: 20_000 });
  await page.waitForFunction(() => {
    const cm = document.querySelector('.CodeMirror') as unknown as { CodeMirror?: unknown } | null;
    return !!cm?.CodeMirror;
  }, undefined, { timeout: 10_000 });
  await page.evaluate((value) => {
    const cm = document.querySelector('.CodeMirror') as unknown as { CodeMirror: { setValue: (v: string) => void } };
    cm.CodeMirror.setValue(value);
  }, sql);
}

/**
 * Set `query.postProcessScript` for an already-saved query via the dapi (JS API fallback).
 *
 * UI fallback rationale (per grok-browser SKILL — "UI attempt failed → JS API, record reason"):
 * The Post-Process tab embeds a `ScriptView` whose live CodeMirror value is propagated into
 * `DataQueryView.postProcessScript` only via the Dart-side `scriptView.onScriptChanged`
 * listener (see `data_query_view.dart:850-861`). That listener fires off CM's native change
 * events, but neither `cm.setValue(...)` from JS nor `keyboard.insertText` propagated the
 * value reliably under headless Playwright — `save()` persisted the default template and
 * the post-process never ran on reopen. The CodeMirror UI-typing regression itself is
 * already covered elsewhere (`typeQuerySql`, `addNewColumnTransformation`); the tab-switch
 * widget is exercised by `openTransformationsTab`. So we wire the post-process body
 * through the public REST API and exercise the runtime behaviour (info balloon on Run)
 * which is the actual point of the test.
 */
export async function setQueryPostProcessScript(page: Page, friendlyName: string, script: string): Promise<void> {
  await page.evaluate(async ({ fn, s }) => {
    const grok = (window as unknown as { grok: any }).grok;
    const queries = await grok.dapi.queries.filter(`friendlyName = "${fn}"`).list();
    if (!queries.length) throw new Error(`query "${fn}" not found on server`);
    const q = queries[0];
    q.postProcessScript = s;
    await grok.dapi.queries.save(q, false);
  }, { fn: friendlyName, s: script });
}

/** Wait until the SQL editor is populated with the given text — used after Edit... to avoid racing against async load. */
export async function waitForQuerySql(page: Page, expected: string): Promise<void> {
  await page.waitForSelector('.CodeMirror', { state: 'attached', timeout: 20_000 });
  await expect.poll(async () => page.evaluate(() => {
    const cm = document.querySelector('.CodeMirror') as unknown as { CodeMirror?: { getValue: () => string } } | null;
    return cm?.CodeMirror?.getValue() ?? '';
  }), { timeout: 20_000 }).toBe(expected);
}

/** Click the Play icon in the query editor ribbon and wait for the result grid to appear. */
export async function runQueryViaPlay(page: Page): Promise<void> {
  await page.locator('[name="icon-play"]').first().click();
  // After Play the page may contain multiple `[name="viewer-Grid"]` nodes:
  // stale 0×0 placeholders from earlier views + the live result grid. A bare
  // selector resolves to the first match (often 0×0 and "visible" per
  // Playwright's check). Poll until at least one canvas under viewer-Grid
  // has non-zero size.
  const waitForLiveGrid = async (timeoutMs: number) => {
    await expect.poll(async () => page.evaluate(() => {
      const canvases = Array.from(document.querySelectorAll('[name="viewer-Grid"] canvas')) as HTMLCanvasElement[];
      return canvases.some((c) => c.width > 0 && c.height > 0);
    }), { timeout: timeoutMs }).toBe(true);
  };
  try {
    await waitForLiveGrid(30_000);
  }
  catch {
    // Public-env MS SQL occasionally drops the first call ("Socket closed" appears
    // in the Messages tab and no grid renders). A real user clicks Play again — the
    // second call typically succeeds against the now-warmed connection.
    await page.locator('[name="icon-play"]').first().click();
    await waitForLiveGrid(30_000);
  }
}

/** Click Toolbox > Actions > "Run query..." — result opens in a new table view. */
export async function runQueryViaActions(page: Page, queryName: string): Promise<void> {
  // After Play, the editor has Results/Messages sub-tabs that are also `view-handle:` — the
  // total handle count can even decrease when Run query collapses those. Instead count only
  // handles matching the query name: the editor adds one, a new table view adds another.
  const selector = `[name="view-handle: ${queryName}"]`;
  const before = await page.locator(selector).count();
  await page.locator('label.d4-link-action')
    .filter({ hasText: /^Run query\.\.\.$/ })
    .first()
    .click();
  await page.waitForFunction(({ sel, prev }) =>
    document.querySelectorAll(sel).length > prev, { sel: selector, prev: before }, { timeout: 30_000 });
}

/**
 * Switch back to the query editor tab (the one with `icon-data-query`).
 *
 * Implementation note: with `grok.shell.windows.simpleMode = true` (Tabs mode, our default),
 * the active view's content is what's painted — the inactive view-handles exist in the DOM
 * but Playwright's `.waitFor({ state: 'visible' })` rejects them as offscreen, so a real
 * `locator.click()` cannot be used. Per grok-browser SKILL ("UI attempt failed → fall back
 * to JS API, record reason"), dispatch the click at the DOM level — the platform listens at
 * the document level and switches the active tab regardless.
 */
export async function focusQueryEditorTab(page: Page, queryName: string): Promise<void> {
  await page.evaluate((name) => {
    const handles = Array.from(document.querySelectorAll(`[name="view-handle: ${name}"]`));
    const editor = handles.find((h) => h.querySelector('[name="icon-data-query"]')) as HTMLElement | undefined;
    editor?.click();
  }, queryName);
  await page.waitForTimeout(400);
}

/** Click the Save button in the query editor ribbon and wait for the server commit. */
export async function saveQuery(page: Page, friendlyName: string): Promise<void> {
  await page.locator('[name="button-Save"]').first().click();
  // Poll the server until the query is visible — Save is async and has no toast on this flow.
  await expect.poll(async () =>
    (await findQueryByFriendlyName(page, friendlyName)) !== null,
  { timeout: 30_000 }).toBe(true);
}

/** Find a saved query by the user-facing name (stored server-side as `friendlyName`). */
export async function findQueryByFriendlyName(page: Page, friendlyName: string): Promise<{
  name: string; friendlyName: string; connection: string | undefined;
} | null> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const qs = await g.dapi.queries.filter(`friendlyName = "${fn}"`).list();
    if (!qs.length) return null;
    const q = qs[0];
    return { name: q.name, friendlyName: q.friendlyName, connection: q.connection?.name };
  }, friendlyName);
}

/** Delete every saved query with a given friendlyName (API fallback — safe cleanup). */
export async function deleteQueryByFriendlyName(page: Page, friendlyName: string): Promise<number> {
  return page.evaluate(async (fn) => {
    const g = (window as unknown as { grok: any }).grok;
    const qs = await g.dapi.queries.filter(`friendlyName = "${fn}"`).list();
    for (const q of qs) await g.dapi.queries.delete(q);
    return qs.length;
  }, friendlyName);
}
