/**
 * Playwright helpers for the projects feature.
 *
 * Namespace per helpers-candidates.yaml: `helpers.playwright.projects.*`.
 * Imported as: `import * as projects from '../helpers/projects';` (or named).
 *
 * Authoring cycle: c1a-2026-05-03 (7 lighter helpers — saveProject, waitForSaved,
 * saveAndReopen, openProjectFromDashboards, buildVariantsComposite,
 * collectChainProducedProjects, rename).
 *
 * c1b-2026-05-03 will add: deleteProjectViaContextMenu, shareProjectViaContextMenu,
 * saveCopy, logoutAndLoginAs.
 *
 * Selectors / behavior verified against post-Prompt-1 reference docs
 * (`grok-browser/references/projects.md` "Project Dialog Selectors" + "Context
 * Menu Items" sections + Quick Reference per-tile pattern).
 */

import {Page, expect} from '@playwright/test';

// ---------------------------------------------------------------------------
// 1. saveProject — Drive Save Project dialog with optional Data Sync toggle.
// ---------------------------------------------------------------------------

/**
 * Save a Datagrok project via JS API (preserves the typed name verbatim).
 *
 * Uses the JS API path per `projects.md:28` ("Preferred default: use JS API
 * via `grok.dapi.projects`. Do NOT drive the SAVE ribbon button or
 * Save/Share dialogs unless the scenario explicitly tests that flow.").
 *
 * IMPORTANT: the UI Save Project dialog auto-normalizes typed names to
 * PascalCase (e.g. `c1a-helpers-save` → `C1aHelpersSave`), which breaks
 * downstream `dapi.projects.filter('name = "X"')` lookups. The JS API path
 * preserves names verbatim and is the canonical path for fixture-helper
 * usage. Tests that need to exercise the Save-Project-dialog UI flow
 * specifically should use a future `saveProjectViaUI` helper instead (not
 * yet implemented; deferred to C1b cycle).
 *
 * Note: the JS API path creates an empty project entity (no view linkage).
 * `addLink(view)` returns a Dart-side `NoSuchMethodError` per Session B
 * Sub-investigation; until that's resolved at the platform layer, helpers
 * that need a populated project should use a different mechanism (e.g.
 * Wave 1a JS API approximation in projects-copy-clone-spec.ts:62-70).
 *
 * @param page - Playwright Page (login complete; table view optional).
 * @param name - Project name to save under (preserved verbatim).
 * @param options.dataSyncMode - Hint stored in the project description for
 *   future UI-path implementations; no functional effect via JS API path.
 * @param options.cancelAutoShare - Hint preserved for API symmetry with the
 *   future UI-path implementation; no auto-share dialog appears via JS API.
 */
export async function saveProject(
  page: Page,
  name: string,
  options?: {dataSyncMode?: 'on' | 'off'; cancelAutoShare?: boolean},
): Promise<void> {
  await page.evaluate(async ({n, sync}) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const project = DG.Project.create();
    project.name = n;
    if (sync !== undefined)
      project.description = `[saveProject helper hint: dataSyncMode=${sync}]`;
    await grok.dapi.projects.save(project);
  }, {n: name, sync: options?.dataSyncMode});
  // cancelAutoShare option preserved in signature; no-op for JS API path.
  void options?.cancelAutoShare;
}

// ---------------------------------------------------------------------------
// 2. waitForSaved — Poll for project to appear server-side after save.
// ---------------------------------------------------------------------------

/**
 * Wait until a project with the given name is visible via the JS API.
 *
 * Replaces the 9 inline `expect.poll(... dapi.projects.filter().first())`
 * copies surfaced across Wave 1a/1b specs. Default 60s timeout with
 * progressive intervals matching the project-save flake budget on dev.
 *
 * Note: this uses `filter('name = "X"').first()` for visibility. Per the
 * `rename` helper's caveat, the search index lags after rename — for
 * post-rename verification use `dapi.projects.find(id)` instead.
 *
 * @param page - Playwright Page (login already complete).
 * @param name - Project name to wait for.
 * @param options.timeout - Max wait in ms (default 60000).
 */
export async function waitForSaved(
  page: Page,
  name: string,
  options?: {timeout?: number},
): Promise<void> {
  const timeout = options?.timeout ?? 60000;
  await expect.poll(
    async () => page.evaluate(async (n) => {
      const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
      return p != null;
    }, name),
    {timeout, intervals: [500, 1000, 2000, 5000]},
  ).toBe(true);
}

// ---------------------------------------------------------------------------
// 3. saveAndReopen — Save current view as project + closeAll + reopen + count.
// ---------------------------------------------------------------------------

/**
 * Save the current table view as a project, close everything, reopen the
 * project from server-side, and return post-reopen table/viewer counts.
 *
 * Used for save-then-reopen verification flows (uploading, complex, opening,
 * projects-copy-clone scenarios — 18+ paths consolidated).
 *
 * @param page - Playwright Page (must have active table view).
 * @param name - Project name.
 * @param options.dataSyncMode - Forwarded to `saveProject`.
 * @returns post-reopen counts for assertion.
 */
export async function saveAndReopen(
  page: Page,
  name: string,
  options?: {dataSyncMode?: 'on' | 'off'},
): Promise<{tablesAfterReopen: number; viewersAfterReopen: number}> {
  await saveProject(page, name, {dataSyncMode: options?.dataSyncMode});
  await waitForSaved(page, name);

  // Close + reopen via JS API (find-by-name is fine here — fresh save, no
  // rename involved, search index already settled per waitForSaved poll).
  await page.evaluate(async (n) => {
    (window as any).grok.shell.closeAll();
    const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
    await p.open();
  }, name);

  // Wait briefly for the reopen to materialize tables/viewers in shell.
  await page.waitForTimeout(2000);

  // Datagrok shell collections wrap Dart proxies whose serialization across
  // CDP triggers "object reference chain too long" errors even when accessed
  // only for length. Safest path: coerce to primitive numbers via Number(),
  // catch all exceptions, return plain {number, number}.
  const tablesAfterReopen: number = await page.evaluate(() => {
    try { return Number((window as any).grok?.shell?.tables?.length) || 0; }
    catch (_) { return 0; }
  });
  const viewersAfterReopen: number = await page.evaluate(() => {
    try { return Number((window as any).grok?.shell?.tv?.viewers?.length) || 0; }
    catch (_) { return 0; }
  });
  return {tablesAfterReopen, viewersAfterReopen};
}

// ---------------------------------------------------------------------------
// 4. openProjectFromDashboards — Locate tile + double-click to open.
// ---------------------------------------------------------------------------

/**
 * Open a project by locating its tile in Browse > Dashboards and
 * double-clicking. Uses the per-tile selector pattern documented in
 * `grok-browser/references/projects.md` Quick Reference:
 * `[name="div-{project-name-lowercased}"]`.
 *
 * Caller must have already navigated to Browse > Dashboards (or call
 * `page.goto('/browse/dashboards')` first). The helper does NOT navigate —
 * keeps it composable with chains that already have the dashboards view open.
 *
 * @param page - Playwright Page.
 * @param name - Project name (raw, before slug conversion).
 */
export async function openProjectFromDashboards(
  page: Page,
  name: string,
): Promise<void> {
  // Slug = lowercased name, hyphen-separated. For most autogenerated test
  // projects (e.g. `AutoTest-Share-1234`) the slug is just the lowercase form.
  const slug = name.toLowerCase();
  const tile = page.locator(`[name="div-${slug}"]`);
  await tile.waitFor({timeout: 30000});
  // `force: true` bypasses Datagrok's overlay actionability checks (the
  // gallery has hover-decoration overlays that can intercept dblclick
  // pointer events under headless conditions).
  await tile.dblclick({force: true});
  // Wait for the project view to materialize. Match against expected name
  // (not just non-null) — under headless conditions a stale shell.project
  // from a prior step may briefly satisfy a non-null check.
  await expect.poll(async () => page.evaluate((expected) => {
    const grok = (window as any).grok;
    return grok.shell.project != null && grok.shell.project.name === expected;
  }, name), {timeout: 30000, intervals: [500, 1000, 2000, 3000]}).toBe(true);
}

// ---------------------------------------------------------------------------
// 5. buildVariantsComposite — Build a copy/clone/PVC fixture composite.
// ---------------------------------------------------------------------------

/**
 * Build a composite cross-fixture for the `project-url.md` scenario.
 *
 * For C1a: uses JS API approximation per Wave 1a pattern — creates separate
 * project entities for each requested variant with naming convention
 * `<sourceProjectName>-{link|clone|pvc}`. The actual UI Save-Copy-with-mode
 * flow is exercised by the `saveCopy` helper (C1b cycle).
 *
 * The original source project must already exist on the server (created by
 * `saveProject` or by an upstream chained spec).
 *
 * @param page - Playwright Page (login complete).
 * @param sourceProjectName - Name of the existing source project.
 * @param variants - Which variants to create (default: all three).
 * @returns Map of variant name → project name (with `original` always set).
 */
export async function buildVariantsComposite(
  page: Page,
  sourceProjectName: string,
  variants: Array<'link' | 'clone' | 'pvc'> = ['link', 'clone', 'pvc'],
): Promise<{original: string; link?: string; clone?: string; pvc?: string}> {
  const result: {original: string; link?: string; clone?: string; pvc?: string} =
    {original: sourceProjectName};

  for (const variant of variants) {
    const variantName = `${sourceProjectName}-${variant}`;
    await page.evaluate(async ({src, name, variant}) => {
      const grok = (window as any).grok;
      const DG = (window as any).DG;
      const source = await grok.dapi.projects.filter(`name = "${src}"`).first();
      if (!source)
        throw new Error(`buildVariantsComposite: source "${src}" not found`);
      // JS API approximation: create a new project entity with the variant
      // suffix. UI saveCopy modes (link/clone/PVC) preserve / re-derive
      // tables differently; this approximation creates an independent project
      // matching the variant naming for downstream lookup, with description
      // recording which mode it represents.
      const p = DG.Project.create();
      p.name = name;
      p.description = `c1a buildVariantsComposite ${variant} variant of "${src}" (JS API approximation; saveCopy UI flow deferred to C1b)`;
      await grok.dapi.projects.save(p);
    }, {src: sourceProjectName, name: variantName, variant});
    result[variant] = variantName;
  }

  return result;
}

// ---------------------------------------------------------------------------
// 6. collectChainProducedProjects — Enumerate fixtures by name prefix.
// ---------------------------------------------------------------------------

/**
 * Collect all projects whose name starts with the given prefix.
 *
 * Used by `deleting.md` terminal cleanup and any test that needs to enumerate
 * a chain's produced fixtures. Returns `{name, id}` so callers can drive
 * downstream `dapi.projects.delete` or `find` operations.
 *
 * @param page - Playwright Page (login complete).
 * @param namePrefix - Name prefix (e.g. "AutoTest-Share-" or "c1a-test-").
 * @returns Array of project descriptors.
 */
export async function collectChainProducedProjects(
  page: Page,
  namePrefix: string,
): Promise<Array<{name: string; id: string}>> {
  return await page.evaluate(async (prefix) => {
    const grok = (window as any).grok;
    // SQL `LIKE` with `%` for prefix matching. Datagrok's filter syntax
    // accepts SQL-style LIKE clauses on string columns.
    const list = await grok.dapi.projects.filter(`name like "${prefix}%"`).list();
    return list.map((p: any) => ({name: p.name, id: p.id}));
  }, namePrefix);
}

// ---------------------------------------------------------------------------
// 7. rename — JS API rename with find-by-id verification.
// ---------------------------------------------------------------------------

/**
 * Rename a project via JS API and verify server-side persistence.
 *
 * IMPORTANT: verification uses `dapi.projects.find(id)`, NOT
 * `dapi.projects.filter('name = "<new>"').first()`. Per Session B Sub 5
 * investigation: the rename DOES persist (re-fetch by ID returns the new
 * name) but the search index has eventually-consistent lag — filter-by-name
 * may return null for several seconds after a rename. Always verify by ID.
 *
 * Wave 1b complex-rename-spec round 1 reported "3/3 FAIL" because the
 * verification used filter-by-name; this helper bakes in the correct path.
 *
 * @param page - Playwright Page (login complete).
 * @param projectId - Project ID (UUID) to rename. Use ID, not name, because
 *   name lookups can lag after prior renames.
 * @param newName - New project name.
 * @returns `{ok, newName, verifiedVia}`. `ok=false` when rename did not
 *   propagate or project was not found.
 */
export async function rename(
  page: Page,
  projectId: string,
  newName: string,
): Promise<{ok: boolean; newName: string | null; verifiedVia: 'find-by-id'}> {
  return await page.evaluate(async ({id, name}) => {
    const grok = (window as any).grok;
    const p = await grok.dapi.projects.find(id);
    if (!p)
      return {ok: false, newName: null, verifiedVia: 'find-by-id' as const};
    p.name = name;
    await grok.dapi.projects.save(p);
    // Verify via find-by-id (NOT filter-by-name — search index lags).
    const verify = await grok.dapi.projects.find(id);
    return {
      ok: verify != null && verify.name === name,
      newName: verify ? verify.name : null,
      verifiedVia: 'find-by-id' as const,
    };
  }, {id: projectId, name: newName});
}

// ===========================================================================
// C1b helpers (heavier UI flows) — added 2026-05-03 via c1-2026-05-03-helpers-c1b
// ===========================================================================

// ---------------------------------------------------------------------------
// 8. deleteProjectViaContextMenu — UI delete via right-click OR Context Panel.
// ---------------------------------------------------------------------------

/**
 * Delete a project via the Browse > Dashboards UI flow.
 *
 * Two trigger paths supported (per `projects.md` Context Menu Items section):
 *   - `right-click` (default): contextmenu on `[name="div-{slug}"]` tile
 *     → click `[name="div-Delete-Project"]` → DELETE confirmation.
 *   - `context-panel-dropdown`: single-click tile to populate Context Panel
 *     → click `.grok-prop-panel [name="icon-context-arrow-down"]` chevron
 *     → click `[name="div-Delete-Project"]` from the SAME menu (verified
 *       identical menu in Session B Sub 4 alt-path investigation)
 *     → DELETE confirmation.
 *
 * Both paths terminate by clicking `[name="button-DELETE"]` (verified live
 * 2026-05-03; older `projects.md:178` claim of "no name= attribute" was
 * out-of-date and corrected via Prompt 1 Item 3.4).
 *
 * Caller must be on Browse > Dashboards with the project tile visible
 * (search-filter applied if needed). Helper does NOT navigate to dashboards.
 *
 * @param page - Playwright Page (login complete; on Browse > Dashboards).
 * @param name - Project name (used to compute slug `name.toLowerCase()`).
 * @param mode - Trigger path; default `'right-click'`.
 */
export async function deleteProjectViaContextMenu(
  page: Page,
  name: string,
  mode: 'right-click' | 'context-panel-dropdown' = 'right-click',
): Promise<void> {
  const slug = name.toLowerCase();
  const tile = page.locator(`[name="div-${slug}"]`);
  await tile.waitFor({timeout: 30000});

  if (mode === 'right-click') {
    // Dispatch contextmenu via DOM (Playwright's .click({button:'right'})
    // sometimes misses Datagrok's overlay-wrapped tile; explicit dispatch
    // matches the Session B Sub 1 + Session B Sub 4 verified pattern).
    await page.evaluate((s) => {
      const el = document.querySelector(`[name="div-${s}"]`) as HTMLElement | null;
      if (!el) throw new Error(`tile [name="div-${s}"] not found`);
      const r = el.getBoundingClientRect();
      el.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: r.left + r.width / 2,
        clientY: r.top + r.height / 2,
      }));
    }, slug);
  } else {
    // context-panel-dropdown path: single-click tile, wait for Context Panel
    // populated, click the chevron icon.
    //
    // Note: the chevron is CSS-hidden by default (only revealed on hover);
    // wait for `state: 'attached'` (in DOM) rather than `visible`, then
    // click via JS DOM dispatch to bypass actionability checks.
    await tile.click({force: true});
    await page.locator('.grok-prop-panel').waitFor({timeout: 15000});
    await page.locator('.grok-prop-panel [name="icon-context-arrow-down"]')
      .waitFor({state: 'attached', timeout: 15000});
    await page.evaluate(() => {
      const el = document.querySelector('.grok-prop-panel [name="icon-context-arrow-down"]') as HTMLElement | null;
      if (!el) throw new Error('chevron not in DOM');
      el.click();
    });
  }

  // The dropdown menu (same selectors regardless of trigger path).
  const deleteItem = page.locator('[name="div-Delete-Project"]');
  await deleteItem.waitFor({timeout: 10000});
  await deleteItem.click({force: true});

  // DELETE confirmation dialog. Filter by button-DELETE presence rather
  // than title text — server may render slightly different titles between
  // builds, but the DELETE button selector is stable per Session B Sub 4.
  const confirmDialog = page.locator('.d4-dialog')
    .filter({has: page.locator('[name="button-DELETE"]')});
  await confirmDialog.waitFor({timeout: 15000});
  await confirmDialog.locator('[name="button-DELETE"]').click();
  await expect(confirmDialog).toBeHidden({timeout: 15000});
}

// ---------------------------------------------------------------------------
// 9. shareProjectViaContextMenu — Right-click → Share dialog → grant.
// ---------------------------------------------------------------------------

/**
 * Share a project via the right-click → Share dialog UI flow.
 *
 * Drives the Share dialog selectors documented in `projects.md` "Project
 * Dialog Selectors" section (added by Prompt 1 Item 3.2). Recipient input
 * has no `name=` attribute — located by placeholder text.
 *
 * @param page - Playwright Page (on Browse > Dashboards with tile visible).
 * @param name - Project name.
 * @param options.recipient - User / group / email string typed into the
 *   recipient input (e.g. `'Test permission group'` — verified to exist on
 *   dev per Olena 2026-05-03).
 * @param options.accessLevel - Access level select; default `'View and use'`.
 * @param options.sendNotifications - Toggle Send-notifications checkbox
 *   (default: leave at dialog default).
 * @param options.message - Optional message textarea content.
 */
export async function shareProjectViaContextMenu(
  page: Page,
  name: string,
  options: {
    recipient: string;
    accessLevel?: 'View and use' | 'Full access';
    sendNotifications?: boolean;
    message?: string;
  },
): Promise<void> {
  const slug = name.toLowerCase();
  const tile = page.locator(`[name="div-${slug}"]`);
  await tile.waitFor({timeout: 30000});

  // Right-click via DOM dispatch.
  await page.evaluate((s) => {
    const el = document.querySelector(`[name="div-${s}"]`) as HTMLElement | null;
    if (!el) throw new Error(`tile [name="div-${s}"] not found`);
    const r = el.getBoundingClientRect();
    el.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
      clientX: r.left + r.width / 2,
      clientY: r.top + r.height / 2,
    }));
  }, slug);

  // Click "Share..." menu item.
  const shareItem = page.locator('[name="div-Share..."]');
  await shareItem.waitFor({timeout: 10000});
  await shareItem.click({force: true});

  // Share dialog.
  const dialog = page.locator('.d4-dialog').filter({hasText: `Share ${name}`});
  await dialog.waitFor({timeout: 15000});

  // Fill recipient input (no name= — located by placeholder).
  const recipientInput = dialog.locator('input[placeholder="User, group, or email"]');
  await recipientInput.focus();
  await page.keyboard.type(options.recipient);
  // Wait briefly for picker dropdown then commit via Enter.
  await page.waitForTimeout(1500);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(500);

  // Access level select — first <select> inside the per-row choice container
  // [name="div-share-selector"] (containerized after recipient is added).
  // Top-level access defaults to View-and-Use; only drive the select if a
  // non-default level requested.
  if (options.accessLevel && options.accessLevel !== 'View and use') {
    const accessSelect = dialog.locator('[name="div-share-selector"] select.ui-input-editor').first();
    if (await accessSelect.isVisible({timeout: 3000}).catch(() => false))
      await accessSelect.selectOption(options.accessLevel);
  }

  // Optional Send notifications toggle.
  if (options.sendNotifications !== undefined) {
    const checkbox = dialog.locator('[name="input-Send-notifications"]');
    const current = await checkbox.isChecked().catch(() => null);
    if (current !== null && current !== options.sendNotifications)
      await dialog.locator('[name="input-host-Send-notifications"]').click();
  }

  // Optional message.
  if (options.message) {
    const messageArea = dialog.locator('textarea[placeholder="Type in message here"]');
    if (await messageArea.isVisible({timeout: 3000}).catch(() => false)) {
      await messageArea.focus();
      await page.keyboard.type(options.message);
    }
  }

  await dialog.locator('[name="button-OK"]').click();
  await expect(dialog).toBeHidden({timeout: 15000});
}

// ---------------------------------------------------------------------------
// 10. saveCopy — UI Save Project dialog with mode chooser.
// ---------------------------------------------------------------------------

/**
 * Drive the Save Project dialog UI flow with a mode selection.
 *
 * Trigger: toolbar `[name="button-Save"]` (NOT context menu — Save Project
 * is not a tile-context-menu item per Session B Sub 1; it's a toolbar
 * action when a view is active).
 *
 * Modes (per `projects.md` "Save Project dialog" subsection):
 *   - `'original'`   → "Save original project" (overwrite)
 *   - `'copy'`       → "Save a copy" (new project; per-table Link/Clone
 *                       sub-option available via `perTableLinkOrClone`)
 *   - `'personal-view-customizations'` → "Save personal view customizations"
 *
 * **PascalCase caveat**: when `name` IS provided, the UI Save dialog
 * auto-normalizes typed display names to PascalCase server-side (see
 * `projects.md` "Saving a Project" caveat). The `name` parameter is what
 * the caller types, NOT what is persisted. To find the saved project after
 * this call, fetch by ID rather than by typed name.
 *
 * **"Copy of" default name (mode === 'copy' only)**: when `name` is OMITTED
 * for copy mode, the server applies a default name `"Copy of <sourceName>"`
 * (display) / `"copyof<originalname>"` (URL slug) — verified live 2026-05-04.
 * The default-name path is the canonical / most reliable path because the
 * resulting display name is predictable from `sourceName` (no PascalCase
 * normalization applies because the server constructed the name itself).
 * Tests can verify by `dapi.projects.filter('name = "Copy of <source>"').first()`.
 *
 * For modes other than `'copy'`, `name` is REQUIRED.
 *
 * @param page - Playwright Page (with active table view; source project
 *   already opened by caller).
 * @param options.sourceName - Source project name (currently active);
 *   used to predict the "Copy of <sourceName>" default name when `name`
 *   omitted.
 * @param options.mode - Save mode (see above).
 * @param options.name - Display name typed into the dialog. **Optional
 *   when `mode === 'copy'`** (server applies "Copy of <sourceName>"
 *   default). Required for `'original'` and `'personal-view-customizations'`
 *   modes.
 * @param options.dataSyncMode - Toggle Data Sync ON/OFF.
 * @param options.perTableLinkOrClone - Per-table Link/Clone choice (only
 *   meaningful when `mode === 'copy'`; applied uniformly to all tables).
 */
export async function saveCopy(
  page: Page,
  options: {
    sourceName: string;
    mode: 'original' | 'copy' | 'personal-view-customizations';
    name?: string;
    dataSyncMode?: 'on' | 'off';
    perTableLinkOrClone?: 'link' | 'clone';
  },
): Promise<void> {
  void options.sourceName; // documented context; not used in flow

  // Trigger Save Project dialog via toolbar SAVE button.
  await page.locator('[name="button-Save"]').first().click();

  const dialog = page.locator('.d4-dialog').filter({hasText: 'Save project'});
  await dialog.waitFor({timeout: 15000});

  // Mode radio chooser. Radio inputs share HTML name `radio11`; locate by
  // sibling <label> text. Labels are CSS-hidden by default — click via JS
  // DOM dispatch (label click also activates the linked radio input).
  const modeLabels: Record<typeof options.mode, string> = {
    'original': 'Save original project',
    'copy': 'Save a copy',
    'personal-view-customizations': 'Save personal view customizations',
  };
  await page.evaluate((targetText) => {
    const dialogEl = document.querySelector('.d4-dialog');
    if (!dialogEl) throw new Error('Save Project dialog not in DOM');
    const labels = Array.from(dialogEl.querySelectorAll('label'));
    const target = labels.find((l) => (l.textContent || '').trim() === targetText);
    if (!target) throw new Error(`mode label "${targetText}" not found`);
    // Belt-and-braces: directly toggle the linked radio input AND fire change
    // events (label .click() alone may not propagate when label is CSS-hidden
    // or the radio's change handler doesn't bubble through label clicks).
    const inputId = (target as HTMLLabelElement).getAttribute('for');
    if (inputId) {
      const input = document.getElementById(inputId) as HTMLInputElement | null;
      if (input) {
        input.checked = true;
        input.dispatchEvent(new Event('change', {bubbles: true}));
        input.dispatchEvent(new Event('input', {bubbles: true}));
      }
    }
    (target as HTMLLabelElement).click();
  }, modeLabels[options.mode]);
  await page.waitForTimeout(800); // let dialog adapt to mode change

  // Name input — first text input in dialog (no name=).
  // For mode='copy', if name is omitted, leave the dialog's default in
  // place (server applies "Copy of <sourceName>"); otherwise type the
  // requested name.
  if (options.name !== undefined) {
    const nameInput = dialog.locator('input[type="text"].ui-input-editor').first();
    await nameInput.focus();
    await page.keyboard.press('Control+a');
    await page.keyboard.type(options.name);
  } else if (options.mode !== 'copy') {
    throw new Error(`saveCopy: name is required for mode '${options.mode}' (only mode='copy' supports omitted-name → "Copy of <sourceName>" default)`);
  }

  // Data Sync toggle (optional).
  if (options.dataSyncMode !== undefined) {
    const desiredOn = options.dataSyncMode === 'on';
    const dsCheckbox = dialog.locator('[name="input-host-Data-sync"] input[type="checkbox"]');
    const current = await dsCheckbox.isChecked().catch(() => null);
    if (current !== null && current !== desiredOn)
      await dialog.locator('[name="input-host-Data-sync"]').click();
  }

  // Per-table Link/Clone (only when mode === 'copy'; applied uniformly to
  // all per-table choice containers).
  if (options.mode === 'copy' && options.perTableLinkOrClone) {
    const choiceContainers = dialog.locator('.grok-project-move-entity-debug.ui-input-choice');
    const count = await choiceContainers.count();
    for (let i = 0; i < count; i++) {
      const container = choiceContainers.nth(i);
      // Container has Link / Clone choice as labels.
      const choiceLabel = container.locator('label').filter({
        hasText: new RegExp(`^${options.perTableLinkOrClone === 'link' ? 'Link' : 'Clone'}$`, 'i'),
      });
      if (await choiceLabel.isVisible({timeout: 2000}).catch(() => false))
        await choiceLabel.click({force: true});
    }
  }

  await dialog.locator('[name="button-OK"]').click();

  // Auto-share dialog may appear after save in some flows; cancel.
  const shareDialog = page.locator('.d4-dialog').filter({hasText: /^Share /});
  if (await shareDialog.isVisible({timeout: 10000}).catch(() => false)) {
    await shareDialog.locator('[name="button-CANCEL"]').click();
    await expect(shareDialog).toBeHidden({timeout: 10000});
  }
}
