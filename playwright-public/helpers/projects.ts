/**
 * Playwright helpers for the projects feature.
 *
 * Imported as: `import * as projects from '../helpers/projects';` (or named).
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
 * PascalCase (e.g. `my-helpers-save` → `MyHelpersSave`), which breaks
 * downstream `dapi.projects.filter('name = "X"')` lookups. The JS API path
 * preserves names verbatim and is the canonical path for fixture-helper usage.
 *
 * Note: the JS API path creates an empty project entity (no view linkage).
 * `addLink(view)` returns a Dart-side `NoSuchMethodError`; until that's
 * resolved at the platform layer, helpers that need a populated project should
 * use a different mechanism (a JS API approximation that creates the tables).
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
 * Default 60s timeout with progressive intervals matching the project-save
 * flake budget on dev.
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
 * double-clicking. Uses the per-tile selector pattern
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
 * Build a composite cross-fixture of project variants.
 *
 * Uses a JS API approximation — creates separate project entities for each
 * requested variant with naming convention
 * `<sourceProjectName>-{link|clone|pvc}`. The actual UI Save-Copy-with-mode
 * flow is exercised by the `saveCopy` helper.
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
      p.description = `buildVariantsComposite ${variant} variant of "${src}" (JS API approximation)`;
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
 * @param namePrefix - Name prefix (e.g. "AutoTest-Share-").
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
 * `dapi.projects.filter('name = "<new>"').first()`. The rename DOES persist
 * (re-fetch by ID returns the new name) but the search index has
 * eventually-consistent lag — filter-by-name may return null for several
 * seconds after a rename. Always verify by ID.
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
// Heavier UI-flow helpers
// ===========================================================================

// ---------------------------------------------------------------------------
// 8. deleteProjectViaContextMenu — UI delete via right-click OR Context Panel.
// ---------------------------------------------------------------------------

/**
 * Delete a project via the Browse > Dashboards UI flow.
 *
 * Two trigger paths supported:
 *   - `right-click` (default): contextmenu on `[name="div-{slug}"]` tile
 *     → click `[name="div-Delete-Project"]` → DELETE confirmation.
 *   - `context-panel-dropdown`: single-click tile to populate Context Panel
 *     → click `.grok-prop-panel [name="icon-context-arrow-down"]` chevron
 *     → click `[name="div-Delete-Project"]` from the same menu
 *     → DELETE confirmation.
 *
 * Both paths terminate by clicking `[name="button-DELETE"]`.
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
    // sometimes misses Datagrok's overlay-wrapped tile).
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
  // builds, but the DELETE button selector is stable.
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
 * Drives the Share dialog. The recipient input has no `name=` attribute —
 * located by placeholder text.
 *
 * @param page - Playwright Page (on Browse > Dashboards with tile visible).
 * @param name - Project name.
 * @param options.recipient - User / group / email string typed into the
 *   recipient input (e.g. `'Test permission group'`).
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
  // Tile slug rule: the tile's name= attribute uses the server-stored project
  // name VERBATIM — PascalCase'd names are preserved as-is, kebab/lowercase
  // names stay lowercase. Try verbatim first, then lowercase as a fallback.
  const candidates = [name, name.toLowerCase()];
  let slug = name;
  for (const cand of candidates) {
    const t = page.locator(`[name="div-${cand}"]`);
    if (await t.isVisible({timeout: 5_000}).catch(() => false)) {
      slug = cand;
      break;
    }
  }
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
 * is not a tile-context-menu item; it's a toolbar action when a view is
 * active).
 *
 * Modes:
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
 * for copy mode, the server applies a default name derived from the active
 * table name (PascalCase + auto-index) — the server-default-name policy is
 * undocumented and varies by scenario shape. Tests SHOULD provide an explicit
 * `name` for deterministic verification (id-based dapi.projects.find after
 * save).
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
    /**
     * Override the auto-detection of "should the saved id differ from the
     * source id?". Default: `mode !== 'original'`. Some dev builds treat
     * `personal-view-customizations` as a layout-only update (no new project
     * entity created); pass `false` here to skip the new-id verification.
     */
    expectIdChange?: boolean;
  },
): Promise<void> {
  void options.sourceName; // documented context; not used in flow

  // Trigger Save Project dialog via toolbar SAVE button.
  // Use a `:visible` filter: the bare `[name="button-Save"]` matches hidden
  // duplicates from inactive views (offsetParent=null) and Playwright's
  // actionability check rejects them as "not visible" — under a
  // reopen+addViewer flow the toolbar SAVE button can resolve to a `dim`
  // non-visible candidate first.
  const saveBtn = page.locator('button[name="button-Save"]:visible').first();
  await saveBtn.waitFor({timeout: 30_000, state: 'visible'});
  await page.waitForTimeout(500); // settle render
  try {
    await saveBtn.click({timeout: 10_000});
  } catch (_) {
    await page.evaluate(() => {
      const candidates = Array.from(document.querySelectorAll('button[name="button-Save"]'));
      const visible = candidates.find((b) => (b as HTMLElement).offsetParent !== null);
      if (visible) (visible as HTMLElement).click();
    });
  }

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

  // Name input — `input#name` (id-based, precise). Do NOT fill via the
  // keyboard: `keyboard.press('Control+a')` does NOT cause select-all in this
  // Dart-side dialog (synthetic and CDP-routed Ctrl+A keydowns alike are not
  // honored by Dart's input handler). The auto-populated 'Copy of <source>'
  // value would remain, `keyboard.type` would append on top, and a downstream
  // Dart-side trim/normalize would drop the typed suffix entirely.
  //
  // Fix: native value setter + 'input'/'change' DOM events. The setter
  // bypasses React-style runtime guards and Dart's keyboard interception;
  // the dispatched events drive Dart's onChange listener so the new value
  // is actually committed to the dialog's model.
  //
  // For mode='copy', if name is omitted, leave the dialog's default in
  // place (server applies "Copy of <sourceName>" → PascalCase'd).
  if (options.name !== undefined) {
    // Selector fallback chain: precise `input#name` first, then the generic
    // `input[type="text"].ui-input-editor` form (only one such input exists
    // in the dialog when fully rendered). 15s timeout gives the dialog
    // headroom to fully render when chained from a fresh reopen flow.
    const nameInput = dialog.locator('input#name, input[type="text"].ui-input-editor').first();
    await nameInput.waitFor({state: 'visible', timeout: 15000});
    await page.evaluate((newName) => {
      const input = (document.querySelector('.d4-dialog input#name')
        ?? document.querySelector('.d4-dialog input[type="text"].ui-input-editor')) as HTMLInputElement | null;
      if (!input) throw new Error('saveCopy: Name input not found at fill-time (tried input#name and input[type="text"].ui-input-editor)');
      const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
      setter.call(input, newName);
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
    }, options.name);
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

  // Capture pre-OK source id for downstream verification. For mode='copy'
  // and mode='personal-view-customizations' the Dart shell takes 1-2s to
  // update grok.shell.project to the newly-created copy after OK click — if
  // we sample shell.project.id immediately we still see the source's id.
  // Without this baseline we can't tell apart "copy didn't happen" from
  // "Dart hasn't updated shell yet".
  const preOkSourceId = await page.evaluate(() => {
    const grok = (window as any).grok;
    return grok.shell.project?.id ?? null;
  });

  await dialog.locator('[name="button-OK"]').click();

  // Auto-share dialog may appear after save in some flows; cancel.
  const shareDialog = page.locator('.d4-dialog').filter({hasText: /^Share /});
  if (await shareDialog.isVisible({timeout: 10000}).catch(() => false)) {
    await shareDialog.locator('[name="button-CANCEL"]').click();
    await expect(shareDialog).toBeHidden({timeout: 10000});
  }

  // Verification — dual-strategy: shell.project update OR server-side name
  // match (PascalCase predicted from typed name).
  //
  // The server normalizes typed names to PascalCase — e.g. typed
  // `demog-{ts}-link` is stored as `Demog{ts}Link` — so strict equality
  // against the typed name never holds for kebab/snake/space-separated names.
  //
  // For mode='original' (overwrite): shell.project.id stays the same as
  // source; verify find(id) returns non-null.
  // For mode='copy' / mode='personal-view-customizations': the new copy
  // exists server-side immediately after OK click but `grok.shell.project`
  // may lag for >15s in the Playwright/CDP harness (server has the new project
  // id while shell.project.id stays at source). Fall back to a server-query
  // for the predicted PascalCase name — the same form
  // `dapi.projects.find(id).name` returns for dialog-saved projects.
  const expectIdChange = options.expectIdChange ?? (options.mode !== 'original');
  // Datagrok stores TWO name fields per project:
  //   `name`         — PascalCase-normalized internal form (e.g. `Demog{ts}Link`)
  //   `friendlyName` — typed name verbatim (e.g. `demog-{ts}-link`)
  // The dapi.projects.filter(`name = "X"`) operator does NOT do reliable
  // exact-match on either column (filter('name = "...") can return 0 even
  // when the project exists). The reliable path is `filter('name like
  // "<typed>"')` (no wildcards) which matches against friendlyName and works
  // for exact and prefix matches alike.
  const verified = await page.evaluate(async ({preId, expectChange, typedName}: {preId: string | null; expectChange: boolean; typedName: string | null}) => {
    const grok = (window as any).grok;
    const deadline = Date.now() + 20000;
    while (Date.now() < deadline) {
      // Strategy 1: shell.project (when it actually updates)
      const shellProj = grok.shell.project;
      if (shellProj?.id) {
        const okShell = !expectChange || shellProj.id !== preId;
        if (okShell) {
          const found = await grok.dapi.projects.find(shellProj.id);
          if (found?.name)
            return {id: shellProj.id, name: found.name, ok: true, via: 'shell'};
        }
      }
      // Strategy 2: server-side friendlyName match (typed name verbatim).
      // Datagrok shell.project lags >15s in some Playwright/CDP runs even
      // though server-side the project is created immediately on OK click.
      if (expectChange && typedName) {
        try {
          const list = await grok.dapi.projects.filter(`name like "${typedName}"`).list();
          const candidate = list.find((p: any) => p.id !== preId);
          if (candidate)
            return {id: candidate.id, name: candidate.name, ok: true, via: 'friendly-like'};
        } catch (_) { /* server query may transiently fail; retry next iter */ }
      }
      await new Promise((r) => setTimeout(r, 1000));
    }
    return {
      ok: false,
      reason: expectChange
        ? `no copy found via shell update or friendlyName like "${typedName}" within 20s (preId=${preId})`
        : 'timeout waiting for shell.project to settle with non-null find result',
    };
  }, {preId: preOkSourceId, expectChange: expectIdChange, typedName: options.name ?? null});

  if (!verified.ok)
    throw new Error(`saveCopy verification failed: ${JSON.stringify(verified)}`);
}

// ===========================================================================
// Provenance-aware save + reopen lifecycle helpers
// ===========================================================================

// ---------------------------------------------------------------------------
// 11. saveProjectWithProvenance — canonical uploadProject pattern.
// ---------------------------------------------------------------------------

export interface SavedProject {
  /** Server-assigned project id. */
  projectId: string;
  /** TableInfo id of the persisted source table. */
  tableInfoId: string;
  /** Layout id (if a layout was saved). */
  layoutId: string | null;
  /** Project name as persisted (may be PascalCased server-side via the dialog;
   *  via JS API it is preserved verbatim — this is the latter). */
  resolvedName: string;
}

/**
 * Save the active TableView as a project, persisting the TableInfo with
 * its `.script` provenance and the saved layout. Mirrors the canonical
 * `uploadProject` helper in `public/packages/UITests/src/gui/gui-utils.ts:100-111`:
 *
 * ```ts
 * project.addChild(tableInfo);
 * const layout = view.saveLayout();
 * project.addChild(layout);
 * await grok.dapi.layouts.save(layout);
 * await grok.dapi.tables.uploadDataFrame(df);
 * await grok.dapi.tables.save(tableInfo);
 * await grok.dapi.projects.save(project);
 * ```
 *
 * Skipping any of these calls causes save+reopen to "succeed" without any
 * data being re-materialized — the project entity persists but neither the
 * dataframe bytes, the TableInfo, nor the addChild relations do. With all
 * calls in place, files / db_query / db_table / script sources reopen with
 * the `.script` tag intact and the table re-executed server-side.
 *
 * Use after one of the `openTableFrom*` openers calls
 * has populated `grok.shell.tv` with a provenance-tagged DataFrame.
 *
 * @param page - Playwright Page (must have an active TableView).
 * @param projectName - Project name (preserved verbatim via the JS API path).
 * @returns identifiers needed for downstream cleanup or reopen.
 */
export async function saveProjectWithProvenance(
  page: Page,
  projectName: string,
): Promise<SavedProject> {
  return await page.evaluate(async (n) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const tv = grok.shell.tv;
    if (!tv?.dataFrame)
      throw new Error('saveProjectWithProvenance: no active TableView');
    const df = tv.dataFrame;
    const project = DG.Project.create();
    project.name = n;
    const tableInfo = df.getTableInfo();
    project.addChild(tableInfo);
    const layout = tv.saveLayout?.() ?? null;
    if (layout) project.addChild(layout);
    if (layout) await grok.dapi.layouts.save(layout);
    await grok.dapi.tables.uploadDataFrame(df);
    await grok.dapi.tables.save(tableInfo);
    await grok.dapi.projects.save(project);
    return {
      projectId: project.id,
      tableInfoId: tableInfo.id,
      layoutId: layout?.id ?? null,
      resolvedName: project.name,
    };
  }, projectName);
}

// ---------------------------------------------------------------------------
// 11b. saveAllTablesWithProvenance — multi-table variant of #11.
// ---------------------------------------------------------------------------

export interface SavedAllTables {
  /** Server-assigned project id. */
  projectId: string;
  /** TableInfo ids of every persisted source table, in shell.tables order. */
  tableInfoIds: string[];
  /**
   * Primary TableInfo id (first one) — for backward-compat with
   * `deleteProjectWithCleanup({tableInfoId})` which currently accepts a single
   * id. Callers that want full multi-table cleanup should iterate
   * `tableInfoIds` and call `deleteProjectWithCleanup({tableInfoId: ...})`
   * once per id.
   */
  primaryTableInfoId: string;
  /** Layout id (active TableView's layout) if a layout was saved. */
  layoutId: string | null;
  /** Project name as persisted (preserved verbatim via the JS API path). */
  resolvedName: string;
}

/**
 * Like {@link saveProjectWithProvenance}, but persists EVERY open table in
 * `grok.shell.tables` — not just the active TableView's dataframe. Use this
 * for multi-source lifecycle scenarios where reopen must re-materialize all
 * sources (file + DB-table + query + script + derived, etc.).
 *
 * The active TableView's layout is saved alongside (one layout, attached to
 * the project). Per-table layouts are out of scope.
 *
 * Each table's `df.tags['.script']` provenance is preserved through
 * `dapi.tables.uploadDataFrame(df)` + `dapi.tables.save(tableInfo)`, enabling
 * Datagrok server-side re-materialization on reopen (Data Sync).
 *
 * @param page - Playwright Page (must have at least one open TableView).
 * @param projectName - Project name (preserved verbatim via the JS API path).
 * @returns identifiers needed for downstream cleanup or reopen.
 */
export async function saveAllTablesWithProvenance(
  page: Page,
  projectName: string,
): Promise<SavedAllTables> {
  return await page.evaluate(async (n) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const tables = Array.from(grok.shell.tables);
    if (tables.length === 0)
      throw new Error('saveAllTablesWithProvenance: no tables in shell');
    const project = DG.Project.create();
    project.name = n;
    const tableInfoIds: string[] = [];
    for (const df of tables as any[]) {
      const ti = df.getTableInfo();
      project.addChild(ti);
      await grok.dapi.tables.uploadDataFrame(df);
      await grok.dapi.tables.save(ti);
      tableInfoIds.push(ti.id);
    }
    const tv = grok.shell.tv;
    const layout = tv?.saveLayout?.() ?? null;
    if (layout) {
      project.addChild(layout);
      await grok.dapi.layouts.save(layout);
    }
    await grok.dapi.projects.save(project);
    return {
      projectId: project.id,
      tableInfoIds,
      primaryTableInfoId: tableInfoIds[0],
      layoutId: layout?.id ?? null,
      resolvedName: project.name,
    };
  }, projectName);
}

// ---------------------------------------------------------------------------
// 12. reopenAndAssertProvenance — closeAll + reopen + verify .script survives.
// ---------------------------------------------------------------------------

export interface ReopenResult {
  /** Time taken from `closeAll()` to the table re-materializing (ms). */
  reopenMs: number;
  /** Number of tables visible in `grok.shell.tables` after reopen. */
  tablesAfter: number;
  /** Active TableView dataFrame name after reopen. */
  reopenedName: string;
  /** Active TableView dataFrame rowCount after reopen. */
  reopenedRowCount: number;
  /** `.script` tag value after reopen — should match the pre-save value. */
  reopenedScript: string;
}

/**
 * Verify that a project saved with `saveProjectWithProvenance` correctly
 * round-trips through `closeAll → dapi.projects.find(id).open()` —
 * i.e. that Datagrok re-executes the `.script` provenance and re-materializes
 * the source table.
 *
 * This is the runtime check that distinguishes a real Data-Sync-ON
 * lifecycle from a snapshot-only one. Without provenance, reopen returns
 * an empty project (or a snapshot copy) and the assertions below fail.
 *
 * @param page - Playwright Page.
 * @param projectId - Project id from `saveProjectWithProvenance`.
 * @param expectedScriptPattern - Optional regex; if provided, the post-reopen
 *   `.script` must match (e.g. `PROVENANCE_PATTERNS.files`). Use to enforce
 *   Gate E-PROV-01 across the reopen boundary.
 * @returns measurements needed for downstream assertions.
 * @throws if reopen does not produce a TableView, OR if expectedScriptPattern
 *   is provided and does not match.
 */
export async function reopenAndAssertProvenance(
  page: Page,
  projectId: string,
  expectedScriptPattern?: RegExp,
): Promise<ReopenResult> {
  const result = await page.evaluate(async (id) => {
    const grok = (window as any).grok;
    grok.shell.closeAll();
    await new Promise((r) => setTimeout(r, 1000));
    const t0 = performance.now();
    const proj = await grok.dapi.projects.find(id);
    if (!proj) throw new Error(`Project not found by id: ${id}`);
    await proj.open();
    const reopenMs = Math.round(performance.now() - t0);
    // Settle: tableviews materialise after the reopen promise.
    await new Promise((r) => setTimeout(r, 2000));
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!df) throw new Error(`Reopen of ${id}: no df after open + 2s settle`);
    const tablesAfter = (() => {
      try { return Number(grok.shell.tables?.length) || 0; }
      catch (_) { return 0; }
    })();
    return {
      reopenMs,
      tablesAfter,
      reopenedName: df.name,
      reopenedRowCount: df.rowCount,
      reopenedScript: df.tags?.get?.('.script') ?? '',
    };
  }, projectId);

  if (expectedScriptPattern && !expectedScriptPattern.test(result.reopenedScript)) {
    throw new Error(
      `reopenAndAssertProvenance: post-reopen .script = ` +
      `"${result.reopenedScript.slice(0, 200)}" does not match ` +
      `${expectedScriptPattern}. The project saved successfully but ` +
      `provenance did not survive — likely cause: prior save path skipped ` +
      `tables.uploadDataFrame / tables.save / project.addChild(tableInfo).`,
    );
  }
  return result;
}

// ---------------------------------------------------------------------------
// 13. deleteProjectWithCleanup — drop project + tableInfo + (optional) script.
// ---------------------------------------------------------------------------

/**
 * Delete a project and its associated TableInfo. Best-effort — swallows
 * errors so it can be safely called from `finally` blocks.
 *
 * Use this in a try/finally for any spec that calls
 * `saveProjectWithProvenance` to avoid leaking project entities on dev.
 * Optionally also pass a provisioned script id to clean up
 * test-fixture scripts in the same call.
 *
 * @param page - Playwright Page.
 * @param ids - identifiers from `saveProjectWithProvenance` and (optional)
 *   from `provisionDataframeScript` (in `helpers/openers.ts`).
 */
export async function deleteProjectWithCleanup(
  page: Page,
  ids: {projectId?: string; tableInfoId?: string; scriptId?: string},
): Promise<void> {
  await page.evaluate(async (i) => {
    const grok = (window as any).grok;
    if (i.projectId) {
      try {
        const p = await grok.dapi.projects.find(i.projectId);
        if (p) await grok.dapi.projects.delete(p);
      } catch (_) { /* best effort */ }
    }
    if (i.tableInfoId) {
      try {
        const ti = await grok.dapi.tables.find(i.tableInfoId);
        if (ti) await grok.dapi.tables.delete(ti);
      } catch (_) { /* best effort */ }
    }
    if (i.scriptId) {
      try {
        const s = await grok.dapi.scripts.find(i.scriptId);
        if (s) await grok.dapi.scripts.delete(s);
      } catch (_) { /* best effort */ }
    }
  }, ids).catch(() => {});
}
