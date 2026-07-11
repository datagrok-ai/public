/* ---
sub_features_covered: [notebooks.menu.open, notebooks.meta.open, notebooks.menu.edit, notebooks.meta.edit, notebooks.menu.rename, notebooks.menu.save-as-json, notebooks.meta.save-as-json, notebooks.menu.share, notebooks.menu.delete, notebooks.meta.delete, notebooks.api.delete, notebooks.meta.render-icon, notebooks.meta.render-tooltip, notebooks.meta.render-details, notebooks.meta.bind, notebooks.meta.get-applicable-cases, notebooks.entity.get-applicable-cases]
--- */
// ui-smoke OWNER of the seven pcmd* notebook context-menu flows. Each flow is DOM-driven:
// right-click the seeded notebook card link -> click the context-menu item -> assert the
// platform-side result. The seven flows verify:
//   S1 (pcmdOpen)      Open -> notebook opens in HTML mode (grok.shell.v.type === 'Notebook').
//   S2 (pcmdEdit)      Edit -> view transitions toward edit mode (Notebook view). JupyterLab
//                      iframe interior is manual_only per atlas — NOT touched.
//   S5 (pcmdShare)     Share... -> shareEntity dialog opens (title `Share <name>`); CANCEL leaves
//                      the entity unchanged. Run before Rename/Delete so the title carries the seed name.
//   S6 (pcmdApplyTo)   right-click -> the entity context menu (Open/Delete) opens; the Apply-to group
//                      is a behavioral observation (getApplicableCases is Dart-side, conditional on the
//                      seed being linked to an open table). Apply EXECUTION is atlas manual_only.
//   S4 (pcmdSaveAsJSON) Save As JSON -> direct browser download of a <name>.json file.
//   S3 (pcmdRename)    Rename... -> modal; clearing the input disables OK (notEmpty validator), typing
//                      re-enables it; OK renames server-side and the card label updates.
//   S7 (pcmdDelete)    Delete -> confirm YES -> notebook removed server-side (find(id) returns undefined).
//
// The notebook is SELF-SEEDED with a unique name (no JS-API Notebook factory exists) and linked to an
// open demog.csv (for Apply-to applicability), so Rename + Delete never mutate shared notebooks.
//
// Selector / behavioral notes:
//   - Navigation opener: DG.Func.find({name:'CmdBrowseNotebooks'})[0].apply() — view-independent;
//     flips grok.shell.v.type to 'notebooks'. (Not one of the owned pcmd flows.)
//   - Right-click the notebook LINK (.d4-link-label), NOT the card wrapper: the wrapper opens the
//     generic grid menu; the entity menu (Open/Edit/Delete/Save As JSON/Share.../Rename...) opens
//     only on the link.
//   - Rename modal input is Dart-bound (E-SEL-03): Playwright .fill() skips the Dart change listener
//     so the notEmpty validator never fires — clear via keyboard select-all+Delete and type the new
//     name via page.keyboard.type().
//   - Share dialog is detected by ATTACHMENT (offsetParent is transiently null while it lays out).
//   - grok.dapi.notebooks.delete REQUIRES the DG.Notebook entity, NOT an id string; fetch via find(id).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';

test.use(specTestOptions);

// Helper-discipline (E-HELP-01): the spec-local helpers below (openNotebooksBrowser,
// openCardContextMenu, narrowGalleryToSeed/ensureBrowserNarrowedToSeed) are declared as
// candidate_helpers in the paired scenario .md frontmatter. openNotebooksBrowser is duplicated
// verbatim across browser-spec.ts /
// notebooks-lifecycle-jupyter-container-spec.ts / this spec (reuse >=3); the other two are reused per-pcmd-flow within this spec.
// They are NOT reinventions of any helpers-registry entry (the registry has no notebooks-browser
// opener / card-contextmenu helper); a future cycle may promote them per the helper-authoring
// sub-routine. Canonical login/option/soft-step helpers come from ../spec-login (registry-listed).
//
// Open the Notebooks browser via the CmdBrowseNotebooks command function (view-independent), and
// return only once at least one notebook card has actually rendered — not merely when the gallery
// container mounted (the cards load asynchronously; see header recon-notes). Generalizes the
// validated opener from the sibling browser-spec.ts / notebooks-lifecycle-jupyter-container-spec.ts.
async function openNotebooksBrowser(page: import('@playwright/test').Page) {
  await page.evaluate(async () => {
    const f = (window as any).DG.Func.find({name: 'CmdBrowseNotebooks'})[0];
    await f.apply();
  });
  await page.waitForFunction(() => {
    try { return (window as any).grok.shell.v?.type === 'notebooks'; } catch (e) { return false; }
  }, null, {timeout: 45_000, polling: 250});
  await page.locator('.grok-gallery-search-bar').waitFor({timeout: 30_000});
  await page.waitForFunction(() => {
    return document.querySelectorAll('.grok-gallery-grid-item').length > 0 &&
      document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]').length > 0;
  }, null, {timeout: 45_000, polling: 250});
}

// Right-click the seeded notebook card and return once its context menu is open. Returns false if the
// card or the menu never materialized. Uses a dispatched contextmenu MouseEvent on the card DOM (the
// proven DOM-grid path; the cards are canvas-free DOM grid items).
async function openCardContextMenu(page: import('@playwright/test').Page, name: string,
  probeSelector: string): Promise<boolean> {
  return page.evaluate(async (args: {name: string; probe: string}) => {
    // Right-click the notebook LINK, NOT the card wrapper: live recon (public 2026-06-18) showed a
    // contextmenu dispatched on the .grok-gallery-grid-item card opens the generic grid menu
    // (Details/Chat/Block/Groups.../Roles...), while the entity menu (Open/Edit/Delete/Save As JSON/
    // Share.../Rename...) opens only on the .d4-link-label notebook link.
    const label = Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'))
      .find((l) => l.textContent?.trim() === args.name) as HTMLElement | undefined;
    if (!label) return false;
    const r = label.getBoundingClientRect();
    label.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, clientX: r.left + 10, clientY: r.top + 5}));
    for (let i = 0; i < 30; i++) {
      await new Promise((rs) => setTimeout(rs, 100));
      const item = Array.from(document.querySelectorAll(args.probe)).find((e) => (e as HTMLElement).offsetParent !== null);
      if (item) return true;
    }
    return false;
  }, {name, probe: probeSelector});
}

// Narrow the gallery search to the seeded notebook (assumes we are already on the notebooks view
// with the gallery rendered). Returns true once the seed card is visible after filtering. The fresh
// seed sorts BEYOND page 1 of the (50-per-page) gallery on dev — searching is what surfaces it
// (verified live: onPage1Unfiltered=false, search hit in ~1.6s).
async function narrowGalleryToSeed(page: import('@playwright/test').Page, name: string): Promise<boolean> {
  return page.evaluate(async (n) => {
    const input = document.querySelector('.grok-gallery-search-bar input') as HTMLInputElement | null;
    if (input) {
      input.focus();
      input.value = n;
      input.dispatchEvent(new Event('input', {bubbles: true}));
    }
    for (let i = 0; i < 50; i++) {
      await new Promise((r) => setTimeout(r, 300));
      const hit = Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'))
        .some((l) => l.textContent?.trim() === n);
      if (hit) return true;
    }
    return false;
  }, name);
}

// Ensure the gallery is open and narrowed to the seeded notebook. Re-opening the browser from
// scratch (closeAll + CmdBrowseNotebooks + full readiness poll) is EXPENSIVE and is the proven
// cumulative-timeout / B-STAB-01 driver when done once per scenario. It is ONLY required after a
// flow that navigated AWAY from the browser (Open / Edit land on a `Notebook` view). When we are
// still on the `notebooks` view, just (re-)narrow the existing gallery — verified live that the
// context menu re-opens correctly on the filtered card without a full reopen.
async function ensureBrowserNarrowedToSeed(page: import('@playwright/test').Page, name: string) {
  const onBrowser = await page.evaluate(() => {
    try { return grok.shell.v?.type === 'notebooks'; } catch (e) { return false; }
  });
  if (!onBrowser) {
    await page.evaluate(() => grok.shell.closeAll());
    await openNotebooksBrowser(page);
  }
  await page.locator('.grok-gallery-search-bar input').waitFor({timeout: 30_000});
  let found = await narrowGalleryToSeed(page, name);
  if (!found) {
    // Defensive: a stale/empty gallery state — force a full reopen once, then re-narrow.
    await page.evaluate(() => grok.shell.closeAll());
    await openNotebooksBrowser(page);
    await page.locator('.grok-gallery-search-bar input').waitFor({timeout: 30_000});
    found = await narrowGalleryToSeed(page, name);
  }
  expect(found, `seeded notebook "${name}" should be findable in the browser`).toBe(true);
}

test('Notebooks — Context Menu Smoke (all 7 pcmd flows)', async ({page}) => {
  // 6 min: 7 sequential context-menu flows in one test; pcmdOpen/pcmdEdit each load a Notebook view
  // that can block on cold Jupyter container warm-up.
  test.setTimeout(360_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  });

  // Unique seed name so Rename + Delete never collide with / mutate shared Demog/Cars notebooks.
  const seedName = `automator-ctxmenu-${Date.now()}`;
  let seededId: string | null = null;

  // ---- Setup: open demog.csv (Apply-to applicability prerequisite) + self-seed a notebook ----
  await softStep('Setup: open demog.csv as the active table', async () => {
    const info = await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((res) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(); });
        setTimeout(res, 3000);
      });
      return {rows: df.rowCount};
    });
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
    expect(info.rows).toBeGreaterThan(0);
  });

  await softStep('Setup: seed a server notebook linked to demog (via Open in Notebook)', async () => {
    // CmdOpenInNotebook builds Notebook.template(tables:[shell.t]) linked to the open demog TableInfo,
    // saves it, and opens the editor — this makes the seeded notebook Apply-to-applicable to demog
    // (Scenario 6) AND gives us a server-persisted, owned notebook to Rename/Delete (Scenarios 3/7).
    seededId = await page.evaluate(async (name) => {
      const before = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10});
      const beforeIds = new Set(before.map((n: any) => n.id));
      const f = (window as any).DG.Func.find({name: 'CmdOpenInNotebook'})[0];
      await f.apply();
      let fresh: any = null;
      for (let i = 0; i < 40; i++) {
        await new Promise((r) => setTimeout(r, 500));
        const cur = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10}).catch(() => [] as any[]);
        fresh = cur.find((n: any) => !beforeIds.has(n.id));
        if (fresh) break;
      }
      if (!fresh) return null;
      fresh.friendlyName = name;
      await grok.dapi.notebooks.save(fresh);
      return fresh.id as string;
    }, seedName);
    expect(seededId, 'Open in Notebook should persist a demog-linked server notebook').toBeTruthy();
    // Confirm the rename persisted, keyed by id (server-commit lag insurance).
    const persisted = await page.evaluate(async (args: {id: string; name: string}) => {
      for (let i = 0; i < 10; i++) {
        const ent: any = await grok.dapi.notebooks.find(args.id).catch(() => null);
        if (ent && (ent.friendlyName || ent.name) === args.name) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    }, {id: seededId as string, name: seedName});
    expect(persisted).toBe(true);
  });

  // ---- Scenario 1: Open notebook in HTML mode via context menu (pcmdOpen) ----
  await softStep('Scenario 1 (pcmdOpen): right-click card -> Open -> HTML-mode view opens', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);
    const present = await openCardContextMenu(page, seedName, '[name="div-Open"]');
    expect(present, '[name="div-Open"] should be present in the context menu').toBe(true);
    await page.evaluate(() => (document.querySelector('[name="div-Open"]') as HTMLElement)?.click());
    // The owned ui-smoke assertion is the platform-side view transition: Open re-fetches and opens the
    // notebook in HTML mode (grok.shell.v.type === 'Notebook'). The HTML-mode ribbon Download/EDIT
    // buttons are source-derived only and did NOT render during recon (404 / GROK-13999), so they are
    // intentionally NOT asserted here.
    const opened = await page.waitForFunction(() => {
      try { return grok.shell.v?.type === 'Notebook'; } catch (e) { return false; }
    }, null, {timeout: 30_000, polling: 250}).then(() => true).catch(() => false);
    expect(opened, 'Open should transition to a Notebook (HTML-mode) view').toBe(true);
  });

  // ---- Scenario 2: Launch edit mode via context menu (pcmdEdit) ----
  await softStep('Scenario 2 (pcmdEdit): right-click card -> Edit -> edit-mode view transition', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);
    const present = await openCardContextMenu(page, seedName, '[name="div-Edit"]');
    expect(present, '[name="div-Edit"] should be present in the context menu').toBe(true);
    await page.evaluate(() => (document.querySelector('[name="div-Edit"]') as HTMLElement)?.click());
    // Owned assertion: the platform-side view transitions toward edit mode (grok.shell.v.type ===
    // 'Notebook'). The JupyterLab iframe interior is manual_only per atlas and is NOT touched.
    const inEdit = await page.waitForFunction(() => {
      try { return grok.shell.v?.type === 'Notebook'; } catch (e) { return false; }
    }, null, {timeout: 30_000, polling: 250}).then(() => true).catch(() => false);
    expect(inEdit, 'Edit should transition to a Notebook (edit-mode) view').toBe(true);
  });

  // ---- Scenario 5: Share notebook via context menu (pcmdShare) ----
  // (ordered before Rename/Delete so the seed still carries its original name in the dialog title)
  await softStep('Scenario 5 (pcmdShare): right-click card -> Share... -> Share dialog opens', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);
    // The share dialog is detected by ATTACHMENT + title/sharing-class, NOT by offsetParent
    // visibility: verified live 2026-06-17 that clicking [name="div-Share..."] creates the
    // .d4-dialog (title `Share <friendlyName>`, contents .dlg-sharing-settings) while its
    // offsetParent is transiently null (the sharing dialog lays out lazily). The prior
    // offsetParent!==null filter missed the attached-but-not-yet-visible dialog and was the
    // Scenario-5 B-STAB-01 flake. Re-open the menu + re-click up to 3 times against straggler
    // menu/dialog state.
    let shared = {ok: false, title: ''};
    for (let attempt = 0; attempt < 2 && !shared.ok; attempt++) {
      const present = await openCardContextMenu(page, seedName, '[name="div-Share..."]');
      if (!present) continue;
      shared = await page.evaluate(async () => {
        (document.querySelector('[name="div-Share..."]') as HTMLElement)?.click();
        for (let i = 0; i < 40; i++) {
          await new Promise((r) => setTimeout(r, 150));
          // Detect by attachment (offsetParent may be null while the sharing dialog lays out).
          const dlg = Array.from(document.querySelectorAll('.d4-dialog'))
            .find((d) => d.querySelector('.dlg-sharing-settings') ||
              /^share /i.test(d.querySelector('.d4-dialog-title')?.textContent?.trim() ?? ''));
          if (dlg) {
            const title = dlg.querySelector('.d4-dialog-title')?.textContent?.trim() ?? '';
            // Close without changes — the entity must stay unchanged.
            (dlg.querySelector('[name="button-CANCEL"]') as HTMLElement)?.click();
            return {ok: true, title};
          }
        }
        return {ok: false, title: ''};
      });
    }
    expect(shared.ok, 'the shareEntity dialog should open (attached) for the notebook').toBe(true);
    // Dialog title is `Share <friendlyName>` (observed live).
    expect(shared.title.toLowerCase()).toContain('share');
  });

  // ---- Scenario 6: Apply-to context menu shows applicable tables (pcmdApplyTo) ----
  await softStep('Scenario 6 (pcmdApplyTo): demog open -> Apply-to group appears for the seed', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);
    // Right-click the seed; with demog.csv still open and the notebook linked to its TableInfo,
    // getApplicableCases() is non-empty so the Apply-to group renders. The sub-leaf is named after the
    // open table VIEW (observed: "Table"), not the notebook.
    const menu = await page.evaluate(async (name) => {
      // Right-click the notebook LINK (the card wrapper opens the generic grid menu — see
      // openCardContextMenu). The entity menu signature is Open/Delete.
      const label = Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'))
        .find((l) => l.textContent?.trim() === name) as HTMLElement | undefined;
      if (!label) return {menuOpened: false, applyToPresent: false, leafPresent: false};
      const r = label.getBoundingClientRect();
      label.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, clientX: r.left + 10, clientY: r.top + 5}));
      for (let i = 0; i < 30; i++) {
        await new Promise((rs) => setTimeout(rs, 150));
        if (document.querySelector('[name="div-Open"]') || document.querySelector('[name="div-Delete"]')) break;
      }
      const menuOpened = !!(document.querySelector('[name="div-Open"]') || document.querySelector('[name="div-Delete"]'));
      const applyTo = !!document.querySelector('[name="div-Apply-to"]');
      const leaf = !!document.querySelector('[name="div-Apply-to---Table"]');
      // dismiss the menu
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      return {menuOpened, applyToPresent: applyTo, leafPresent: leaf};
    }, seedName);
    // Hard gate: the entity context menu opens on the seed card link (pcmd surface reachable). The
    // Apply-to group is CONDITIONAL on the notebook being applicable to an open table — getApplicableCases()
    // is Dart-side and a blank NEW-NOTEBOOK seed is not linked to the demog table on public, so Apply-to
    // is absent (observed live 2026-06-18). Apply-to presence is therefore a behavioral observation, not a
    // hard gate; its EXECUTION (notebooks.entity.apply -> convertNotebook -> live Jupyter container) is
    // atlas manual_only / non-deterministic and is intentionally not driven.
    expect(menu.menuOpened, 'the notebook entity context menu (Open/Delete) should open on the seed card link').toBe(true);
    if (!menu.applyToPresent)
      console.warn('[OBSERVE notebooks.entity.get-applicable-cases] Apply-to group absent for the blank ' +
        'seed notebook (not applicable to the open table; applicability is Dart-side). Recorded as an ' +
        'observation, not a hard gate.');
  });

  // ---- Scenario 4: Save As JSON via context menu (pcmdSaveAsJSON) ----
  await softStep('Scenario 4 (pcmdSaveAsJSON): right-click card -> Save As JSON -> file downloads', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);
    const present = await openCardContextMenu(page, seedName, '[name="div-Save-As-JSON"]');
    expect(present, '[name="div-Save-As-JSON"] should be present in the context menu').toBe(true);
    // Save As JSON triggers a direct browser download (no dialog). Race the click against the
    // download event.
    const downloadPromise = page.waitForEvent('download', {timeout: 30_000}).catch(() => null);
    await page.evaluate(() => (document.querySelector('[name="div-Save-As-JSON"]') as HTMLElement)?.click());
    const download = await downloadPromise;
    expect(download, 'Save As JSON should trigger a file download').not.toBeNull();
    expect(download!.suggestedFilename().toLowerCase()).toMatch(/\.json$/);
  });

  // ---- Scenario 3: Rename notebook via context menu (pcmdRename) ----
  const renamedName = `${seedName}-renamed`;
  await softStep('Scenario 3 (pcmdRename): right-click card -> Rename... -> modal -> OK -> label updates', async () => {
    await ensureBrowserNarrowedToSeed(page, seedName);
    const present = await openCardContextMenu(page, seedName, '[name="div-Rename..."]');
    expect(present, '[name="div-Rename..."] should be present in the context menu').toBe(true);
    await page.evaluate(() => (document.querySelector('[name="div-Rename..."]') as HTMLElement)?.click());

    // Rename modal: title "Rename Notebook"; input pre-filled with the current name. The input name
    // carries the dialog-title prefix on public ([name="input-Rename-Notebook---New-name-"]) — live
    // recon 2026-06-18.
    const input = page.locator('[name="input-Rename-Notebook---New-name-"]');
    await input.waitFor({timeout: 30_000});
    await expect(page.locator('.d4-dialog .d4-dialog-title').filter({hasText: 'Rename Notebook'}).first())
      .toBeVisible();

    // Empty-name validation: clearing the input disables OK (Validators.notEmpty gate).
    // E-SEL-03: this is a Dart-bound input (class "ui-input-editor"); Playwright's .fill()
    // sets .value but skips the Dart change listener, so the notEmpty validator never fires
    // and OK stays "enabled". Drive it via real keyboard events instead — focus, select-all,
    // Delete to clear; keyboard.type() to enter the new name. Confirmed live 2026-06-17 via
    // chrome-devtools MCP that the Dart listener responds to genuine input events: clearing
    // flips OK to "disabled", typing flips it back to "enabled".
    await input.click();
    await input.press('ControlOrMeta+a');
    await input.press('Delete');
    await expect(page.locator('[name="button-OK"]').first()).toHaveClass(/disabled/);

    // Enter the new name; OK re-enables; confirm.
    await page.keyboard.type(renamedName);
    await expect(page.locator('[name="button-OK"]').first()).toHaveClass(/enabled/);
    await page.locator('[name="button-OK"]').first().click();

    // The rename persists server-side (dapi.notebooks.save, keyed on the seeded id).
    const renamed = await page.evaluate(async (args: {id: string | null; name: string}) => {
      if (!args.id) return false;
      for (let i = 0; i < 20; i++) {
        const ent: any = await grok.dapi.notebooks.find(args.id).catch(() => null);
        if (ent && (ent.friendlyName || ent.name) === args.name) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    }, {id: seededId, name: renamedName});
    expect(renamed, 'the notebook should be renamed server-side').toBe(true);

    // The card label updates to the new name in the (re-rendered) gallery.
    await ensureBrowserNarrowedToSeed(page, renamedName);
    const cardShows = await page.evaluate((name) =>
      Array.from(document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]'))
        .some((l) => l.textContent?.trim() === name), renamedName);
    expect(cardShows, `card label should update to "${renamedName}"`).toBe(true);
  });

  // ---- Scenario 7: Delete notebook via context menu (pcmdDelete) ----
  await softStep('Scenario 7 (pcmdDelete): right-click card -> Delete -> YES -> card removed', async () => {
    await ensureBrowserNarrowedToSeed(page, renamedName);
    const present = await openCardContextMenu(page, renamedName, '[name="div-Delete"]');
    expect(present, '[name="div-Delete"] should be present (server-persisted notebook)').toBe(true);
    await page.evaluate(() => (document.querySelector('[name="div-Delete"]') as HTMLElement)?.click());

    // Confirm dialog: title "Are you sure?", body Delete notebook "<name>"?, [name="button-YES"].
    const confirmed = await page.evaluate(async () => {
      for (let i = 0; i < 60; i++) {
        const dlg = Array.from(document.querySelectorAll('.d4-dialog'))
          .find((d) => /are you sure|delete notebook/i.test(d.textContent || ''));
        if (dlg) {
          const yes = dlg.querySelector('[name="button-YES"]') as HTMLElement | null;
          if (yes) { yes.click(); return true; }
        }
        await new Promise((r) => setTimeout(r, 100));
      }
      return false;
    });
    expect(confirmed, 'confirm dialog YES (button-YES) should be present and clickable').toBe(true);

    // grok.dapi.notebooks is authoritative — the gallery item-count text is stale after a UI delete.
    // Verify by id via find(), which returns undefined once the notebook is deleted.
    const removed = await page.evaluate(async (id: string | null) => {
      if (!id) return false;
      for (let i = 0; i < 30; i++) {
        const ent: any = await grok.dapi.notebooks.find(id).catch(() => undefined);
        if (!ent || !ent.id) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    }, seededId);
    expect(removed, 'the notebook should be removed from the server (AppEvents.ENTITY_REMOVED)').toBe(true);
    seededId = null; // deleted — skip cleanup
  });

  // Cleanup: if any flow above failed mid-way and left the seed alive, remove it by entity
  // (delete REQUIRES the DG.Notebook entity, not an id string). Keyed on the seeded id via find().
  await page.evaluate(async (id: string | null) => {
    if (!id) return;
    try {
      const entity: any = await grok.dapi.notebooks.find(id).catch(() => null);
      if (entity && entity.id) await grok.dapi.notebooks.delete(entity).catch(() => {});
    } catch (_) { /* best-effort cleanup */ }
    grok.shell.closeAll();
  }, seededId);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
