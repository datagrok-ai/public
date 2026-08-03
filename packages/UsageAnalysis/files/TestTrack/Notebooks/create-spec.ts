// Notebooks / Create — UI smoke covering Open in Notebook (Scenario 1), View as HTML
// (Scenario 2), and Delete via the Context Panel Actions (Scenario 3). Scenarios 4-5
// (Download formats, .ipynb drag-drop import) live in the create-ui.md manual companion.
//
// Known open bugs (both unfixed; asserted as today's invariant, broken behavior only observed):
//   GROK-16296 — Open in Notebook logs a console error during init (findProjectByView type cast).
//   GROK-13999 — clicking HTML throws a localStorage DOMException / 404 from a data: URL context.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// scenario: Create Notebook — Scenarios 1 (Open demog table in Notebook) + 2 (View notebook as HTML)
test('Notebooks / Create (UI Smoke): open demog in Notebook -> view as HTML', async ({page}) => {
  // 240s: cold `grok test` boots the Notebooks editor ribbon asynchronously (~10-13s post-click)
  // and the spec drives three sequential server round-trips (save, HTML-mode, delete).
  test.setTimeout(240_000);
  stepErrors.length = 0;

  let notebookId: string | null = null;

  await loginToDatagrok(page);
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });
  // Wait for closeAll to settle (no open table views) instead of a blind sleep.
  await page.waitForFunction(() => {
    try { return ((window as any).grok.shell.tableViews?.length ?? 0) === 0; }
    catch (e) { return false; }
  }, null, {timeout: 15_000}).catch(() => {});

  // ---- Setup: open demog.csv as the active table ----
  // The scenario opens via the demo-files panel; readCsv reaches the same TableView state the
  // Open-in-Notebook command consumes (it reads grok.shell.t). Open is JS-API setup, not an
  // owned UI flow, so this substitution is sanctioned for a ui-smoke spec.
  await softStep('Setup: open demog.csv and make it the active table', async () => {
    const info = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise<void>((res) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(); });
        setTimeout(res, 3000);
      });
      return {rows: df.rowCount, type: tv.type};
    });
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
    expect(info.rows).toBeGreaterThan(0);
  });

  // ---- Scenario 1: ML | Notebooks | Open in Notebook (UI top-menu) ----
  await softStep('Scenario 1: ML | Notebooks | Open in Notebook -> notebook view opens + persists', async () => {
    // Open the ML top menu, hover the Notebooks GROUP to expand its flyout, then click the leaf.
    // Live recon (public 2026-06-18): clicking the ML header builds the submenu DOM, but the
    // Open-in-Notebook leaf stays zero-size (offsetParent null, 0x0) — it does NOT become
    // clickable until the Notebooks group [name="div-ML---Notebooks"] is hovered, which lays out
    // the flyout and the leaf renders 160x24 / offsetParent set. (Behaviour identical with
    // simpleMode true/false; the prior raw-click-only opener was the cause of the 15s "resolved to
    // hidden" timeout.)
    await page.locator('[name="div-ML"]').first().click({timeout: 10_000});
    // Wait for the submenu DOM to build before hovering the group (was a blind 400ms sleep).
    const nbGroup = page.locator('[name="div-ML---Notebooks"]').first();
    await nbGroup.waitFor({state: 'attached', timeout: 8_000});
    await nbGroup.hover({timeout: 8_000});
    const leaf = page.locator('[name="div-ML---Notebooks---Open-in-Notebook"]').first();
    await leaf.waitFor({state: 'visible', timeout: 15_000});
    await leaf.click();

    // The command builds a template notebook from shell.t, saves it (notebooks.api.save), and
    // opens the Notebook editor view. Poll for the Notebook view to become current AND for the
    // notebook entity (grok.shell.o.id) to be populated. Live MCP recon (2026-06-17) showed the
    // view-type flips to 'Notebook' and o.id appears together at ~850ms warm, but on a cold
    // `grok test` run the Notebooks package boot can split those transitions — gating only on
    // v.type==='Notebook' lets Scenario 1 "complete" before the entity is current and races the
    // ribbon boot, which is the B-RUN-PASS/B-STAB-01 cold flake. Require BOTH before proceeding.
    const res = await page.waitForFunction(() => {
      const grok = (window as any).grok;
      let v; let o;
      try { v = grok.shell.v; o = grok.shell.o; } catch (e) { return false; } // getters can throw during boot
      if (v && v.type === 'Notebook' && o && o.id)
        return {viewName: v.name, oId: o.id};
      return false;
    }, null, {timeout: 90_000, polling: 250});
    const opened = await res.jsonValue() as {viewName: string; oId: string | null};

    // Invariant: a Notebook view opened (named after the source table — observed "Table_NN").
    expect(opened.viewName).toBeTruthy();
    const viewType = await page.evaluate(() => (window as any).grok.shell.v?.type);
    expect(viewType).toBe('Notebook');

    // Invariant: the notebook entity persisted to the server (notebooks.api.save). Resolve it
    // through dapi to confirm it is on the server, not just the in-memory current object. The
    // Open-in-Notebook command saves before returning (find() resolved on first try in warm
    // recon), but a cold server-commit lag is cheap to absorb with a short retry loop.
    notebookId = opened.oId;
    const persisted = await page.evaluate(async (id) => {
      const grok = (window as any).grok;
      if (!id) return {found: false};
      for (let i = 0; i < 10; i++) {
        try {
          const nb = await grok.dapi.notebooks.find(id);
          if (nb && nb.id) return {found: true};
        } catch (e) { /* retry */ }
        await new Promise((r) => setTimeout(r, 500));
      }
      return {found: false};
    }, notebookId);
    expect(persisted.found).toBe(true);

    // Behavioral remark (NOT a hard gate): GROK-16296 — Open in Notebook logs a console error
    // during init on dev (findProjectByView type cast). The view opens regardless. We surface
    // the observation but do not fail on it, because the bug is unfixed (regression-risk).
    const initErrors = await page.evaluate(() => {
      // grok.shell.lastError is not reliable across builds; this is a best-effort note only.
      try { return (window as any).grok?.shell?.lastError ?? null; } catch (e) { return null; }
    });
    if (initErrors)
      console.warn(`[OBSERVE GROK-16296] notebook-init surfaced an error (expected while unfixed): ${initErrors}`);
  });

  // ---- Scenario 2: View notebook as HTML (UI ribbon button) ----
  await softStep('Scenario 2: click HTML in the notebook ribbon -> HTML mode is reachable', async () => {
    // [name="button-HTML"] is the edit-mode ribbon button that switches to HTML rendering
    // (class-2 selector — see recon-notes header). Driving it exercises notebooks.editor.html-mode
    // + notebooks.editor.ribbon.html-mode. The ribbon mounts asynchronously AFTER the Notebook
    // view becomes current: live MCP recon (2026-06-17) measured the button absent until ~2.4s
    // after v.type==='Notebook' while the Notebooks package finishes booting the editor. On a
    // cold `grok test` run that lag is larger, so wait generously (the prior 30s/30s split was
    // the B-STAB-01 risk under cold boot). Wait attached then visible with a single 60s budget.
    const htmlBtn = page.locator('[name="button-HTML"]').first();
    await htmlBtn.waitFor({timeout: 60_000, state: 'attached'});
    await htmlBtn.waitFor({timeout: 60_000, state: 'visible'});
    await htmlBtn.click();

    // Invariant: the view stayed a Notebook view and the mode switch actually happened. Clicking
    // [name="button-HTML"] toggles the editor into HTML mode, and the ribbon button flips to its
    // inverse [name="button-EDIT"] (HTML-mode -> edit) — confirmed live on public 2026-06-18.
    // Asserting button-HTML is *still* present after the click is wrong (it has become button-EDIT);
    // asserting button-EDIT is the positive evidence the HTML mode switch was driveable. The
    // auto-retrying toBeVisible below waits for the switch — no fixed sleep needed. We do NOT
    // assert the rendered HTML content loaded: GROK-13999 (localStorage DOMException / 404 in the
    // data: URL render context) is unfixed, so the content area may stay empty.
    await expect(page.locator('[name="button-EDIT"]').first()).toBeVisible();
    const viewType = await page.evaluate(() => (window as any).grok.shell.v?.type);
    expect(viewType).toBe('Notebook');

    // Behavioral remark (NOT a hard gate): GROK-13999 — HTML render fails on dev with a 404 /
    // DOMException; the Download combo therefore never mounts. Recorded as an observation.
    const contentChildren = await page.evaluate(() => {
      const el = document.querySelector('[data-source="Notebooks:Notebook"]');
      return el ? el.children.length : -1;
    });
    if (contentChildren <= 0)
      console.warn('[OBSERVE GROK-13999] HTML-mode content did not render (expected while unfixed); ' +
        'Download (As HTML / As PDF) is consequently unavailable — covered manually in create-ui.md.');
  });

  // ---- Scenario 3: Delete the created notebook via the Context Panel Actions > Delete (UI) ----
  // Deletion coverage (notebooks.menu.delete / notebooks.meta.delete / notebooks.api.delete) is folded
  // into the create flow here: the standalone gallery context-menu delete is environmentally flaky on
  // public (the "find a freshly-created card in the gallery" round-trip races the server search index
  // and the gallery re-render). The just-created notebook is the CURRENT object (grok.shell.o), so its
  // Context Panel accordion exposes Actions > "Delete" directly — no gallery round-trip needed (live
  // recon public 2026-06-18: the Actions pane lists Open/Edit/Delete/Save As JSON/Share.../Rename... as
  // .d4-link-action links). Drive the UI Delete, confirm, and verify server-side removal via dapi.
  await softStep('Scenario 3: delete the notebook via Context Panel Actions > Delete -> removed from server', async () => {
    await page.evaluate(() => {
      const grok = (window as any).grok;
      try { grok.shell.windows.showContextPanel = true; } catch (e) { /* older builds */ }
    });
    // The accordion sections (incl. Actions) render for the current notebook object; the section toggle
    // is [name="div-section--Actions"] (there is no [name="pane-Actions"] on public — see notebooks.md).
    await page.locator('[name="div-section--Actions"]').first().waitFor({state: 'attached', timeout: 30_000});
    await page.locator('[name="div-section--Actions"]').first().click().catch(() => {});
    // Poll for the Delete link to render inside the Actions pane (was a blind 800ms sleep).
    await page.waitForFunction(() => {
      const sec = document.querySelector('[name="div-section--Actions"]');
      const pane = sec?.closest('.d4-accordion-pane');
      if (!pane) return false;
      return Array.from(pane.querySelectorAll('.d4-link-action'))
        .some((e) => (e.textContent ?? '').trim() === 'Delete');
    }, null, {timeout: 30_000}).catch(() => {});
    const clickedDelete = await page.evaluate(() => {
      const sec = document.querySelector('[name="div-section--Actions"]');
      const pane = sec?.closest('.d4-accordion-pane');
      const del = pane
        ? Array.from(pane.querySelectorAll('.d4-link-action')).find((e) => (e.textContent ?? '').trim() === 'Delete') as HTMLElement | undefined
        : undefined;
      if (del) { del.click(); return true; }
      return false;
    });
    expect(clickedDelete, 'Context Panel Actions > "Delete" link should be present and clickable').toBe(true);

    // Delete opens a confirm dialog (button-YES) — click it if present (best-effort: some builds delete
    // without a prompt). Removal is verified authoritatively via dapi below, so this is not hard-gated.
    await page.evaluate(async () => {
      for (let i = 0; i < 60; i++) {
        const dlg = Array.from(document.querySelectorAll('.d4-dialog'))
          .find((d) => /are you sure|delete notebook/i.test(d.textContent || ''));
        if (dlg) { const yes = dlg.querySelector('[name="button-YES"]') as HTMLElement | null; if (yes) { yes.click(); return; } }
        await new Promise((r) => setTimeout(r, 100));
      }
    });

    // Invariant: the notebook no longer resolves on the server (notebooks.api.delete). dapi.find is
    // authoritative and returns undefined once deleted; retry to absorb server-commit lag.
    const gone = await page.evaluate(async (id) => {
      const grok = (window as any).grok;
      if (!id) return false;
      for (let i = 0; i < 30; i++) {
        const ent: any = await grok.dapi.notebooks.find(id).catch(() => undefined);
        if (!ent || !ent.id) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    }, notebookId);
    expect(gone, 'the deleted notebook should no longer resolve on the server').toBe(true);
  });

  // ---- Cleanup: safety-net removal in case the UI delete above did not remove the entity ----
  // dapi.notebooks.delete requires the resolved entity object (passing the id string throws
  // "reading 'gaf'") — observed live 2026-06-17. Best-effort so the test never leaks a notebook.
  await page.evaluate(async (id) => {
    const grok = (window as any).grok;
    try {
      if (id) {
        const nb = await grok.dapi.notebooks.find(id);
        if (nb) await grok.dapi.notebooks.delete(nb);
      }
      grok.shell.closeAll();
    } catch (e) { /* best-effort cleanup */ }
  }, notebookId).catch(() => {});

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
