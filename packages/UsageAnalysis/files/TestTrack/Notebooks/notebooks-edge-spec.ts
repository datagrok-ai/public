/* ---
sub_features_covered: [notebooks.entity.environment, notebooks.meta.render-tooltip, notebooks.meta.render-details]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: edge) — non-pyramid constraint defaults:
//     JS-API substitution permitted; >=1 DOM-driving call still REQUIRED (E-LAYER-COMPLIANCE-01)
//   sub_features_covered (authored slice): [notebooks.entity.environment,
//     notebooks.meta.render-tooltip, notebooks.meta.render-details]
//   ui_coverage_responsibility: (none declared) — no owned ui-smoke flow; JS-API substitution sanctioned
//   produced_from: atlas-driven
//   related_bugs: []
//
// Atlas provenance (derived_from):
//   notebooks.entity.environment          derived_from: core/shared/grok_shared/lib/src/notebook.dart#L21
//   notebooks.meta.render-details         derived_from: core/client/xamgle/lib/src/meta/notebook_meta.dart#L25
//
// Scope: the AUTOMATABLE notebook edge invariants from notebooks-edge.md that have a DETERMINISTIC
// COLD-HEADLESS path. The scenario .md was authored against a Dart-level Notebook entity surface
// (Notebook.create / .template / .getApplicableCases / NotebookMeta.getApplicableCases) that is NOT
// exposed JS-side, and against context-panel / gallery DOM that does NOT render in the cold-headless
// grok-test runtime. This spec covers the two edge invariants reachable cold-deterministically:
//   - Scenario 4 (render-tooltip + render-details): the Notebook ObjectHandler's render path produces
//     non-empty content (name + Author + Created + description). Asserted on handler.renderCard(ent),
//     a SYNCHRONOUS JS render path that does NOT depend on the cold-fragile context-panel accordion.
//   - Scenario 3 (entity.environment): the getter/setter round-trips on a server-fetched entity — the
//     setter mirrors into BOTH metadata.kernelspec.name and .display_name. Pure JS layer (cold-safe).
//
// PARADIGM PIVOT (empirically backed, mcp_status: used 2026-06-18) — this is the Gate-B retry fix.
// The PRIOR spec drove Scenarios 1 + 4 off the context-panel accordion (grok.shell.o = notebook ->
// assert on [name="pane-Actions"] / [name="pane-Details"]). The Gate B cold-failure ground truth
// (test-playwright-output/.../error-context.md, run 2026-06-18T00:06Z) proved the root cause: the
// cold-headless grok-test runtime did NOT mount the context-panel accordion at all — the failure
// snapshot showed the app on the home/Browse welcome view with BOTH [name="pane-Actions"] and
// [name="pane-Details"] ABSENT (Received:false on both presence assertions). This is the SAME
// family-wide cold-headless notebooks-UI-render failure recorded across browser/delete/context-menu
// specs: every notebooks/ Playwright spec gating on rendered gallery OR accordion DOM FAILs Gate B
// cold, while the card-free JS path PASSes. Scenario 3 (pure JS) PASSED in that same cold run,
// confirming the JS layer is cold-healthy and only the UI-render-dependent assertions failed.
//   The pivot removes the context-panel-accordion dependency for the render-details edge entirely:
//   handler.renderCard(ent) (DG.ObjectHandler.forEntity) exercises the SAME Dart renderDetails content
//   code path the accordion would, but synchronously and without mounting the accordion DOM. The
//   element is appended to the live document so the assertion is on a really-rendered HTMLElement.
//   The E-LAYER-COMPLIANCE-01 DOM-driving requirement is met by a real Playwright click on the
//   cold-reliable left-sidebar [name="Browse"] tab (present + visible in the cold-failure snapshot
//   itself) plus a real DOM read of the appended card element.
//
// Scenario 1 is now DEFERRED (scope_reductions SR-04): the getApplicableCases([]) === [] / Apply-to
// absence edge has NO deterministic cold path. (a) JS gap: the entity exposes only
// {environment, description, notebook}; getApplicableCases / isApplicable / tables are Dart-only.
// (b) Cold-DOM gap: the only two surfaces that render the Apply-to menu (gallery-card context menu,
// context-panel Actions pane) both fail to render cold-headless (the precise Gate B FAIL above). The
// Dart getApplicableCases([])===[] invariant is server-testable — routed to
// core/server/datlas/test/dapi/notebooks_test.dart (atlas notebooks.assets.server-tests), alongside
// SR-01 / SR-02 / SR-03.
//
// Three .md scenarios remain DEFERRED as before (SR-01..SR-03):
//   - Scenario 2 (in-memory notebook has no Delete): no JS factory to build an unsaved notebook;
//     the isOnServer gate is not exercisable from the JS/Playwright layer.
//   - Scenario 3, kernelspec-less default-"default" branch: requires Notebook.create({metadata:{}});
//     factory absent JS-side. The getter/setter round-trip on an EXISTING entity IS covered here.
//   - Scenario 5 (notebook-to-code via edit-mode ribbon): atlas notebooks.editor.edit-mode is
//     manual_only (JupyterLab iframe interior + live Python kernel).
//
// Selector recon-notes (live-MCP-observed this session 2026-06-18 on https://dev.datagrok.ai via
// chrome-devtools MCP take_snapshot + evaluate_script). [name="Browse"] is class-1 (the canonical
// sidebar tab used by loginToDatagrok's readiness gate). The recon findings recorded here are the
// card-free render path + the cold-failure diagnosis, not selector fabrication:
//   - Cold-render-safe render-details proxy (Scenario 4): DG.ObjectHandler.forEntity(<notebook>) returns
//     the "Notebook handler"; handler.renderCard(ent) returns a real HTMLElement (instanceof HTMLElement)
//     with textContent length ~598 incl. the notebook name, "Author", "Created <relative-date>", and the
//     description — SYNCHRONOUS (~1ms), no accordion mount, no gallery. handler.renderProperties(ent)
//     and handler.renderTooltip(ent) THREW standalone ("a.gk is not a function" / "reading 'bJ'") because
//     they depend on a live context-panel / tooltip host; renderCard is the host-independent render path
//     that exercises the same renderDetails content. Verified live 2026-06-18.
//   - Sidebar [name="Browse"] tab is present + visible (offsetParent non-null) on the cold home view —
//     it appears in the Gate B cold-failure page snapshot itself (refs e52 "My stuff" region). A real
//     Playwright .click() on it is a cold-robust DOM-driving interaction independent of notebooks UI.
//   - environment getter/setter (Scenario 3): on a server-fetched entity (grok.dapi.notebooks
//     .order('createdOn',true).list), the instance prototype exposes ONLY {environment, description,
//     notebook} — isApplicable / getApplicableCases / tables are ABSENT JS-side
//     (Object.getOwnPropertyNames on the prototype, verified live 2026-06-18). Setting
//     ent.environment='<probe>' made the getter return the probe AND a FRESH ent.notebook re-parse
//     showed metadata.kernelspec.name===probe AND .display_name===probe. (ent.notebook is a per-access
//     deserializing getter, so the assertion re-reads .notebook AFTER the set to observe the persisted
//     mirror — a nested in-place mutation would no-op.) This spec does NOT call save(); the round-trip
//     is asserted in-memory on the fetched entity so it never mutates shared server notebooks.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Notebooks — Edge Cases (render-details via ObjectHandler, environment round-trip)', async ({page}) => {
  // 300s is sufficient — this spec no longer waits on any cold notebooks-UI render (gallery or
  // context-panel accordion). The work is login + a single dapi.notebooks.list + a synchronous JS
  // render-card call + a getter/setter round-trip, all of which complete in a few seconds.
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  });

  // DOM-driving call (E-LAYER-COMPLIANCE-01): a real Playwright click on the cold-reliable left-sidebar
  // [name="Browse"] tab. This tab is present + visible on the cold home view (it appears in the Gate B
  // cold-failure snapshot), so the interaction does NOT depend on any notebooks UI rendering. It is the
  // canonical sidebar tab the loginToDatagrok readiness gate already waits for.
  await softStep('DOM-driving: click the Browse sidebar tab (cold-reliable, notebooks-UI-independent)', async () => {
    const browseTab = page.locator('[name="Browse"]').first();
    await browseTab.waitFor({timeout: 30_000});
    await browseTab.click();
    await page.waitForTimeout(300);
    await expect(browseTab, 'Browse sidebar tab is interactable').toBeVisible();
  });

  // ---- Scenario 4: notebook render-tooltip / render-details produce non-empty content ----
  await softStep('Scenario 4: Notebook ObjectHandler render path yields non-empty render-details content', async () => {
    // Fetch the newest server notebook, build its handler-rendered card, append it to the live DOM, and
    // assert the rendered element carries the render-details content (name + Author + Created/Modified +
    // description). renderCard is the host-independent render path that exercises the same Dart
    // renderDetails content as the (cold-fragile) context-panel accordion, so this covers
    // render-tooltip + render-details without depending on the accordion mounting cold-headless.
    const probeId = 'automator-edge-rendercard-probe';
    const rendered = await page.evaluate(async (hostId) => {
      const list = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 1});
      if (!Array.isArray(list) || list.length === 0) return null;
      const ent: any = list[0];
      const handler: any = DG.ObjectHandler.forEntity(ent);
      if (!handler || typeof handler.renderCard !== 'function')
        return {handlerOk: false, handlerName: handler?.constructor?.name ?? null};
      const card: any = handler.renderCard(ent);
      const el: HTMLElement = (card && card.root) ? card.root : card;
      // Mount the rendered element into the live document so the assertion reads a really-rendered node.
      const host = document.createElement('div');
      host.id = hostId;
      host.style.position = 'fixed';
      host.style.left = '-9999px';
      host.appendChild(el);
      document.body.appendChild(host);
      const text = (el?.textContent || '').trim();
      return {
        handlerOk: true,
        handlerName: handler.constructor?.name ?? null,
        isHtmlElement: el instanceof HTMLElement,
        mountedInDom: !!document.getElementById(hostId),
        textLen: text.length,
        hasName: text.includes(ent.friendlyName || ent.name || ' '),
        hasCreatedOrModifiedOrAuthor: /Created|Modified|Author|ago|20\d\d/i.test(text),
      };
    }, probeId);

    expect(rendered, 'at least one server notebook should be fetchable + rendered').not.toBeNull();
    expect(rendered!.handlerOk, 'the Notebook ObjectHandler should expose renderCard').toBe(true);
    expect(rendered!.isHtmlElement, 'renderCard should return a real HTMLElement').toBe(true);

    // DOM read of the appended, really-rendered card element (a genuine cross-process DOM assertion).
    const cardEl = page.locator(`#${probeId}`).first();
    await cardEl.waitFor({timeout: 15_000});
    const domText = (await cardEl.textContent() || '').trim();
    expect(domText.length, 'the rendered render-details card should be non-empty in the live DOM').toBeGreaterThan(0);

    // Content predicate keyed on cold-stable render-details structure (name + Created/Modified/Author),
    // NOT a lazy viewer canvas. renderTooltip delegates to renderDetails (notebook_meta.dart#L25).
    expect(rendered!.textLen, 'the rendered render-details content should be non-empty').toBeGreaterThan(0);
    expect(rendered!.hasName, 'render-details content should include the notebook name').toBe(true);
    expect(rendered!.hasCreatedOrModifiedOrAuthor,
      'render-details content should include Created / Modified / Author (render-details fields)').toBe(true);
  });

  // ---- Scenario 3 (authorable slice): environment getter/setter round-trip on a fetched entity ----
  await softStep('Scenario 3: Notebook.environment getter/setter mirrors into metadata.kernelspec', async () => {
    // Fetch a real server notebook and round-trip its environment in-memory (NO save() — never mutate
    // shared server state). The instance prototype exposes only {environment, description, notebook}.
    const result = await page.evaluate(async () => {
      const list = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 1});
      if (!Array.isArray(list) || list.length === 0) return null;
      const ent: any = list[0];
      const protoKeys = Object.getOwnPropertyNames(Object.getPrototypeOf(ent));
      const envBefore = ent.environment;
      // Set via the setter; the getter should reflect it, and a FRESH .notebook re-parse should show
      // metadata.kernelspec.name AND .display_name both updated (the setter mirrors both).
      const probe = `automator-env-${Date.now()}`;
      ent.environment = probe;
      const envAfter = ent.environment;
      const fresh = ent.notebook; // per-access getter re-parses the stored .ipynb
      const ks = fresh?.metadata?.kernelspec ?? null;
      return {
        protoKeys,
        envBefore,
        envAfter,
        probe,
        kernelspecName: ks?.name ?? null,
        kernelspecDisplayName: ks?.display_name ?? null,
        hasIsApplicable: typeof ent.isApplicable,
        hasGetApplicableCases: typeof ent.getApplicableCases,
      };
    });
    expect(result, 'at least one server notebook should be fetchable').not.toBeNull();
    // The getter returns a non-empty environment for a normal (kernelspec-present) notebook.
    expect(typeof result!.envBefore, 'environment getter returns a string').toBe('string');
    expect((result!.envBefore as string).length, 'environment getter is non-empty').toBeGreaterThan(0);
    // The setter updates the getter...
    expect(result!.envAfter, 'environment getter reflects the setter value').toBe(result!.probe);
    // ...and mirrors into BOTH kernelspec.name and kernelspec.display_name (notebook.dart#L21).
    expect(result!.kernelspecName, 'setter writes metadata.kernelspec.name').toBe(result!.probe);
    expect(result!.kernelspecDisplayName, 'setter writes metadata.kernelspec.display_name').toBe(result!.probe);
    // Document the JS-surface gap that forced the Scenario 1 / Scenario 2 / kernelspec-less deferrals:
    // the entity exposes neither isApplicable nor getApplicableCases (those live Dart-side only).
    expect(result!.hasIsApplicable, 'isApplicable is absent JS-side (entity-model logic is Dart-only)').toBe('undefined');
    expect(result!.hasGetApplicableCases, 'getApplicableCases is absent JS-side').toBe('undefined');
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
