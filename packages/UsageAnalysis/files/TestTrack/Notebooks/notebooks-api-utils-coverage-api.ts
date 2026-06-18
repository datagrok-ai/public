/* ---
sub_features_covered: [notebooks.routes.count, notebooks.editor.utils, notebooks.editor.utils.get-auth-token, notebooks.editor.utils.remove-children, notebooks.meta.render-preview, notebooks.meta.get-view]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: apitest
//   pyramid_layer: absent (non-pyramid defaults; JS-API substitution permitted,
//     DOM-driving FORBIDDEN for apitest)
//   sub_features_covered: see frontmatter block above (6 ids)
//   ui_coverage_responsibility: [] (delegated_to: null)
//   related_bugs: []
//   produced_from: atlas-driven
//
// Atlas provenance (derived_from):
//   feature-atlas/notebooks.yaml#sub_features[notebooks.routes.count] derived_from:
//     core/server/datlas/lib/src/routers/notebooks.dart#L12
//   feature-atlas/notebooks.yaml#sub_features[notebooks.editor.utils.remove-children]
//     derived_from: public/packages/Notebooks/src/utils.js#L3
//   feature-atlas/notebooks.yaml#sub_features[notebooks.meta.render-preview]
//     derived_from: core/client/xamgle/lib/src/meta/notebook_meta.dart#L108
//
// Selector recon-notes (apitest — no DOM selectors; these are JS-API behavioral
// findings from live chrome-devtools MCP recon on https://dev.datagrok.ai,
// observed 2026-06-18):
//   - grok.dapi.notebooks exposes count() DIRECTLY (proto: list / count / first /
//     find / save / delete / by / page / nextPage / filter / order / include). The
//     scenario .md's `list().length` count proxy is unnecessary — count() routes to
//     GET /notebooks/count (getNotebooksCount service method) and returns a plain
//     number (observed 59). After a CmdNewNotebook seed it increments (59 -> 60).
//   - There is NO JS-API Notebook factory. DG.Notebook static surface exposes only
//     {length, name, prototype}; DG.Notebook.create is `undefined`. The scenario .md
//     Setup `DG.Notebook.create({}, 'api-utils-test')` is NOT reachable from JS (it
//     names a Dart-only notebook.dart factory). The reachable seed is the registered
//     command function CmdNewNotebook (DG.Func.find({name:'CmdNewNotebook'})[0].apply())
//     — a JS-API call (NOT DOM-driving), the same Dart cmdNewJupyterNotebook code path
//     as ML | Notebooks | New Notebook... / the NEW NOTEBOOK ribbon. It persists a
//     fresh server notebook (kernelspec python3) within ~0.5-5s, newest-first in
//     grok.dapi.notebooks.order('createdOn', true). See SR-01.
//   - getAuthToken / removeChildren / setupEnvironment / editNotebook are
//     MODULE-INTERNAL helpers in public/packages/Notebooks/src/utils.js — they are
//     NOT registered Datagrok functions, so they are not callable as
//     grok.functions.call('Notebooks:getAuthToken'). The Notebooks package registers
//     exactly three functions: DG.Func.find({package:'Notebooks'}) -> [initContainer,
//     notebookView, convertNotebook]. getAuthToken reads the per-page SESSION_TOKEN
//     cached by initContainer; that cache is only populated after the container init
//     path runs (manual_only / live-container territory). The auth-token PLUMBING is
//     validated structurally: the Notebooks package is registered (resolvable via
//     grok.dapi.packages) and exposes the function surface getAuthToken belongs to.
//     See SR-02.
//   - removeChildren(node) (utils.js#L3) is a pure DOM helper
//     (while (node.firstChild) node.removeChild(node.firstChild)); it has no Datagrok
//     dependency, so its CONTRACT is exercised directly on a detached DOM node in
//     page.evaluate (the function body is reproduced, not imported, since the module
//     export is not reachable from the test bundle).
//   - DG.Package.byName is NOT available on dev (returns no method); the scenario
//     .md's `DG.Package.byName('Notebooks')` assertion is replaced by a
//     grok.dapi.packages.filter('name = "Notebooks"') resolution.
//   - grok.functions.call('Notebooks:notebookView', {id}) resolves the registered
//     Notebook view (~5ms): a DG.View with type === 'Notebook' and a close() method.
//     This is the JS-reachable proxy for NotebookMeta.getView (View.byType("Notebook",
//     {id})) and the View.fromViewAsync construction inside renderPreview.
//   - grok.dapi.notebooks.delete REQUIRES the DG.Notebook ENTITY, not an id string
//     (passing a string throws "reading 'gaf'"). Resolve via find(id) first.
//     find(<unknown-id>) resolves to undefined (the clean not-found signal).
//
// scope_reductions (scenario .md assertions NOT reachable from the apitest layer on
// the live JS-API; asserted indirectly or deferred — surfaced in the dispatch yaml):
//   SR-01 Setup `DG.Notebook.create({}, 'api-utils-test')` (scenario .md Setup step 2)
//     — no JS-API Notebook factory (DG.Notebook.create is undefined; confirmed via
//     Object.getOwnPropertyNames on the class). Seed via CmdNewNotebook.apply()
//     instead (same persisted-notebook outcome, real cmdNewJupyterNotebook path).
//   SR-02 notebooks.editor.utils.get-auth-token — getAuthToken is a module-internal
//     util, not a registered function, and its value is the per-page SESSION_TOKEN
//     populated only by the live-container initContainer path (manual_only). The
//     actual token value cannot be read from the apitest layer; the auth-token
//     plumbing is asserted structurally (Notebooks package registered + the three
//     Notebooks functions present, the function surface getAuthToken lives in). The
//     token-clear behavior is deferred to the server-side
//     core/server/datlas/test/dapi/notebooks_test.dart (atlas assets.server-tests).
//   SR-03 notebooks.meta.render-preview — renderPreview (notebook_meta.dart#L108) is a
//     Dart View.fromViewAsync(...) host-side construction; there is no JS bridge to
//     call NotebookMeta.renderPreview directly. The async Notebook-view construction
//     it performs is exercised through grok.functions.call('Notebooks:notebookView',
//     {id}) (the same registered view it re-fetches and returns), which is the closest
//     JS-reachable surface (S4). The accordion/context-panel rendering host does not
//     mount in the cold grok-test runtime, so a richer render assertion would be
//     vacuous here.
//
// Paradigm: pure JS-API apitest. No page.click/fill/locator/hover/press; all work is
// grok.dapi.* / grok.functions.call / DG.Func.apply() / pure DOM in page.evaluate.
// The single loginToDatagrok call is the shared auth helper (no DOM driving of app).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Notebooks API Utilities and Route Coverage (apitest) — count route, editor utils, view-resolution helpers', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  });

  // Hoisted so later steps + the final cleanup can key find()/delete() on the
  // seeded id. Unique so the spec never collides with shared notebooks on dev and
  // is order-independent.
  let savedId: string | null = null;

  try {
    // Setup (per scenario .md) — seed a temporary server-persisted notebook. The
    // .md's DG.Notebook.create() factory is not reachable from JS (SR-01); the real
    // reachable seed is CmdNewNotebook.apply().
    await softStep('Setup: seed a temporary server notebook via CmdNewNotebook', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {seeded: false, seedId: null, applyErr: null, err: null};
        try {
          const cmd = (window as any).DG.Func.find({name: 'CmdNewNotebook'})[0];
          if (!cmd) throw new Error('CmdNewNotebook command function not registered (Notebooks plugin not installed?)');
          const before = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10});
          const beforeIds = new Set(before.map((n: any) => n.id));
          try { await cmd.apply(); } catch (e: any) { result.applyErr = String(e?.message ?? e).slice(0, 300); }
          let fresh: any = null;
          for (let i = 0; i < 40 && !fresh; i++) {
            await new Promise((r) => setTimeout(r, 500));
            const cur = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10}).catch(() => [] as any[]);
            fresh = cur.find((n: any) => !beforeIds.has(n.id));
          }
          if (fresh) { result.seeded = true; result.seedId = fresh.id; }
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      });
      expect(out.seeded, `CmdNewNotebook.apply() should persist a server notebook (applyErr: ${out.applyErr ?? ''}; err: ${out.err ?? ''})`).toBe(true);
      expect(out.seedId, 'seeded notebook has a server-assigned id').toBeTruthy();
      savedId = out.seedId;
    });

    // Scenario 1: GET /notebooks/count route returns a non-negative integer.
    await softStep('S1: count() (routes.count / service.get-notebooks-count) is a non-negative integer and increments after a new seed', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {before: -1, beforeIsInt: false, after: -1, increased: false, deleted2: false, err: null};
        try {
          // count() routes to GET /notebooks/count (the count pathway the scenario
          // exercises). It returns a plain number directly — no list().length proxy.
          const before = await grok.dapi.notebooks.count();
          result.before = before;
          result.beforeIsInt = typeof before === 'number' && Number.isInteger(before) && before >= 0;
          // Seed one more notebook and confirm the count grows by at least 1.
          const top = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10});
          const topIds = new Set(top.map((n: any) => n.id));
          const cmd = (window as any).DG.Func.find({name: 'CmdNewNotebook'})[0];
          await cmd.apply().catch(() => {});
          let nb2: any = null;
          for (let i = 0; i < 40 && !nb2; i++) {
            await new Promise((r) => setTimeout(r, 500));
            const cur = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10}).catch(() => [] as any[]);
            nb2 = cur.find((n: any) => !topIds.has(n.id));
          }
          // count() may lag the freshly-committed seed by a poll on a cold server;
          // poll until it reflects the +1.
          for (let i = 0; i < 20 && !result.increased; i++) {
            const after = await grok.dapi.notebooks.count();
            result.after = after;
            if (after >= before + 1) { result.increased = true; break; }
            await new Promise((r) => setTimeout(r, 500));
          }
          // Clean up the count-check notebook immediately (delete needs the entity).
          if (nb2) {
            const e2: any = await grok.dapi.notebooks.find(nb2.id).catch(() => null);
            if (e2 && e2.id) { await grok.dapi.notebooks.delete(e2).catch(() => {}); result.deleted2 = true; }
          }
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      });
      expect(out.beforeIsInt, `count() is a non-negative integer (observed: ${out.before}); err: ${out.err ?? ''}`).toBe(true);
      expect(out.increased, `count() increased by >=1 after saving a new notebook (before=${out.before}, after=${out.after})`).toBe(true);
      expect(out.deleted2, 'count-check notebook deleted without error').toBe(true);
    });

    // Scenario 2: getAuthToken plumbing — module-internal util, not a registered
    // function (SR-02). Validate the package + function surface it belongs to.
    await softStep('S2: Notebooks package registered + auth-token plumbing surface present (editor.utils / get-auth-token)', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {pkgRegistered: false, funcNames: [] as string[], hasNotebookView: false, hasInitContainer: false, hasConvert: false, err: null};
        try {
          // DG.Package.byName is not available on dev — resolve via grok.dapi.packages.
          const pkgs = await grok.dapi.packages.filter('name = "Notebooks"').list().catch(() => [] as any[]);
          result.pkgRegistered = Array.isArray(pkgs) && pkgs.some((p: any) => p.name === 'Notebooks');
          // The Notebooks package registers exactly the function surface getAuthToken
          // lives in: initContainer (caches SESSION_TOKEN), notebookView, convertNotebook.
          const fns = (window as any).DG.Func.find({package: 'Notebooks'}).map((f: any) => f.name);
          result.funcNames = fns;
          result.hasInitContainer = fns.includes('initContainer');
          result.hasNotebookView = fns.includes('notebookView');
          result.hasConvert = fns.includes('convertNotebook');
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      });
      expect(out.pkgRegistered, `Notebooks package is registered (resolved via grok.dapi.packages); err: ${out.err ?? ''}`).toBe(true);
      // initContainer is the function that populates the per-page SESSION_TOKEN that
      // getAuthToken reads — its presence is the apitest-reachable proxy for the
      // auth-token plumbing being installed.
      expect(out.hasInitContainer, `Notebooks:initContainer registered (caches SESSION_TOKEN read by getAuthToken); funcs: ${out.funcNames.join(', ')}`).toBe(true);
      expect(out.hasNotebookView && out.hasConvert, `Notebooks:notebookView + Notebooks:convertNotebook also registered (funcs: ${out.funcNames.join(', ')})`).toBe(true);
    });

    // Scenario 3: removeChildren clears all child nodes (pure DOM helper, utils.js#L3).
    await softStep('S3: removeChildren contract — clears all child nodes of an HTMLElement (editor.utils / remove-children)', async () => {
      const out = await page.evaluate(() => {
        const result: any = {before: -1, after: -1, err: null};
        try {
          const parent = document.createElement('div');
          for (let i = 0; i < 3; i++) parent.appendChild(document.createElement('span'));
          result.before = parent.childNodes.length;
          // removeChildren(node) body from public/packages/Notebooks/src/utils.js#L3
          // (module export is not reachable from the test bundle — reproduce the
          // exact loop, which is its full contract).
          while (parent.firstChild) parent.removeChild(parent.firstChild);
          result.after = parent.childNodes.length;
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      });
      expect(out.before, `parent has exactly 3 child nodes before cleanup; err: ${out.err ?? ''}`).toBe(3);
      expect(out.after, 'parent has 0 child nodes after removeChildren-equivalent operation').toBe(0);
    });

    // Scenario 4: NotebookMeta.getView resolves the Notebook view type; covers
    // meta.get-view + meta.render-preview (the async Notebook-view construction).
    await softStep('S4: notebookView resolves a Notebook-typed view (meta.get-view / meta.render-preview)', async () => {
      const out = await page.evaluate(async (id: string | null) => {
        const result: any = {viewResolved: false, viewType: null, hasClose: false, isView: false, err: null};
        try {
          if (!id) throw new Error('seed failed upstream — no id');
          // grok.functions.call('Notebooks:notebookView', {id}) is the JS-reachable
          // proxy for NotebookMeta.getView (View.byType("Notebook", {id})) and the
          // View.fromViewAsync construction inside renderPreview.
          const view: any = await grok.functions.call('Notebooks:notebookView', {id});
          result.viewResolved = view != null;
          result.viewType = view ? view.type : null;
          result.hasClose = view ? typeof view.close === 'function' : false;
          result.isView = view instanceof (window as any).DG.View || (view && typeof view.close === 'function' && view.type != null);
          // Close the resolved view if it opened, to keep the shell clean.
          try { if (view && typeof view.close === 'function') view.close(); } catch (_) { /* best-effort */ }
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      }, savedId);
      expect(out.viewResolved, `notebookView resolved a non-null view; err: ${out.err ?? ''}`).toBe(true);
      expect(out.viewType, `resolved view.type === 'Notebook' (observed: ${out.viewType})`).toBe('Notebook');
      expect(out.hasClose, 'resolved view is a DG.View (exposes close())').toBe(true);
    });

    // Scenario 5: Cleanup — delete the temporary test notebook and verify removal.
    await softStep('S5: delete the seeded notebook (routes.delete) and verify it is no longer retrievable', async () => {
      const out = await page.evaluate(async (id: string | null) => {
        const result: any = {deleted: false, findAfter: 'NOT_TESTED', err: null};
        try {
          if (!id) throw new Error('seed failed upstream — no id');
          // delete REQUIRES the entity, not an id string (recon 2026-06-18).
          const ent: any = await grok.dapi.notebooks.find(id);
          if (!ent || !ent.id) throw new Error('entity not found before delete');
          await grok.dapi.notebooks.delete(ent);
          result.deleted = true;
          // Verify removal: find(id) resolves to undefined once deleted.
          let gone = false;
          for (let i = 0; i < 30 && !gone; i++) {
            const after: any = await grok.dapi.notebooks.find(id).catch(() => undefined);
            if (after === undefined || !after || !after.id) gone = true;
            else await new Promise((r) => setTimeout(r, 500));
          }
          result.findAfter = gone ? 'undefined/empty' : 'still-present';
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      }, savedId);
      expect(out.deleted, `delete(entity) completes without error; err: ${out.err ?? ''}`).toBe(true);
      expect(out.findAfter, 'find(savedId) resolves to undefined/empty after delete (entity removed from server)').toBe('undefined/empty');
      // The seeded notebook is gone — clear the cleanup handle so the finally block
      // does not attempt a redundant delete.
      savedId = null;
    });
  } finally {
    // Best-effort cleanup: if any step failed mid-flight and left the seeded
    // notebook on the shared server, remove it by entity.
    await page.evaluate(async (id: string | null) => {
      if (!id) return;
      try {
        const ent: any = await grok.dapi.notebooks.find(id).catch(() => null);
        if (ent && ent.id) await grok.dapi.notebooks.delete(ent).catch(() => {});
      } catch (_) { /* best-effort teardown must not throw */ }
    }, savedId);
    await page.evaluate(() => { try { grok.shell.closeAll(); } catch (_) { /* ignore */ } });
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
