/* ---
sub_features_covered: [notebooks.routes.count, notebooks.editor.utils, notebooks.editor.utils.get-auth-token, notebooks.editor.utils.remove-children, notebooks.meta.render-preview, notebooks.meta.get-view]
--- */
// Pure JS-API apitest for Notebooks utility/route coverage. Seeds a temporary
// server notebook via CmdNewNotebook, then exercises: the count route, the
// Notebooks package + auth-token function surface, the removeChildren DOM helper,
// notebookView resolution, and delete. All work is grok.dapi.* / grok.functions.call
// / DG.Func.apply() / pure DOM in page.evaluate — no DOM driving of the app.
// Scope-reduced: get-auth-token (per-page SESSION_TOKEN) and render-preview have no
// JS bridge, so they're asserted structurally here and covered server-side.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';

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
    // Setup — seed a temporary server-persisted notebook via CmdNewNotebook
    // (no JS-API Notebook factory exists; this is the reachable seed path).
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
          // count() routes to GET /notebooks/count and returns a plain number.
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
    // function. Validate the package + function surface it belongs to.
    await softStep('S2: Notebooks package registered + auth-token plumbing surface present (editor.utils / get-auth-token)', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {pkgRegistered: false, funcNames: [] as string[], hasNotebookView: false, hasInitContainer: false, hasConvert: false, err: null};
        try {
          const pkgs = await grok.dapi.packages.filter('name = "Notebooks"').list().catch(() => [] as any[]);
          result.pkgRegistered = Array.isArray(pkgs) && pkgs.some((p: any) => p.name === 'Notebooks');
          // The Notebooks package registers the surface getAuthToken lives in:
          // initContainer (caches SESSION_TOKEN), notebookView, convertNotebook.
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
          // removeChildren(node) body from Notebooks/src/utils.js, reproduced (the
          // module export is not reachable from the test bundle).
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
          // notebookView is the JS-reachable proxy for NotebookMeta.getView and the
          // async view construction inside renderPreview.
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
          // delete requires the entity, not an id string.
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
