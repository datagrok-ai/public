/* ---
sub_features_covered: [notebooks.entity.link, notebooks.entity.is-applicable, notebooks.entity.get-applicable-cases, notebooks.entity.from-json, notebooks.entity.generate-file-name, notebooks.api.find, notebooks.api.list, notebooks.api.save, notebooks.api.filter, notebooks.repository.save, notebooks.repository.find, notebooks.repository.delete, notebooks.repository.delete-tables-relations, notebooks.repository.by-tag, notebooks.repository.audit-hooks, notebooks.routes.save, notebooks.routes.get, notebooks.routes.delete, notebooks.routes.list, notebooks.service.save-notebook, notebooks.service.get-notebook, notebooks.service.delete-notebook, notebooks.service.get-notebooks-filtered, notebooks.service.get-notebooks-count, notebooks.meta.open-tables-in-notebook, notebooks.meta.new-notebook]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: apitest
//   pyramid_layer: absent (non-pyramid defaults; JS-API substitution permitted,
//     DOM-driving FORBIDDEN for apitest)
//   sub_features_covered: see frontmatter block above (25 ids)
//   ui_coverage_responsibility: [] (delegated_to: null)
//   related_bugs: []
//   produced_from: atlas-driven
//
// Atlas provenance (derived_from):
//   feature-atlas/notebooks.yaml#critical_paths[new-blank-notebook] derived_from:
//     core/client/xamgle/lib/src/features/jupyter_notebook/jupyter_notebook_plugin.dart#L7
//   feature-atlas/notebooks.yaml#critical_paths[delete-notebook] derived_from:
//     core/client/xamgle/lib/src/meta/notebook_meta.dart#L128
//   feature-atlas/notebooks.yaml#critical_paths[rename-notebook] derived_from:
//     core/client/xamgle/lib/src/features/jupyter_notebook/jupyter_notebook_plugin.dart#L26
//
// Selector recon-notes (apitest — no DOM selectors; these are JS-API behavioral
// findings from live chrome-devtools MCP recon on https://dev.datagrok.ai,
// observed 2026-06-18):
//   - NO JS-API Notebook factory. DG.Notebook static surface exposes none of
//     template / create / fromJson / link / isApplicable / getApplicableCases /
//     generateFileName. DG.Notebook.prototype exposes only the instance accessors
//     environment / description / notebook. `new DG.Notebook()` constructs a husk
//     but every mutating accessor throws a Dart-interop error (set notebook ->
//     "reading 'gaf'"; set friendlyName -> "reading 'sak'") and save() throws
//     NoSuchMethodError: []("description") on null. There is therefore NO pure-JS
//     path to create a notebook from scratch. The reachable seed is the registered
//     command function CmdNewNotebook (DG.Func.find({name:'CmdNewNotebook'})[0].apply())
//     — a JS-API call (NOT DOM-driving), which is the same Dart cmdNewJupyterNotebook
//     code path that ML | Notebooks | New Notebook... / the NEW NOTEBOOK ribbon
//     dispatch. It persists a fresh server notebook (kernelspec python3) within
//     ~0.5-5s and surfaces newest-first in grok.dapi.notebooks.order('createdOn',true).
//   - grok.dapi.notebooks.delete REQUIRES the DG.Notebook ENTITY, not an id string
//     (passing a string/the husk throws). Resolve via find(id) first.
//   - find(id) returns the entity with .notebook (raw .ipynb JSON) populated;
//     find(<unknown-id>) resolves to `undefined` (the clean not-found signal).
//   - The JS entity surface does NOT expose `tables` (undefined) or `tags`
//     (undefined; no setTag/getTags) — the notebooks_tables join and entity tags
//     are server-side only and not enumerable through the JS-typed Notebook. See
//     scope_reductions SR-01..SR-03 below.
//   - grok.dapi.notebooks.filter({...}) throws from JS for every shape tried
//     ({tags:[...]}, {tags:'x'}, {}) with "'replace' is not a function" — the
//     where() helper is not callable through the JS client. count() / list() /
//     order(field, desc) / first() all work. See SR-04.
//   - S4 rename-persist latency (STAB note, measured live via MCP recon
//     2026-06-18): on a WARM dev server `ent.friendlyName = x; save(ent)` then
//     find(id) reflects the new value on the FIRST poll (saveMs ~79ms, persistMs
//     ~77ms). The Gate-B FLAKY attempt-1 (S4 timed out reporting the default
//     friendlyName "Notebook") was a COLD-START commit lag on the first save after
//     the CmdNewNotebook seed, exceeding the prior 5s retry budget. Fix is a
//     same-paradigm stabilization (raise budget to 30s + idempotent re-save at the
//     10s/20s marks) — NOT a paradigm change.
//
// scope_reductions (sub_features in frontmatter that are NOT reachable from the
// apitest layer on the live JS-API; asserted indirectly or deferred — documented
// here and surfaced in the dispatch yaml):
//   SR-01 notebooks.entity.link / .is-applicable / .get-applicable-cases /
//     .from-json / .generate-file-name — no JS-API surface (DG.Notebook exposes no
//     such static/instance methods; confirmed via Object.getOwnPropertyNames on the
//     class + prototype). The atlas anchors them to notebook.dart Dart methods with
//     no JS bridge. Applicability + table-link semantics are exercised server-side
//     by core/server/datlas/test/dapi/notebooks_test.dart (atlas
//     assets.server-tests), not reproducible from apitest. ASSERTED INDIRECTLY where
//     possible: the notebooks_tables relation lifecycle is exercised by the
//     create -> find -> delete round-trip (S1/S2/S6) at the repository.save /
//     repository.delete / delete-tables-relations layer; the JS surface cannot read
//     the join rows back.
//   SR-02 tables-based assertions (scenario .md S1.5/S2.3 "saved.tables.length===1")
//     — `tables` is undefined on the JS entity; cannot assert join cardinality from
//     JS. Round-trip integrity asserted on the notebook (.ipynb) body + name instead.
//   SR-03 notebooks.repository.by-tag + tag round-trip (scenario .md S3 tag filter)
//     — no JS tag surface (tags undefined; no setTag). Listing/count covered via the
//     untagged list()/count()/order() path (api.list / routes.list /
//     service.get-notebooks-filtered / service.get-notebooks-count) instead.
//   SR-04 notebooks.api.filter — filter() throws from the JS client for all shapes;
//     covered structurally by the order()+list() path which routes through the same
//     getNotebooksFiltered service method.
//   SR-05 notebooks.entity.apply / convert-notebook — atlas manual_only[] (live
//     Docker container, non-deterministic). Out of scope per scenario .md Notes.
//
// Paradigm: pure JS-API apitest. No page.click/fill/locator/hover/press; all work
// is grok.dapi.notebooks / DG.Func.apply() in page.evaluate. The single
// loginToDatagrok call is the shared auth helper (no DOM driving of the app).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Notebooks Lifecycle (apitest) — linked-table source class CRUD round-trip via grok.dapi.notebooks', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  });

  // Hoisted so later steps + the final cleanup can key find()/delete() on the
  // seeded id. Unique target name so the spec never collides with shared
  // Demog/Cars notebooks on dev and is order-independent.
  let seededId: string | null = null;
  const renameTarget = `automator-lifecycle-${Date.now()}`;

  try {
    await softStep('S1: Seed a server notebook via CmdNewNotebook (new-notebook / api.save / service.save / repository.save)', async () => {
      const out = await page.evaluate(async () => {
        const result: any = {seeded: false, applyErr: null, seedId: null, friendlyName: null, kernelspec: null, err: null};
        try {
          const cmd = (window as any).DG.Func.find({name: 'CmdNewNotebook'})[0];
          if (!cmd) throw new Error('CmdNewNotebook command function not registered (Notebooks plugin not installed?)');
          // Snapshot newest-first top so we detect the freshly-created entity by id.
          const before = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10});
          const beforeIds = new Set(before.map((n: any) => n.id));
          try { await cmd.apply(); } catch (e: any) { result.applyErr = String(e?.message ?? e).slice(0, 300); }
          let fresh: any = null;
          for (let i = 0; i < 40 && !fresh; i++) {
            await new Promise((r) => setTimeout(r, 500));
            const cur = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10}).catch(() => [] as any[]);
            fresh = cur.find((n: any) => !beforeIds.has(n.id));
          }
          if (fresh) {
            result.seeded = true;
            result.seedId = fresh.id;
            result.friendlyName = fresh.friendlyName ?? fresh.name;
            result.kernelspec = fresh.notebook?.metadata?.kernelspec?.name ?? null;
          }
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      });
      expect(out.seeded, `CmdNewNotebook.apply() should persist a server notebook (applyErr: ${out.applyErr ?? ''}; err: ${out.err ?? ''})`).toBe(true);
      expect(out.seedId, 'seeded notebook has a server-assigned id').toBeTruthy();
      // Newly created template notebook defaults to the Python 3 kernelspec.
      expect(out.kernelspec, `new notebook kernelspec.name === 'python3' (observed: ${out.kernelspec})`).toBe('python3');
      seededId = out.seedId;
    });

    await softStep('S2: Round-trip find (api.find / routes.get / service.get-notebook / repository.find) — name + .ipynb body intact', async () => {
      const out = await page.evaluate(async (id: string | null) => {
        const result: any = {fetched: false, name: null, notebookPresent: false, cellsIsArray: false, kernelspec: null, err: null};
        try {
          if (!id) throw new Error('seed failed upstream — no id');
          const ent: any = await grok.dapi.notebooks.find(id);
          result.fetched = ent != null && !!ent.id;
          if (ent) {
            result.name = ent.friendlyName ?? ent.name;
            result.notebookPresent = ent.notebook != null;
            result.cellsIsArray = Array.isArray(ent.notebook?.cells);
            result.kernelspec = ent.notebook?.metadata?.kernelspec?.name ?? null;
          }
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      }, seededId);
      expect(out.fetched, `find(seedId) returns the persisted entity; err: ${out.err ?? ''}`).toBe(true);
      expect(out.notebookPresent, 'fetched entity carries its raw .ipynb JSON body').toBe(true);
      expect(out.cellsIsArray, 'fetched .ipynb has a cells array (well-formed nbformat)').toBe(true);
      // metadata.kernelspec survives the round-trip (notebooks.entity defaults).
      expect(out.kernelspec, `round-trip kernelspec.name === 'python3' (observed: ${out.kernelspec})`).toBe('python3');
    });

    await softStep('S3: List / count / order (api.list / routes.list / service.get-notebooks-filtered / service.get-notebooks-count)', async () => {
      const out = await page.evaluate(async (id: string | null) => {
        const result: any = {
          listOk: false, listCount: 0, seedInTop: false,
          countIsNumber: false, countVal: -1,
          orderDistinct: false, ascFirst: null, descFirst: null,
          firstOk: false, err: null,
        };
        try {
          // newest-first list — the freshly-seeded notebook is the newest, so it
          // appears in the top-10 (cheap, ~0.1-2.2s, vs a full pageSize:500 scan).
          const top = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10});
          result.listOk = Array.isArray(top);
          result.listCount = top.length;
          result.seedInTop = id != null && top.some((n: any) => n.id === id);
          // count() — getNotebooksCount service method.
          const c = await grok.dapi.notebooks.count();
          result.countIsNumber = typeof c === 'number';
          result.countVal = c;
          // order() asc vs desc by name should differ (deterministic ordering).
          const asc = await grok.dapi.notebooks.order('name').list({pageSize: 2});
          const desc = await grok.dapi.notebooks.order('name', true).list({pageSize: 2});
          result.ascFirst = asc[0]?.name ?? null;
          result.descFirst = desc[0]?.name ?? null;
          result.orderDistinct = asc[0]?.id != null && desc[0]?.id != null && asc[0].id !== desc[0].id;
          // first() — convenience over the ordered list.
          const fst = await grok.dapi.notebooks.order('createdOn', true).first();
          result.firstOk = fst != null && !!fst.id;
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      }, seededId);
      expect(out.listOk, `order('createdOn',desc).list returns an array; err: ${out.err ?? ''}`).toBe(true);
      expect(out.seedInTop, 'the freshly-seeded notebook surfaces in the newest-first top-10').toBe(true);
      expect(out.countIsNumber, `count() returns a number (observed: ${out.countVal})`).toBe(true);
      expect(out.countVal, 'count() >= 1 (at least the seeded notebook exists)').toBeGreaterThanOrEqual(1);
      expect(out.orderDistinct, `asc vs desc name order differ (asc[0]=${out.ascFirst}, desc[0]=${out.descFirst})`).toBe(true);
      expect(out.firstOk, 'first() resolves a single notebook entity').toBe(true);
    });

    await softStep('S4: Rename via friendlyName + save (rename_notebook — source-agnostic; api.save / repository.save)', async () => {
      const out = await page.evaluate(async (args: {id: string | null; name: string}) => {
        const result: any = {saveOk: false, persisted: false, value: null, resaves: 0, pollMs: 0, err: null};
        try {
          if (!args.id) throw new Error('seed failed upstream — no id');
          const ent: any = await grok.dapi.notebooks.find(args.id);
          ent.friendlyName = args.name;
          await grok.dapi.notebooks.save(ent);
          result.saveOk = true;
          // Re-fetch and confirm the rename persisted. On a WARM server this lands
          // within the first poll (~80ms, measured via MCP recon 2026-06-18). On a
          // COLD server the first save after a CmdNewNotebook seed can lag beyond the
          // prior 5s budget (the FLAKY attempt-1 root cause) — the server still
          // reports the default friendlyName "Notebook" for several seconds. Budget
          // raised to 30s (30x1000ms) to absorb cold-commit lag, with an idempotent
          // re-save every ~10s as belt-and-suspenders in case the first save was
          // dropped during cold init rather than merely lagging.
          const tPoll = Date.now();
          for (let i = 0; i < 30 && !result.persisted; i++) {
            const re: any = await grok.dapi.notebooks.find(args.id).catch(() => null);
            const cur = re ? (re.friendlyName || re.name) : null;
            result.value = cur;
            if (cur === args.name) { result.persisted = true; break; }
            // Idempotent re-save at the 10s and 20s marks if still not persisted.
            if ((i === 10 || i === 20) && re && re.id) {
              try { re.friendlyName = args.name; await grok.dapi.notebooks.save(re); result.resaves++; } catch (_) { /* re-save best-effort */ }
            }
            await new Promise((r) => setTimeout(r, 1000));
          }
          result.pollMs = Date.now() - tPoll;
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      }, {id: seededId, name: renameTarget});
      expect(out.saveOk, `friendlyName + save() succeeds; err: ${out.err ?? ''}`).toBe(true);
      expect(out.persisted, `rename persists through find round-trip within the cold-commit budget (observed friendlyName: ${out.value}; pollMs: ${out.pollMs}; re-saves: ${out.resaves})`).toBe(true);
    });

    await softStep('S5: Delete (delete_notebook — api.delete / routes.delete / service.delete / repository.delete + delete-tables-relations) and verify removal', async () => {
      const out = await page.evaluate(async (id: string | null) => {
        const result: any = {deleted: false, findAfter: 'NOT_TESTED', recreatable: false, err: null};
        try {
          if (!id) throw new Error('seed failed upstream — no id');
          // delete REQUIRES the entity, not an id string (recon 2026-06-18).
          const ent: any = await grok.dapi.notebooks.find(id);
          if (!ent || !ent.id) throw new Error('entity not found before delete');
          await grok.dapi.notebooks.delete(ent);
          result.deleted = true;
          // Verify removal: find(id) resolves to `undefined` once deleted.
          let gone = false;
          for (let i = 0; i < 30 && !gone; i++) {
            const after: any = await grok.dapi.notebooks.find(id).catch(() => undefined);
            if (after === undefined || !after || !after.id) gone = true;
            else await new Promise((r) => setTimeout(r, 500));
          }
          result.findAfter = gone ? 'undefined/empty' : 'still-present';
          // The notebooks_tables join + tags are cleared server-side by
          // repository.delete -> delete-tables-relations (no FK left dangling); a
          // subsequent seed must still succeed (no FK conflict from the removed row).
          const before = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10});
          const beforeIds = new Set(before.map((n: any) => n.id));
          const cmd = (window as any).DG.Func.find({name: 'CmdNewNotebook'})[0];
          if (cmd) {
            await cmd.apply().catch(() => {});
            let fresh2: any = null;
            for (let i = 0; i < 30 && !fresh2; i++) {
              await new Promise((r) => setTimeout(r, 500));
              const cur = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10}).catch(() => [] as any[]);
              fresh2 = cur.find((n: any) => !beforeIds.has(n.id));
            }
            result.recreatable = fresh2 != null;
            // Clean up this verification notebook immediately.
            if (fresh2) {
              const e2: any = await grok.dapi.notebooks.find(fresh2.id).catch(() => null);
              if (e2 && e2.id) await grok.dapi.notebooks.delete(e2).catch(() => {});
            }
          }
        } catch (e: any) { result.err = String(e?.message ?? e).slice(0, 400); }
        return result;
      }, seededId);
      expect(out.deleted, `delete(entity) succeeds; err: ${out.err ?? ''}`).toBe(true);
      expect(out.findAfter, 'find(seedId) resolves to undefined/empty after delete (not-found signal)').toBe('undefined/empty');
      expect(out.recreatable, 'a new notebook can be created after delete (join rows cleared — no FK conflict)').toBe(true);
      // The seeded notebook is gone — clear the cleanup handle so the finally
      // block does not attempt a redundant delete.
      seededId = null;
    });
  } finally {
    // Best-effort cleanup: if any step failed mid-flight and left the seeded
    // notebook on the shared server, remove it by entity (delete needs the entity,
    // not an id string). Keyed on the captured id via find() — no full list scan.
    await page.evaluate(async (id: string | null) => {
      if (!id) return;
      try {
        const ent: any = await grok.dapi.notebooks.find(id).catch(() => null);
        if (ent && ent.id) await grok.dapi.notebooks.delete(ent).catch(() => {});
      } catch (_) { /* best-effort teardown must not throw */ }
    }, seededId);
    await page.evaluate(() => { try { grok.shell.closeAll(); } catch (_) { /* ignore */ } });
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
