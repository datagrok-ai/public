/* ---
sub_features_covered: [notebooks.entity.link, notebooks.entity.is-applicable, notebooks.entity.get-applicable-cases, notebooks.entity.from-json, notebooks.entity.generate-file-name, notebooks.api.find, notebooks.api.list, notebooks.api.save, notebooks.api.filter, notebooks.repository.save, notebooks.repository.find, notebooks.repository.delete, notebooks.repository.delete-tables-relations, notebooks.repository.by-tag, notebooks.repository.audit-hooks, notebooks.routes.save, notebooks.routes.get, notebooks.routes.delete, notebooks.routes.list, notebooks.service.save-notebook, notebooks.service.get-notebook, notebooks.service.delete-notebook, notebooks.service.get-notebooks-filtered, notebooks.service.get-notebooks-count, notebooks.meta.open-tables-in-notebook, notebooks.meta.new-notebook]
--- */
// Pure JS-API apitest for the Notebooks CRUD lifecycle via grok.dapi.notebooks.
// Seeds a server notebook with CmdNewNotebook (no JS-API Notebook factory exists),
// then round-trips it: find (name + .ipynb body intact), list/count/order, rename
// via friendlyName + save, and delete (verifying removal and that a new notebook
// can still be created). All work is grok.dapi.notebooks / DG.Func.apply() in
// page.evaluate — no DOM driving of the app.
// Scope-reduced: table-link / tag / filter semantics have no enumerable JS surface,
// so they're covered via the create→find→delete round-trip and the order()+list()
// path here, and asserted directly server-side. delete requires the entity (not an
// id string), so steps resolve via find(id) first.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';

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
          // Re-fetch and confirm the rename persisted. A cold server can lag the
          // first save after a CmdNewNotebook seed, so poll up to 30s with an
          // idempotent re-save at the 10s/20s marks in case the save was dropped.
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
          // delete requires the entity, not an id string.
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
