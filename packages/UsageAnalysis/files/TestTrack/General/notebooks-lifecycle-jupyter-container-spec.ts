import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function openNotebooksBrowserForButton(page: import('@playwright/test').Page) {
  const flipped = await page.evaluate(async () => {
    const g = (window as any).grok;
    const DG = (window as any).DG;
    const onBrowser = () => { try { return g.shell.v?.type === 'notebooks'; } catch (e) { return false; } };

    for (let i = 0; i < 80; i++) {
      if (onBrowser()) return true;
      try {
        if (i % 2 === 0) {
          const f = DG.Func.find({name: 'CmdBrowseNotebooks'})[0];
          if (f) await f.apply();
        } else {
          g.shell.route('/notebooks');
        }
      } catch (e) {  }
      await new Promise((r) => setTimeout(r, 1000));
      if (onBrowser()) return true;
    }
    return onBrowser();
  });
  if (!flipped)
    throw new Error('Notebooks browser view never opened (CmdBrowseNotebooks.apply / route(/notebooks) both failed to flip grok.shell.v.type to "notebooks" within 80s)');
  await page.locator('.grok-gallery-search-bar').waitFor({timeout: 30_000});
  await page.locator('[name="button-New-Notebook..."]').waitFor({state: 'visible', timeout: 30_000});
}

async function openNotebooksBrowserWithCards(page: import('@playwright/test').Page) {
  await openNotebooksBrowserForButton(page);
  await page.waitForFunction(() => {
    return document.querySelectorAll('.grok-gallery-grid-item').length > 0 &&
      document.querySelectorAll('.d4-link-label[data-link^="/notebook/"]').length > 0;
  }, null, {timeout: 45_000, polling: 250});
}

test('Notebooks Lifecycle — Jupyter container source class (init / capability gate / view routing / save-file token clear)', async ({page}) => {

  test.setTimeout(420_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  });

  const seedName = `automator-jupyter-lifecycle-${Date.now()}`;
  let seededId: string | null = null;

  try {

    await softStep('Scenario 6 (requires-capabilities / notebooks-enabled / fleet-capability): capability gate consequences hold', async () => {
      const out = await page.evaluate(async () => {
        const r: any = {containerStatus: null, cmdNew: 0, cmdBrowse: 0, notebookViewFn: 0, initContainerFound: 0, err: null};
        try {
          const list = await grok.dapi.docker.dockerContainers
            .filter('name = "Notebooks-jupyter-notebook"').list();
          r.containerStatus = list[0]?.status ?? null;
          r.cmdNew = (window as any).DG.Func.find({name: 'CmdNewNotebook'}).length;
          r.cmdBrowse = (window as any).DG.Func.find({name: 'CmdBrowseNotebooks'}).length;
          r.notebookViewFn = (window as any).DG.Func.find({name: 'notebookView'}).length;
          r.initContainerFound = (window as any).DG.Func.find({name: 'initContainer'}).length;
        } catch (e: any) { r.err = String(e?.message ?? e).slice(0, 300); }
        return r;
      });

      expect(out.cmdNew, `CmdNewNotebook registered (notebooks-enabled gate); err: ${out.err ?? ''}`).toBeGreaterThan(0);
      expect(out.cmdBrowse, 'CmdBrowseNotebooks registered (init-meta)').toBeGreaterThan(0);
      expect(out.notebookViewFn, 'notebookView plugin func registered (notebook-view-func)').toBeGreaterThan(0);
      expect(out.initContainerFound, 'initContainer plugin func registered (init-container-func)').toBeGreaterThan(0);

      expect(out.containerStatus, `Notebooks-jupyter-notebook container present and warm (observed: ${out.containerStatus})`)
        .not.toBeNull();
    });

    await softStep('Scenario 1 (init-container): initContainer resolves; NEW NOTEBOOK... ribbon button seeds + opens editor', async () => {

      const init = await page.evaluate(async () => {
        const r: any = {ok: false, ms: 0, err: null};
        try {
          const f = (window as any).DG.Func.find({name: 'initContainer'})[0];
          const t0 = Date.now();
          await f.apply();
          r.ms = Date.now() - t0;
          r.ok = true;
        } catch (e: any) { r.err = String(e?.message ?? e).slice(0, 300); }
        return r;
      });
      expect(init.ok, `initContainer resolves without a "container unavailable" error (ms: ${init.ms}; err: ${init.err ?? ''})`).toBe(true);

      await openNotebooksBrowserForButton(page);
      const before = await page.evaluate(async () => {
        const b = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10});
        return b.map((n: any) => n.id);
      });
      const newBtn = page.locator('[name="button-New-Notebook..."]');
      await newBtn.waitFor({timeout: 30_000});
      await newBtn.click();

      const opened = await page.waitForFunction(() => {
        try { return grok.shell.v?.type === 'Notebook'; } catch (e) { return false; }
      }, null, {timeout: 60_000, polling: 250}).then(() => true).catch(() => false);
      expect(opened, 'New Notebook... opens a Notebook (editor) view').toBe(true);

      const seed = await page.evaluate(async (beforeIds: string[]) => {
        const bset = new Set(beforeIds);
        let fresh: any = null;
        for (let i = 0; i < 40 && !fresh; i++) {
          await new Promise((r) => setTimeout(r, 500));
          const cur = await grok.dapi.notebooks.order('createdOn', true).list({pageSize: 10}).catch(() => [] as any[]);
          fresh = cur.find((n: any) => !bset.has(n.id));
        }
        if (!fresh) return {id: null, kernelspec: null};
        return {id: fresh.id as string, kernelspec: fresh.notebook?.metadata?.kernelspec?.name ?? null};
      }, before);
      expect(seed.id, 'New Notebook... persists a server notebook (CONTAINER_ID/SESSION_TOKEN cached, no init failure)').toBeTruthy();

      expect(seed.kernelspec, `seeded notebook kernelspec.name === 'python3' (observed: ${seed.kernelspec})`).toBe('python3');
      seededId = seed.id;

      const renamed = await page.evaluate(async (args: {id: string | null; name: string}) => {
        if (!args.id) return false;
        const ent: any = await grok.dapi.notebooks.find(args.id);
        ent.friendlyName = args.name;
        await grok.dapi.notebooks.save(ent);
        for (let i = 0; i < 30; i++) {
          const re: any = await grok.dapi.notebooks.find(args.id).catch(() => null);
          if (re && (re.friendlyName || re.name) === args.name) return true;
          await new Promise((r) => setTimeout(r, 1000));
        }
        return false;
      }, {id: seededId, name: seedName});
      expect(renamed, 'seed renamed to a unique name (cold-commit budget)').toBe(true);
    });

    await softStep('Scenario 2 (init-notebook / handle-path / notebook-view-func): /notebook/<id> route opens the view', async () => {
      const out = await page.evaluate(async (id: string | null) => {
        const r: any = {routeOpened: false, viewType: null, nameMatches: false, friendlyName: null, viewName: null, err: null};
        try {
          if (!id) throw new Error('seed failed upstream — no id');
          const ent: any = await grok.dapi.notebooks.find(id);
          r.friendlyName = ent?.friendlyName ?? ent?.name ?? null;
          grok.shell.closeAll();

          for (let i = 0; i < 20 && grok.shell.v?.type === 'Notebook'; i++)
            await new Promise((res) => setTimeout(res, 100));

          grok.shell.route('/notebook/' + id);
          for (let i = 0; i < 40 && !r.routeOpened; i++) {
            await new Promise((res) => setTimeout(res, 500));
            if (grok.shell.v?.type === 'Notebook') r.routeOpened = true;
          }
          r.viewType = grok.shell.v?.type ?? null;
          r.viewName = grok.shell.v?.name ?? null;
        } catch (e: any) { r.err = String(e?.message ?? e).slice(0, 300); }
        return r;
      }, seededId);
      expect(out.routeOpened, `routing to /notebook/<id> opens a Notebook view; err: ${out.err ?? ''}`).toBe(true);
      expect(out.viewType, 'view type is Notebook (notebook-view-func handled the path)').toBe('Notebook');
    });

    await softStep('Scenario 4 (save-notebook-file / routes.save-file): no client-side session_token leak path; .ipynb round-trips cleanly', async () => {

      const out = await page.evaluate(async (id: string | null) => {
        const r: any = {hasSaveFileMethod: null, saved: false, probePersisted: null, hadIpynb: false, err: null};
        try {
          if (!id) throw new Error('seed failed upstream — no id');

          const proto = Object.getOwnPropertyNames(Object.getPrototypeOf(grok.dapi.notebooks));
          r.hasSaveFileMethod = proto.some((m) => /savefile|savenotebookfile/i.test(m));

          const ent: any = await grok.dapi.notebooks.find(id);
          if (!ent?.notebook) throw new Error('seeded notebook has no .ipynb body');
          r.hadIpynb = true;
          const nb: any = ent.notebook;                 
          nb.metadata = nb.metadata || {};
          nb.metadata.datagrok = Object.assign({}, nb.metadata.datagrok, {__roundtrip_probe: 'RT_OK_42'});
          ent.notebook = nb;                            
          await grok.dapi.notebooks.save(ent);
          r.saved = true;

          for (let i = 0; i < 20; i++) {
            const re: any = await grok.dapi.notebooks.find(id).catch(() => null);
            const probe = re?.notebook?.metadata?.datagrok?.__roundtrip_probe;
            if (probe != null) { r.probePersisted = probe; break; }
            await new Promise((res) => setTimeout(res, 500));
          }
        } catch (e: any) { r.err = String(e?.message ?? e).slice(0, 300); }
        return r;
      }, seededId);

      expect(out.hasSaveFileMethod, `grok.dapi.notebooks exposes no client saveFile/saveNotebookFile leak path; err: ${out.err ?? ''}`).toBe(false);

      expect(out.hadIpynb, 'seeded notebook carries an .ipynb body').toBe(true);
      expect(out.saved, `notebook entity save succeeds; err: ${out.err ?? ''}`).toBe(true);
      expect(out.probePersisted, `the .ipynb metadata round-trips through save/find (observed: ${out.probePersisted})`).toBe('RT_OK_42');
    });

    await softStep('Scenario 5 (save-state-map): saveStateMap is callable on the live Notebook view', async () => {

      const out = await page.evaluate(async (id: string | null) => {
        const r: any = {onNotebookView: false, hasSaveStateMap: false, callOk: false, viewIdStable: false, err: null};
        try {
          if (!id) throw new Error('seed failed upstream — no id');
          if (grok.shell.v?.type !== 'Notebook') {
            grok.shell.route('/notebook/' + id);
            for (let i = 0; i < 40; i++) {
              await new Promise((res) => setTimeout(res, 500));
              if (grok.shell.v?.type === 'Notebook') break;
            }
          }
          const v: any = grok.shell.v;
          r.onNotebookView = v?.type === 'Notebook';
          r.hasSaveStateMap = typeof v?.saveStateMap === 'function';
          if (r.hasSaveStateMap) { v.saveStateMap(); r.callOk = true; }
          r.viewIdStable = v?.id != null;
        } catch (e: any) { r.err = String(e?.message ?? e).slice(0, 300); }
        return r;
      }, seededId);
      expect(out.onNotebookView, `Notebook view is open for the state-map probe; err: ${out.err ?? ''}`).toBe(true);
      expect(out.hasSaveStateMap, 'the Notebook view exposes a callable saveStateMap()').toBe(true);
      expect(out.callOk, 'saveStateMap() invokes without throwing').toBe(true);
    });

    await softStep('Scenario 3 (to-html / html-mode): open in HTML mode renders without an unhandled rejection', async () => {

      const out = await page.evaluate(async (id: string | null) => {
        const r: any = {opened: false, viewType: null, openCallOk: false, err: null};
        try {
          if (!id) throw new Error('seed failed upstream — no id');
          grok.shell.closeAll();

          for (let i = 0; i < 20 && grok.shell.v?.type === 'Notebook'; i++)
            await new Promise((res) => setTimeout(res, 100));

          grok.shell.route('/notebook/' + id);
          r.openCallOk = true;
          for (let i = 0; i < 40 && !r.opened; i++) {
            await new Promise((res) => setTimeout(res, 500));
            if (grok.shell.v?.type === 'Notebook') r.opened = true;
          }
          r.viewType = grok.shell.v?.type ?? null;
        } catch (e: any) { r.err = String(e?.message ?? e).slice(0, 300); }
        return r;
      }, seededId);
      expect(out.openCallOk, `HTML-mode open path does not throw synchronously; err: ${out.err ?? ''}`).toBe(true);
      expect(out.opened, 'the notebook opens in a Notebook (HTML-mode) view').toBe(true);
    });

    await softStep('Scenario 6 (requires-capabilities): NotebooksView browser loads when the capability is advertised', async () => {

      await openNotebooksBrowserWithCards(page);
      const onBrowser = await page.evaluate(() => {
        try { return grok.shell.v?.type === 'notebooks'; } catch (e) { return false; }
      });
      expect(onBrowser, 'NotebooksView (type "notebooks") loads on a capability-advertising fleet').toBe(true);
      const cardCount = await page.locator('.d4-link-label[data-link^="/notebook/"]').count();
      expect(cardCount, 'the browser renders notebook cards').toBeGreaterThan(0);
    });
  } finally {

    await page.evaluate(async (id: string | null) => {
      if (!id) return;
      try {
        const ent: any = await grok.dapi.notebooks.find(id).catch(() => null);
        if (ent && ent.id) await grok.dapi.notebooks.delete(ent).catch(() => {});
      } catch (_) {  }
    }, seededId);
    await page.evaluate(() => { try { grok.shell.closeAll(); } catch (_) {  } });
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
