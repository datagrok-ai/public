/* ---
sub_features_covered: [bio.api.get-monomer-lib-helper, bio.lifecycle.init, bio.manage.libraries-view, bio.manage.monomer-collections-app]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {
  saveAllTablesWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';
test.use(specTestOptions);
test('Bio monomer_collection source-class lifecycle: write collection → reload via app → save project (collection FileShare entry survives) → reopen and verify', async ({page}) => {
  test.setTimeout(180_000);
  stepErrors.length = 0;
  const stamp = Date.now();
  const workingCollection = `bio-lifecycle-monomer-collection-${stamp}.json`;
  const workingCollectionPath = `System:AppData/Bio/monomer-collections/${workingCollection}`;
  const projectName = `bio-lifecycle-monomer-collection-project-${stamp}`;
  const syntheticSymbols = ['A', 'G', 'T', 'C'];
  const syntheticDescription = `bio-lifecycle proactive test collection (${stamp})`;
  const syntheticTags = ['bio-lifecycle-test', 'automation'];
  let saved: {projectId: string; primaryTableInfoId: string; tableInfoIds: string[]; layoutId: string | null} | null = null;
  let workingCollectionWritten = false;
  await loginToDatagrok(page);
  // Setup — open the HELM dataset so the renderer touches the monomer collection paths.
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 4000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
    if (hasMacromolecule) {
      for (let i = 0; i < 150; i++) {
        const cv: any = document.querySelector('[name="viewer-Grid"] canvas');
        if (cv && cv.width > 0 && cv.height > 0) break;
        await new Promise((r) => setTimeout(r, 200));
      }
    }
  }, 'System:AppData/Bio/tests/filter_HELM.csv');
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getMonomerLibHelper', 'Bio:getSeqHelper', 'Bio:getBioLib'];
    const deadline = Date.now() + 15_000;
    while (Date.now() < deadline) {
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
      }
      await new Promise((r) => setTimeout(r, 500));
    }
  });
  try {
    // Scenario 1 — Write a new collection via the shared lib-manager path
    await softStep('S1.1-1.3: getMonomerLibHelper returns singleton + writeCollection lands on FileShare with content round-trip', async () => {
      const result = await page.evaluate(async ({path, name, symbols, desc, tags}) => {
        const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        const hasHelper = helper != null;
        try {
          if (helper && typeof helper.awaitLoaded === 'function') {
            try { await helper.awaitLoaded(30_000); }
            catch (_) { try { await helper.awaitLoaded(); } catch (__) { /* non-fatal */ } }
          }
        } catch (_) { /* timeout non-fatal */ }
        let writeErr: string | null = null;
        try {
          if (typeof helper.addOrUpdateMonomerCollection === 'function')
            await helper.addOrUpdateMonomerCollection(name, symbols, desc, tags);
          else
            writeErr = 'helper.addOrUpdateMonomerCollection is not a function on the resolved singleton';
        } catch (e) {
          writeErr = String(e).slice(0, 200);
        }
        if (writeErr) return {hasHelper, writeErr, exists: false, readEqual: false};
        const exists = await grok.dapi.files.exists(path);
        let readBack: string | null = null;
        try { readBack = await grok.dapi.files.readAsText(path); }
        catch (e) { /* leave null */ }
        let parsed: any = null;
        try { parsed = readBack ? JSON.parse(readBack) : null; }
        catch (e) { /* leave null */ }
        const monomersMatch = parsed != null && Array.isArray(parsed.monomerSymbols) &&
          parsed.monomerSymbols.length === symbols.length &&
          symbols.every((s: string) => parsed.monomerSymbols.includes(s));
        return {
          hasHelper,
          writeErr: null,
          exists,
          readEqual: parsed != null,
          monomersMatch,
          parsedDescription: parsed?.description ?? null,
          parsedTags: Array.isArray(parsed?.tags) ? parsed.tags : null,
          parsedUpdatedBy: parsed?.updatedBy ?? null,
          parsedUpdatedOn: parsed?.updatedOn ?? null,
          fileLen: typeof readBack === 'string' ? readBack.length : 0,
        };
      }, {
        path: workingCollectionPath,
        name: workingCollection.replace(/\.json$/, ''),
        symbols: syntheticSymbols,
        desc: syntheticDescription,
        tags: syntheticTags,
      });
      expect(result.hasHelper).toBe(true);
      expect(result.writeErr, `addOrUpdateMonomerCollection error: ${result.writeErr}`).toBeNull();
      expect(result.exists,
        `expected ${workingCollectionPath} to exist on FileShare after addOrUpdateMonomerCollection`).toBe(true);
      expect(result.fileLen).toBeGreaterThan(0);
      expect(result.monomersMatch,
        `expected synthetic symbols [${syntheticSymbols.join(', ')}] in monomerSymbols round-trip`).toBe(true);
      expect(result.parsedDescription).toBe(syntheticDescription);
      expect(Array.isArray(result.parsedTags)).toBe(true);
      expect(result.parsedTags!.length).toBe(syntheticTags.length);
      expect((result.parsedUpdatedBy || '').length).toBeGreaterThan(0);
      expect((result.parsedUpdatedOn || '').length).toBeGreaterThan(0);
      workingCollectionWritten = true;
    });
    // Scenario 2 — Reload via the Monomer Collections app
    await softStep('S2.1-2.2: Monomer Collections app opens via registered function; working collection card appears in listing', async () => {
      const installDiag = await page.evaluate(async () => {
        const diag: any = {
          callOk: false,
          callErr: null as string | null,
          viewCaptured: false,
          viewNameAfterCall: null as string | null,
          addViewTried: false,
          addViewErr: null as string | null,
          setCurrentTried: false,
          setCurrentErr: null as string | null,
          viewsWalkFound: false,
          viewsWalkPromoted: false,
        };
        let view: any = null;
        try {
          view = await (grok as any).functions.call('Bio:monomerCollectionsApp', {});
          diag.callOk = true;
          diag.viewCaptured = view != null;
          try { diag.viewNameAfterCall = view?.name ?? null; } catch (_) { /* ignore */ }
        } catch (e: any) {
          diag.callErr = String(e).slice(0, 200);
        }
        if (view) {
          diag.addViewTried = true;
          try { (grok.shell as any).addView(view); }
          catch (e: any) { diag.addViewErr = String(e).slice(0, 200); }
          diag.setCurrentTried = true;
          try { (grok.shell as any).v = view; }
          catch (e: any) { diag.setCurrentErr = String(e).slice(0, 200); }
        }
        try {
          const views: any[] = Array.from((grok.shell as any).views || []);
          const target: any = views.find((v: any) => v?.name === 'Monomer Collections');
          if (target) {
            diag.viewsWalkFound = true;
            if ((grok.shell as any).v?.name !== 'Monomer Collections') {
              try { (grok.shell as any).v = target; diag.viewsWalkPromoted = true; }
              catch (_) { /* ignore */ }
            }
          }
        } catch (_) { /* ignore */ }
        return diag;
      });
      const installed = await page.waitForFunction(() => {
        try {
          const cur = (window as any).grok?.shell?.v?.name;
          if (cur === 'Monomer Collections') return true;
          const views: any[] = Array.from((window as any).grok?.shell?.views || []);
          return views.some((v: any) => v?.name === 'Monomer Collections');
        } catch (_) { return false; }
      }, null, {timeout: 30_000}).catch(() => null);
      // Settle: poll until loadCollections has built the working-collection card.
      const stem = workingCollection.replace(/\.json$/, '');
      await page.waitForFunction((s) => {
        try {
          const shell: any = (window as any).grok?.shell;
          let v: any = shell?.v;
          if (v?.name !== 'Monomer Collections') {
            const views: any[] = Array.from(shell?.views || []);
            v = views.find((x: any) => x?.name === 'Monomer Collections') || v;
          }
          const grid: any = v?.root?.querySelector?.('.monomer-collections-grid');
          if (!grid) return false;
          const cards = Array.from(grid.querySelectorAll('.monomer-collection-card'));
          return cards.some((c: any) => (c.dataset?.collectionName ?? '').includes(s));
        } catch (_) { return false; }
      }, stem, {timeout: 30_000}).catch(() => null);
      const installDiagStr = `installDiag=${JSON.stringify(installDiag)}`;
      const result = await page.evaluate((fileName) => {
        let v: any = grok.shell.v;
        if (v?.name !== 'Monomer Collections') {
          const views: any[] = Array.from((grok.shell as any).views || []);
          const found: any = views.find((x: any) => x?.name === 'Monomer Collections');
          if (found) v = found;
        }
        const root: any = v?.root;
        const viewBody: any = root?.querySelector?.('.monomer-collections-view');
        const grid: any = root?.querySelector?.('.monomer-collections-grid');
        const cards = grid ? Array.from(grid.querySelectorAll('.monomer-collection-card')) : [];
        const cardCollectionNames: string[] = cards.map((c: any) =>
          c.dataset?.collectionName ?? '');
        const stem = fileName.replace(/\.json$/, '');
        const hasWorkingCollection = cardCollectionNames.some((n: string) =>
          n === fileName || n === stem || n.includes(stem));
        const allViewNames: string[] = Array.from((grok.shell as any).views || [])
          .map((x: any) => x?.name ?? '<unnamed>');
        return {
          viewName: v?.name ?? null,
          currentViewName: (grok.shell as any).v?.name ?? null,
          allViewNames,
          rootPresent: !!root,
          viewBodyPresent: !!viewBody,
          gridPresent: !!grid,
          cardCount: cards.length,
          cardCollectionNames,
          hasWorkingCollection,
        };
      }, workingCollection);
      expect(result.viewName,
        `expected to locate a 'Monomer Collections' view (current or in shell.views); ${installDiagStr}; observed current=${result.currentViewName}, allViews=[${result.allViewNames.join(', ')}], installed-marker=${installed != null}`).toBe('Monomer Collections');
      expect(result.rootPresent).toBe(true);
      expect(result.viewBodyPresent).toBe(true);
      expect(result.gridPresent).toBe(true);
      expect(result.cardCount).toBeGreaterThanOrEqual(1);
      expect(result.hasWorkingCollection,
        `expected working collection '${workingCollection}' in card dataset names; observed: [${result.cardCollectionNames.join(', ')}]`).toBe(true);
    });
    // Pre-S2.3: promote the HELM TableView to current so the TableView-gated Bio top-menu re-mounts.
    await page.evaluate(async () => {
      const tvs: any[] = Array.from((grok.shell as any).tableViews || []);
      const helm = tvs.find((tv: any) => {
        const df = tv?.dataFrame;
        if (!df) return false;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.some((c: any) => c.semType === 'Macromolecule');
      });
      if (helm) {
        try { (grok.shell as any).v = helm; } catch (_) { /* setter read-only on some builds */ }
      }
      await new Promise((r) => setTimeout(r, 500));
    });
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 15_000});
    await softStep('S2.3: Bio | Manage | Monomer Libraries reachable as cross-surface entry point (top-menu DOM drive)', async () => {
      await page.evaluate(() => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      });
      await page.locator('[name="div-Bio---Manage"]').waitFor({state: 'attached', timeout: 15_000});
      await page.evaluate(() => {
        const manage = document.querySelector('[name="div-Bio---Manage"]')!;
        manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      });
      await page.locator('[name="div-Bio---Manage---Monomer-Libraries"]')
        .waitFor({state: 'attached', timeout: 15_000});
      await page.evaluate(() => {
        (document.querySelector(
          '[name="div-Bio---Manage---Monomer-Libraries"]') as HTMLElement).click();
      });
      await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Manage Monomer Libraries',
        null, {timeout: 30_000});
      const info = await page.evaluate(() => {
        const v = grok.shell.v;
        const root: any = v?.root;
        const form: any = root?.querySelector?.('.monomer-lib-controls-form');
        return {
          viewName: v?.name ?? null,
          formPresent: !!form,
        };
      });
      expect((info.viewName || '').toLowerCase()).toContain('monomer librar');
      expect(info.formPresent).toBe(true);
    });
    // Scenario 3 — Save project; verify the save/reopen path does not corrupt or
    // remove the collection's FileShare entry (no project→collection linkage API exists,
    // so the .md invariant is FileShare-entry survival, not a stored project reference).
    await softStep('S3.1: HELM dataset remains open; Macromolecule column renderer is dispatchable post-manage-view', async () => {
      const info = await page.evaluate(() => {
        const tvs: any[] = Array.from((grok.shell as any).tableViews || []);
        let chosen: any = null;
        let chosenMacro: any = null;
        for (const tv of tvs) {
          const df = tv?.dataFrame;
          if (!df) continue;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
          if (macro) {
            chosen = tv;
            chosenMacro = macro;
            break;
          }
        }
        if (!chosen) {
          return {
            hasDf: false,
            hasMacro: false,
            units: null,
            rowCount: 0,
            tableViewCount: tvs.length,
            viewNames: tvs.map((tv: any) => tv?.name ?? '<unnamed>'),
          };
        }
        const df = chosen.dataFrame;
        return {
          hasDf: true,
          hasMacro: !!chosenMacro,
          units: chosenMacro?.meta?.units ?? null,
          rowCount: df.rowCount,
          tableViewCount: tvs.length,
          viewNames: tvs.map((tv: any) => tv?.name ?? '<unnamed>'),
        };
      });
      expect(info.hasDf,
        `expected a TableView with a Macromolecule column among open table views; observed: count=${info.tableViewCount}, names=[${info.viewNames.join(', ')}]`).toBe(true);
      expect(info.hasMacro).toBe(true);
      expect(info.units).toBe('helm');
      expect(info.rowCount).toBeGreaterThan(0);
    });
    // Pre-S3.2: close auxiliary views + bring the HELM TableView forward for layout save.
    await page.evaluate(async () => {
      const views: any[] = Array.from(grok.shell.views || []);
      for (const name of ['Manage Monomer Libraries', 'Monomer Collections']) {
        const v: any = views.find((x: any) => x?.name === name);
        if (v && typeof v.close === 'function') {
          v.close();
          await new Promise((r) => setTimeout(r, 300));
        }
      }
      const tvs: any[] = Array.from((grok.shell as any).tableViews || []);
      const helm = tvs.find((tv: any) => {
        const df = tv?.dataFrame;
        if (!df) return false;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.some((c: any) => c.semType === 'Macromolecule');
      });
      if (helm && typeof (grok.shell as any).v !== 'undefined') {
        try { (grok.shell as any).v = helm; } catch (_) { /* setter may be read-only on some shell builds */ }
      }
      await new Promise((r) => setTimeout(r, 500));
    });
    await softStep('S3.2: save project with provenance (JS API path; mirrors monomer-library / macromolecule-column siblings)', async () => {
      const openTableCount = await page.evaluate(() => Array.from((grok.shell as any).tables).length);
      expect(openTableCount).toBeGreaterThan(0);
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.primaryTableInfoId).toBeTruthy();
      expect(saved.tableInfoIds.length,
        `expected every open shell table (${openTableCount}) persisted, not a partial multi-table save`).toBe(openTableCount);
    });
    await softStep('S3.3: reopen project — HELM survives + collection catalogue stable + Monomer Collections app shows working copy', async () => {
      if (!saved) throw new Error('S3.2 did not produce a saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      const post = await page.evaluate(async ({fileName}) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        let collectionNames: string[] = [];
        try {
          if (typeof helper?.listMonomerCollections === 'function')
            collectionNames = await helper.listMonomerCollections();
        } catch (_) { /* leave empty */ }
        const stem = fileName.replace(/\.json$/, '');
        const hasWorkingCollectionInList = collectionNames.some((n: string) =>
          n === fileName || n === stem || n.includes(stem));
        let view: any = null;
        try { view = await (grok as any).functions.call('Bio:monomerCollectionsApp', {}); }
        catch (_) { /* surface via viewName check below */ }
        if (view) {
          try { (grok.shell as any).addView(view); } catch (_) { /* ignore */ }
          try { (grok.shell as any).v = view; } catch (_) { /* ignore */ }
        }
        try {
          const views: any[] = Array.from((grok.shell as any).views || []);
          const target: any = views.find((vv: any) => vv?.name === 'Monomer Collections');
          if (target && (grok.shell as any).v?.name !== 'Monomer Collections')
            (grok.shell as any).v = target;
        } catch (_) { /* ignore */ }
        let installed = false;
        for (let i = 0; i < 30; i++) {
          try {
            const cur = (grok.shell as any).v?.name;
            if (cur === 'Monomer Collections') { installed = true; break; }
            const views: any[] = Array.from((grok.shell as any).views || []);
            if (views.some((vv: any) => vv?.name === 'Monomer Collections')) {
              installed = true;
              break;
            }
          } catch (_) { /* keep polling */ }
          await new Promise((r) => setTimeout(r, 200));
        }
        // Poll until the reopened Monomer Collections grid has built the working-collection card.
        for (let i = 0; i < 45; i++) {
          let vv: any = (grok.shell as any).v;
          if (vv?.name !== 'Monomer Collections') {
            const views: any[] = Array.from((grok.shell as any).views || []);
            vv = views.find((x: any) => x?.name === 'Monomer Collections') || vv;
          }
          const grid: any = vv?.root?.querySelector?.('.monomer-collections-grid');
          if (grid) {
            const cards = Array.from(grid.querySelectorAll('.monomer-collection-card'));
            if (cards.some((c: any) => (c.dataset?.collectionName ?? '').includes(stem))) break;
          }
          await new Promise((r) => setTimeout(r, 200));
        }
        let v: any = (grok.shell as any).v;
        if (v?.name !== 'Monomer Collections') {
          const views: any[] = Array.from((grok.shell as any).views || []);
          const found: any = views.find((x: any) => x?.name === 'Monomer Collections');
          if (found) v = found;
        }
        let cardCollectionNames: string[] = [];
        const viewName: string | null = v?.name ?? null;
        if (viewName === 'Monomer Collections') {
          const root: any = v?.root;
          const grid: any = root?.querySelector?.('.monomer-collections-grid');
          const cards = grid ? Array.from(grid.querySelectorAll('.monomer-collection-card')) : [];
          cardCollectionNames = cards.map((c: any) => c.dataset?.collectionName ?? '');
        }
        const hasWorkingCollectionInCards = cardCollectionNames.some((n: string) =>
          n === fileName || n === stem || n.includes(stem));
        const allViewNames: string[] = Array.from((grok.shell as any).views || [])
          .map((x: any) => x?.name ?? '<unnamed>');
        return {
          hasMacro: !!macro,
          units: macro?.meta?.units ?? null,
          helperResolved: helper != null,
          hasWorkingCollectionInList,
          collectionNamesCount: collectionNames.length,
          collectionNames,
          appViewName: viewName,
          currentViewName: (grok.shell as any).v?.name ?? null,
          allViewNames,
          installed,
          cardCount: cardCollectionNames.length,
          cardCollectionNames,
          hasWorkingCollectionInCards,
        };
      }, {fileName: workingCollection});
      expect(post.hasMacro).toBe(true);
      expect(post.units).toBe('helm');
      expect(post.helperResolved).toBe(true);
      expect(post.hasWorkingCollectionInList,
        `expected working collection '${workingCollection}' (or stem) in post-reopen listMonomerCollections(); observed: [${post.collectionNames.join(', ')}]`).toBe(true);
      expect(post.collectionNamesCount).toBeGreaterThanOrEqual(1);
      expect(post.appViewName,
        `expected to locate a 'Monomer Collections' view (current or in shell.views) post-reopen; observed currentView=${post.currentViewName}, installed=${post.installed}, allViews=[${post.allViewNames.join(', ')}]`).toBe('Monomer Collections');
      expect(post.cardCount).toBeGreaterThanOrEqual(1);
      expect(post.hasWorkingCollectionInCards,
        `expected working collection '${workingCollection}' (or stem) in post-reopen Monomer Collections app cards; observed: [${post.cardCollectionNames.join(', ')}]`).toBe(true);
    });
  } finally {
    // Scenario 4 — Cleanup (runs regardless of earlier failures)
    if (workingCollectionWritten) {
      await page.evaluate(async ({path, name}) => {
        try {
          const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
          if (helper && typeof helper.deleteMonomerCollection === 'function') {
            try { await helper.deleteMonomerCollection(name); }
            catch (_) { /* fall through to raw delete */ }
          }
          try {
            if (await grok.dapi.files.exists(path))
              await grok.dapi.files.delete(path);
          } catch (_) { /* best effort */ }
        } catch (_) { /* best effort */ }
      }, {path: workingCollectionPath, name: workingCollection.replace(/\.json$/, '')}).catch(() => {});
    }
    if (saved) {
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.primaryTableInfoId,
      });
    }
    await page.evaluate(async () => {
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      for (const d of dialogs)
        d.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      try {
        const views = Array.from(grok.shell.views || []);
        for (const name of ['Monomer Collections', 'Manage Monomer Libraries']) {
          const v: any = views.find((x: any) => x?.name === name);
          if (v && typeof v.close === 'function') v.close();
        }
      } catch (_) { /* best effort */ }
    }).catch(() => {});
  }
  finishSpec();
});
