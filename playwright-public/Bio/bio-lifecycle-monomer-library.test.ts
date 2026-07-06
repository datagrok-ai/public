/* ---
sub_features_covered: [bio.api.get-monomer-lib-helper, bio.lifecycle.init, bio.manage.libraries-app, bio.manage.libraries-dialog, bio.manage.libraries-view]
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

// Opens ribbon Bio | Manage | Monomer Libraries, waiting on each submenu instead of blind sleeps.
async function openManageMonomerLibrariesRibbon(page: any): Promise<void> {
  await page.evaluate(() => (document.querySelector('[name="div-Bio"]') as HTMLElement).click());
  await page.locator('[name="div-Bio---Manage"]').waitFor({state: 'visible', timeout: 15_000});
  await page.evaluate(() => {
    const manage = document.querySelector('[name="div-Bio---Manage"]')!;
    manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
    manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
  });
  await page.locator('[name="div-Bio---Manage---Monomer-Libraries"]').waitFor({state: 'visible', timeout: 15_000});
  await page.evaluate(() =>
    (document.querySelector('[name="div-Bio---Manage---Monomer-Libraries"]') as HTMLElement).click());
  await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Manage Monomer Libraries',
    null, {timeout: 30_000});
  await page.locator('.monomer-lib-controls-form .ui-input-bool').first().waitFor({timeout: 15_000});
}

// Deterministic per-monomer color signature for the first N peptide symbols, read via the JS color API.
async function captureMonomerColors(page: any): Promise<{sig: string; count: number}> {
  return await page.evaluate(async () => {
    const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
    const lib: any = helper.getMonomerLib();
    const symbols: string[] = (typeof lib.getMonomerSymbolsByType === 'function'
      ? lib.getMonomerSymbolsByType('PEPTIDE') : []).slice().sort().slice(0, 20);
    const map: {[s: string]: string} = {};
    for (const s of symbols) map[s] = JSON.stringify(lib.getMonomerColors('HELM_AA', s));
    return {sig: JSON.stringify(map), count: symbols.length};
  });
}
test('Bio monomer_library source-class lifecycle: load → edit/save round-trip → save project with library reference', async ({page}) => {
  test.setTimeout(240_000);
  stepErrors.length = 0;
  const stamp = Date.now();
  const workingCopy = `bio-lifecycle-monomer-library-${stamp}.json`;
  const workingCopyPath = `System:AppData/Bio/monomer-libraries/${workingCopy}`;
  const projectName = `bio-lifecycle-monomer-library-project-${stamp}`;
  const canonicalLibCandidates = [
    'System:AppData/Bio/monomer-libraries/HELMCoreLibrary.json',
    'System:AppData/Bio/monomer-libraries/polytool-lib.json',
    'System:AppData/Bio/monomer-libraries/sample-lib.json',
  ];
  const syntheticSymbol = `XYZ_TEST_${stamp}`;
  let saved: {projectId: string; primaryTableInfoId: string; layoutId: string | null} | null = null;
  let workingCopyWritten = false;
  let preColors: {sig: string; count: number} | null = null;
  await loginToDatagrok(page);
  // Setup — open the HELM dataset so the renderer touches the library color-coding path.
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    const hasMacro = () => Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Macromolecule');
    for (let i = 0; i < 100; i++) {
      if (hasMacro()) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    if (hasMacro()) {
      for (let i = 0; i < 75; i++) {
        const cnv = document.querySelector('[name="viewer-Grid"] canvas') as HTMLCanvasElement | null;
        if (cnv && cnv.width > 0 && cnv.height > 0) break;
        await new Promise((r) => setTimeout(r, 200));
      }
    }
  }, 'System:AppData/Bio/tests/filter_HELM.csv');
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    for (let i = 0; i < 30; i++) {
      try {
        const helper = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        if (helper != null) return;
      } catch { /* retry */ }
      await new Promise((r) => setTimeout(r, 1000));
    }
  });
  try {
    // Scenario 1 — Load library via service surface
    await softStep('S1.1-1.2: getMonomerLibHelper returns singleton + canonical lib readable via FileShare', async () => {
      const result = await page.evaluate(async (candidates) => {
        const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        const helper2: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        const hasHelper = helper != null;
        const sameHelper = helper != null && helper === helper2;
        try {
          if (helper && typeof helper.awaitLoaded === 'function') {
            try { await helper.awaitLoaded(30_000); }
            catch (_) { try { await helper.awaitLoaded(); } catch (__) { /* non-fatal */ } }
          }
        } catch (_) { /* timeout is non-fatal here */ }
        let chosenPath: string | null = null;
        let sourceJson: string | null = null;
        let readErr: string | null = null;
        for (const p of candidates) {
          try {
            sourceJson = await grok.dapi.files.readAsText(p);
            chosenPath = p;
            break;
          } catch (e) {
            readErr = String(e).slice(0, 200);
            sourceJson = null;
          }
        }
        if (!sourceJson || !chosenPath)
          throw new Error(`S1.2: none of the canonical libraries resolved; last error: ${readErr}`);
        let parsed: any;
        try {
          parsed = JSON.parse(sourceJson);
        } catch (e) {
          throw new Error(`S1.2: ${chosenPath} did not parse as JSON: ${String(e).slice(0, 200)}`);
        }
        const monomers: any[] = Array.isArray(parsed) ? parsed : [];
        const allHaveSymbol = monomers.length > 0 && monomers.every((m: any) => typeof m?.symbol === 'string');
        const allHaveStructure = monomers.length > 0 &&
          monomers.every((m: any) => typeof m?.molfile === 'string' || typeof m?.smiles === 'string');
        return {
          hasHelper,
          sameHelper,
          chosenPath,
          monomersArrayPresent: Array.isArray(parsed),
          monomerCount: monomers.length,
          allHaveSymbol,
          allHaveStructure,
          sourceJsonLen: sourceJson.length,
        };
      }, canonicalLibCandidates);
      expect(result.hasHelper).toBe(true);
      expect(result.sameHelper, 'getMonomerLibHelper must return the same singleton instance on repeat calls').toBe(true);
      expect(result.chosenPath).not.toBeNull();
      expect(result.sourceJsonLen).toBeGreaterThan(0);
      expect(result.monomersArrayPresent).toBe(true);
      expect(result.monomerCount).toBeGreaterThan(0);
      expect(result.allHaveSymbol).toBe(true);
      expect(result.allHaveStructure).toBe(true);
    });
    await softStep('S1.3: ribbon Bio | Manage | Monomer Libraries opens View; canonical library appears in listing', async () => {
      await openManageMonomerLibrariesRibbon(page);
      const listing = await page.evaluate(() => {
        const v = grok.shell.v;
        const root: any = v?.root;
        const form: any = root?.querySelector?.('.monomer-lib-controls-form');
        const rows = form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [];
        const labels: string[] = rows.map((r: any) => {
          const span = r.querySelector('span');
          return span ? (span.textContent || '').trim() : '';
        });
        return {
          viewName: v?.name || null,
          formPresent: !!form,
          rowCount: rows.length,
          labels,
        };
      });
      expect((listing.viewName || '').toLowerCase()).toContain('monomer librar');
      expect(listing.formPresent).toBe(true);
      expect(listing.rowCount).toBeGreaterThanOrEqual(1);
      const canonicalStems = ['HELMCoreLibrary', 'polytool', 'sample-lib'];
      const hasCanonical = listing.labels.some((l: string) =>
        canonicalStems.some((stem) => l.toLowerCase().includes(stem.toLowerCase())));
      expect(hasCanonical, `expected one of [${canonicalStems.join(', ')}] in labels; observed: [${listing.labels.join(', ')}]`).toBe(true);
    });
    await softStep('S1.4: alternate Bio:manageMonomerLibraries dispatch yields a dialog with the same catalogue', async () => {
      const result = await page.evaluate(async () => {
        const root: any = grok.shell.v?.root;
        const viewForm: any = root?.querySelector?.('.monomer-lib-controls-form');
        const viewRows = viewForm ? Array.from(viewForm.querySelectorAll('.ui-input-bool')) : [];
        const viewLabels: string[] = viewRows.map((r: any) => {
          const span = r.querySelector('span');
          return span ? (span.textContent || '').trim() : '';
        }).filter((s: string) => s.length > 0);
        try {
          // Fire-and-forget — awaiting would deadlock (promise resolves on dialog close).
          (grok as any).functions.call('Bio:manageMonomerLibraries', {}).catch(() => {});
        } catch (e) {
          return {dispatchErr: String(e).slice(0, 200), viewLabels};
        }
        let dialog: Element | null = null;
        for (let i = 0; i < 50; i++) {
          const candidates = Array.from(document.querySelectorAll('.d4-dialog'));
          for (const d of candidates) {
            if (d.querySelector('.monomer-lib-controls-form')) {
              dialog = d;
              break;
            }
          }
          if (dialog) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        let dialogLabels: string[] = [];
        let dialogRowCount = 0;
        if (dialog) {
          const form: any = dialog.querySelector('.monomer-lib-controls-form');
          const rows = form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [];
          dialogRowCount = rows.length;
          dialogLabels = rows.map((r: any) => {
            const span = r.querySelector('span');
            return span ? (span.textContent || '').trim() : '';
          }).filter((s: string) => s.length > 0);
          dialog.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
          for (let i = 0; i < 40; i++) {
            if (!document.querySelector('.d4-dialog .monomer-lib-controls-form')) break;
            await new Promise((r) => setTimeout(r, 100));
          }
        }
        return {
          dispatchErr: null,
          dialogOpened: !!dialog,
          viewRowCount: viewLabels.length,
          dialogRowCount,
          viewLabels,
          dialogLabels,
          cataloguesAgree:
            viewLabels.length === dialogLabels.length &&
            viewLabels.every((l: string) => dialogLabels.includes(l)),
        };
      });
      expect(result.dispatchErr, `dispatch error: ${result.dispatchErr}`).toBeNull();
      expect(result.dialogOpened,
        `expected dialog with .monomer-lib-controls-form within 10s; view labels: [${result.viewLabels.join(', ')}]`).toBe(true);
      expect(result.dialogRowCount).toBeGreaterThanOrEqual(1);
      expect(result.cataloguesAgree,
        `view labels [${result.viewLabels.join(', ')}] disagree with dialog labels [${result.dialogLabels.join(', ')}]`).toBe(true);
    });
    // Scenario 2 — Save edited library back to FileShare
    await softStep('S2.1: working copy lands under System:AppData/Bio/monomer-libraries via writeAsText', async () => {
      const result = await page.evaluate(async ({src, dst}) => {
        const sourceJson = await grok.dapi.files.readAsText(src);
        await grok.dapi.files.writeAsText(dst, sourceJson);
        const readBack = await grok.dapi.files.readAsText(dst);
        return {
          srcLen: sourceJson.length,
          dstLen: readBack.length,
          contentEqual: readBack === sourceJson,
        };
      }, {src: canonicalLibCandidates[0], dst: workingCopyPath}).catch(async () => {
        return await page.evaluate(async ({candidates, dst}) => {
          let sourceJson: string | null = null;
          for (const p of candidates) {
            try {
              sourceJson = await grok.dapi.files.readAsText(p);
              break;
            } catch (_) { /* try next */ }
          }
          if (!sourceJson) throw new Error('S2.1: no canonical library readable for working-copy seed');
          await grok.dapi.files.writeAsText(dst, sourceJson);
          const readBack = await grok.dapi.files.readAsText(dst);
          return {
            srcLen: sourceJson.length,
            dstLen: readBack.length,
            contentEqual: readBack === sourceJson,
          };
        }, {candidates: canonicalLibCandidates, dst: workingCopyPath});
      });
      expect(result.srcLen).toBeGreaterThan(0);
      expect(result.dstLen).toBe(result.srcLen);
      expect(result.contentEqual).toBe(true);
      workingCopyWritten = true;
    });
    await softStep('S2.2: edit JSON in memory + write back; synthetic monomer round-trips through FileShare', async () => {
      const result = await page.evaluate(async ({path, symbol}) => {
        const before = await grok.dapi.files.readAsText(path);
        const parsed: any = JSON.parse(before);
        const monomers: any[] = Array.isArray(parsed) ? parsed : [];
        const beforeCount = monomers.length;
        const synthetic: any = {
          symbol,
          name: symbol,
          molfile: '\n     RDKit          2D\n\n  7  6  0  0  0  0  0  0  0  0999 V2000\n    1.6702    1.3929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3712    0.6429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9279    1.3929    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.2269    0.6429    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    0.3712   -0.8571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9279   -1.6071    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6702   -1.6071    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1\n  2  3  1  0\n  3  4  1  0\n  2  5  1  0\n  5  6  2  0\n  5  7  1  0\nM  RGP  2   4   1   7   2\nM  END\n',
          smiles: 'C[C@H](N[*:1])C(=O)[*:2]',
          polymerType: 'PEPTIDE',
          monomerType: 'Backbone',
          naturalAnalog: 'X',
          id: 0,
          rgroups: [
            {
              alternateId: 'R1-H',
              capGroupName: 'H',
              capGroupSMILES: '[*:1][H]',
              label: 'R1',
            },
            {
              alternateId: 'R2-OH',
              capGroupName: 'OH',
              capGroupSMILES: 'O[*:2]',
              label: 'R2',
            },
          ],
        };
        monomers.push(synthetic);
        const after = JSON.stringify(monomers);
        await grok.dapi.files.writeAsText(path, after);
        const readBack = await grok.dapi.files.readAsText(path);
        const reparsed: any = JSON.parse(readBack);
        const reMonomers: any[] = Array.isArray(reparsed) ? reparsed : [];
        const hasSynthetic = reMonomers.some((m: any) => m?.symbol === symbol);
        return {
          beforeCount,
          afterCount: reMonomers.length,
          hasSynthetic,
        };
      }, {path: workingCopyPath, symbol: syntheticSymbol});
      expect(result.afterCount).toBe(result.beforeCount + 1);
      expect(result.hasSynthetic).toBe(true);
    });
    await softStep('S2.3: reload via loadMonomerLib(true) — working copy available AND synthetic monomer in in-memory cache', async () => {
      const result = await page.evaluate(async ({fileName, stem, symbol}) => {
        const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        try {
          // Refresh provider file lists first so the freshly written working copy is picked up by the reload.
          if (typeof helper.getAvaliableLibraryNames === 'function')
            await helper.getAvaliableLibraryNames(true);
          if (typeof helper.loadMonomerLib === 'function')
            await helper.loadMonomerLib(true);
          else if (typeof helper.loadLibraries === 'function')
            await helper.loadLibraries(true);
        } catch (e) {
          return {reloadErr: String(e).slice(0, 200), inAvailable: false, inMemory: false, availableNames: [], monomerCount: 0};
        }
        try {
          if (typeof helper.awaitLoaded === 'function') {
            try { await helper.awaitLoaded(20_000); }
            catch (_) { try { await helper.awaitLoaded(); } catch (__) { /* ignore */ } }
          }
        } catch (_) { /* timeout non-fatal */ }
        let availableNames: string[] = [];
        try {
          try { availableNames = await helper.getAvaliableLibraryNames(true); }
          catch (_) { availableNames = await helper.getAvaliableLibraryNames(); }
        } catch (e) {
          return {reloadErr: String(e).slice(0, 200), inAvailable: false, inMemory: false, availableNames: [], monomerCount: 0};
        }
        const inAvailable = availableNames.some((n: string) => n === fileName || n.includes(stem));
        let inMemory = false;
        let monomerCount = 0;
        try {
          const lib: any = helper.getMonomerLib();
          const all: any[] = typeof lib.toJSON === 'function' ? lib.toJSON() : [];
          monomerCount = all.length;
          inMemory = all.some((m: any) => m?.symbol === symbol) ||
            (typeof lib.getMonomer === 'function' && lib.getMonomer('PEPTIDE', symbol) != null);
        } catch (_) { /* leave false */ }
        return {reloadErr: null, inAvailable, inMemory, availableNames, monomerCount};
      }, {fileName: workingCopy, stem: workingCopy.replace(/\.json$/, ''), symbol: syntheticSymbol});
      expect(result.reloadErr, `reload error: ${result.reloadErr}`).toBeNull();
      expect(result.inAvailable,
        `expected working copy '${workingCopy}' (or stem) in available libraries; observed: [${result.availableNames.join(', ')}]`).toBe(true);
      expect(result.monomerCount).toBeGreaterThan(0);
      expect(result.inMemory,
        `expected synthetic monomer '${syntheticSymbol}' present in reloaded in-memory library cache (${result.monomerCount} monomers) — stale-singleton regression`).toBe(true);
    });
    await softStep('S2.4: reopen Manage Monomer Libraries view — working copy appears in listing', async () => {
      await page.evaluate(() => {
        const v = grok.shell.v;
        if (v && v.name === 'Manage Monomer Libraries' && typeof v.close === 'function') v.close();
      });
      await openManageMonomerLibrariesRibbon(page);
      const result = await page.evaluate(({fileName, stem}) => {
        const root: any = grok.shell.v?.root;
        const form: any = root?.querySelector?.('.monomer-lib-controls-form');
        const rows = form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [];
        const labels: string[] = rows.map((r: any) => {
          const span = r.querySelector('span');
          return span ? (span.textContent || '').trim() : '';
        });
        const hasWorkingCopy = labels.some((l: string) =>
          l === fileName || l.includes(stem));
        return {
          formPresent: !!form,
          rowCount: rows.length,
          labels,
          hasWorkingCopy,
        };
      }, {fileName: workingCopy, stem: workingCopy.replace(/\.json$/, '')});
      expect(result.formPresent).toBe(true);
      expect(result.rowCount).toBeGreaterThanOrEqual(1);
      expect(result.hasWorkingCopy,
        `expected working copy '${workingCopy}' (or stem) in manage-view labels; observed: [${result.labels.join(', ')}]`).toBe(true);
    });
    // Scenario 3 — Save project that references the library
    await softStep('S3.1: HELM dataset remains open; Macromolecule column renderer is dispatchable', async () => {
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
      // Snapshot the monomer color set from the working-copy library before the save (scenario headline invariant).
      preColors = await captureMonomerColors(page);
      expect(preColors.count, 'expected peptide monomer colors captured before save').toBeGreaterThan(0);
    });
    // Close the Manage view so the layout-save helper sees the HELM TableView.
    await page.evaluate(async () => {
      const views: any[] = Array.from(grok.shell.views || []);
      const manage: any = views.find((v: any) => v?.name === 'Manage Monomer Libraries');
      if (manage && typeof manage.close === 'function') {
        manage.close();
        await new Promise((r) => setTimeout(r, 500));
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
    await softStep('S3.2: save project with provenance (JS API path; mirrors macromolecule-column sibling)', async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.primaryTableInfoId).toBeTruthy();
    });
    await softStep('S3.3: reopen project — HELM survives + library catalogue stable + Macromolecule semType holds', async () => {
      if (!saved) throw new Error('S3.2 did not produce a saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      const post = await page.evaluate(async ({fileName, stem}) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        let availableNames: string[] = [];
        try {
          if (typeof helper?.getAvaliableLibraryNames === 'function') {
            try { availableNames = await helper.getAvaliableLibraryNames(true); }
            catch (_) { availableNames = await helper.getAvaliableLibraryNames(); }
          }
        } catch (_) { /* leave empty */ }
        const hasWorkingCopy = availableNames.some((n: string) =>
          n === fileName || n.includes(stem));
        return {
          hasMacro: !!macro,
          units: macro?.meta?.units ?? null,
          helperResolved: helper != null,
          hasWorkingCopy,
          availableNames,
          availableCount: availableNames.length,
        };
      }, {fileName: workingCopy, stem: workingCopy.replace(/\.json$/, '')});
      expect(post.hasMacro).toBe(true);
      expect(post.units).toBe('helm');
      expect(post.helperResolved).toBe(true);
      expect(post.hasWorkingCopy,
        `expected working copy '${workingCopy}' (or stem) in post-reopen available libraries; observed: [${post.availableNames.join(', ')}]`).toBe(true);
      expect(post.availableCount).toBeGreaterThanOrEqual(1);
      // Color-stability: the monomer color set must be unchanged across save/reopen (no color reset on reload).
      const postColors = await captureMonomerColors(page);
      expect(postColors.count, 'expected peptide monomer colors after reopen').toBeGreaterThan(0);
      expect(postColors.sig,
        'Macromolecule monomer color set changed across save/reopen — color-stability regression')
        .toBe(preColors!.sig);
    });
  } finally {
    // Scenario 4 — Cleanup (runs regardless of earlier failures)
    if (workingCopyWritten) {
      await page.evaluate(async (p) => {
        try { await grok.dapi.files.delete(p); } catch (_) { /* best effort */ }
      }, workingCopyPath).catch(() => {});
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
        const manageView: any = views.find((v: any) => v?.name === 'Manage Monomer Libraries');
        if (manageView && typeof manageView.close === 'function') manageView.close();
      } catch (_) { /* best effort */ }
    }).catch(() => {});
  }
  finishSpec();
});
