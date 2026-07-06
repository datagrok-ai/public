/* ---
sub_features_covered: [biostructure.file-open.importWithNgl, biostructure.file-preview.ngl-density, biostructure.file-preview.ngl-structure, biostructure.file-preview.ngl-surface, biostructure.grid-context-menu.show-ngl-viewer, biostructure.ngl-viewer, biostructure.ngl-viewer.props, biostructure.panel.pdb-id-ngl]
--- */
// NGL viewer extension: mount + Style/Data/Behaviour props + file-handler routing + grid context menu
// + PDB id panel. NGL representation uses 'ball+stick' (Mol* uses 'ball-and-stick'). File-routing
// (Scenarios 2-3) is asserted via DG.Func.find registry probes (no NGL-format fixtures on disk).
// Canvas pixels are not asserted; mount + props round-trip + menu labels + widget roots are.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const sample1bdq = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

test('BiostructureViewer / NGL viewer extension (mount + props + file-routing + grid-context + PDB id panel)', async ({page}) => {
  test.setTimeout(180_000);
  stepErrors.length = 0;

  const pageErrors: string[] = [];
  page.on('pageerror', (err) => { pageErrors.push(err.message); });

  await loginToDatagrok(page);

  // Baseline environment setup.
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  });

  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  // Prerequisite guard: the NGL viewer + PDB-id NGL panel live in BiostructureViewer. If the
  // NGL extension isn't registered on this deployment, skip fast with a clear message instead
  // of burning the full 180s budget on a viewer/pane that can never mount.
  const nglReady = await page.evaluate(() => {
    const has = (name: string) => (DG.Func.find({package: 'BiostructureViewer', name}) || []).length > 0;
    return has('pdbIdNglPanelWidget') && has('importWithNgl');
  });
  test.skip(!nglReady, 'BiostructureViewer NGL extension not registered on this deployment');

  try {
    // SCENARIO 1 — NGL viewer add via tv.addViewer('NGL') + Style/Data/Behaviour props round-trip.
    let scenario1Mounted = false;

    await softStep('Scenario 1 step 1-3 — Open Molecule3D table; tv.addViewer("NGL"); canvas mount + container DOM', async () => {
      const res = await page.evaluate(async (pdbPath) => {
        const pollUntil = async (pred: () => boolean, timeoutMs = 3000, stepMs = 75) => {
          const t0 = Date.now();
          while (Date.now() - t0 < timeoutMs) {
            try { if (pred()) return true; } catch (_) { /* predicate not ready */ }
            await new Promise((r) => setTimeout(r, stepMs));
          }
          return false;
        };
        grok.shell.closeAll();
        const content = await grok.dapi.files.readAsText(pdbPath);
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['1bdq']),
          DG.Column.fromStrings('structure', [content]),
        ]);
        const col = df.col('structure');
        col.semType = 'Molecule3D';
        try { col.setTag('cell.renderer', 'Molecule3D'); } catch (_) { /* tag set */ }
        df.name = 'ngl-extension-fixture';
        const tv = grok.shell.addTableView(df);
        await pollUntil(() => !!document.querySelector('[name="viewer-Grid"]'));
        const v = tv.addViewer('NGL');
        await pollUntil(() => !!document.querySelector('[name="viewer-NGL"]') &&
          v?.props?.get?.('representation') === 'cartoon', 30_000);
        return {
          rowCount: df.rowCount,
          hasNglContainer: !!document.querySelector('[name="viewer-NGL"]'),
          vType: v?.type,
          defaultRep: v?.props?.get?.('representation') ?? null,
          viewerTypes: tv.viewers ? Array.from(tv.viewers).map((x: any) => x.type) : [],
        };
      }, sample1bdq);

      await expect(page.locator('[name="viewer-NGL"]')).toBeVisible({timeout: 30_000});

      expect(res.hasNglContainer).toBe(true);
      expect(res.vType).toBe('NGL');
      expect(res.defaultRep).toBe('cartoon');
      expect(res.viewerTypes).toContain('NGL');
      scenario1Mounted = true;
    });

    await softStep('Scenario 1 step 5 — Style: representation cartoon -> ball+stick (NGL choice; sibling Mol* differs)', async () => {
      if (!scenario1Mounted) return;
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'NGL') as any : null;
        if (!v) return {ok: false};
        const pollProp = async (prop: string, expected: any) => {
          const t0 = Date.now();
          while (Date.now() - t0 < 3000) {
            if (v.props.get(prop) === expected) return true;
            await new Promise((r) => setTimeout(r, 75));
          }
          return false;
        };
        // NGL uses 'ball+stick' (Mol* uses 'ball-and-stick').
        v.setOptions({representation: 'ball+stick'});
        await pollProp('representation', 'ball+stick');
        const afterRep = v.props.get('representation');
        v.setOptions({representation: 'cartoon'});
        await pollProp('representation', 'cartoon');
        const restoredRep = v.props.get('representation');
        return {ok: true, afterRep, restoredRep};
      });
      expect(res.ok).toBe(true);
      expect(res.afterRep).toBe('ball+stick');
      expect(res.restoredRep).toBe('cartoon');
    });

    await softStep('Scenario 1 step 6 — Data: set ligandColumnName to the Molecule3D column', async () => {
      if (!scenario1Mounted) return;
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'NGL') as any : null;
        if (!v) return {ok: false};
        v.setOptions({ligandColumnName: 'structure'});
        const t0 = Date.now();
        while (Date.now() - t0 < 3000 && v.props.get('ligandColumnName') !== 'structure')
          await new Promise((r) => setTimeout(r, 75));
        const ligandColAfter = v.props.get('ligandColumnName');
        return {ok: true, ligandColAfter};
      });
      expect(res.ok).toBe(true);
      expect(res.ligandColAfter).toBe('structure');
    });

    await softStep('Scenario 1 steps 7-8 — Behaviour: toggle showCurrentRowLigand + showMouseOverRowLigand round-trip', async () => {
      if (!scenario1Mounted) return;
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'NGL') as any : null;
        if (!v) return {ok: false};
        const pollProp = async (prop: string, expected: any) => {
          const t0 = Date.now();
          while (Date.now() - t0 < 3000) {
            if (v.props.get(prop) === expected) return true;
            await new Promise((r) => setTimeout(r, 75));
          }
          return false;
        };
        // Defaults: showCurrentRowLigand=true, showMouseOverRowLigand=true.
        const initCurrent = v.props.get('showCurrentRowLigand');
        const initMouseOver = v.props.get('showMouseOverRowLigand');
        v.setOptions({showCurrentRowLigand: false});
        await pollProp('showCurrentRowLigand', false);
        const curOff = v.props.get('showCurrentRowLigand');
        v.setOptions({showCurrentRowLigand: true});
        await pollProp('showCurrentRowLigand', true);
        const curOn = v.props.get('showCurrentRowLigand');
        v.setOptions({showMouseOverRowLigand: false});
        await pollProp('showMouseOverRowLigand', false);
        const moOff = v.props.get('showMouseOverRowLigand');
        v.setOptions({showMouseOverRowLigand: true});
        await pollProp('showMouseOverRowLigand', true);
        const moOn = v.props.get('showMouseOverRowLigand');
        return {ok: true, initCurrent, initMouseOver, curOff, curOn, moOff, moOn};
      });
      expect(res.ok).toBe(true);
      expect(res.initCurrent).toBe(true);
      expect(res.initMouseOver).toBe(true);
      expect(res.curOff).toBe(false);
      expect(res.curOn).toBe(true);
      expect(res.moOff).toBe(false);
      expect(res.moOn).toBe(true);
    });

    await softStep('Scenario 1 step 9 — Behaviour: showSelectedRowsLigands + select rows via dataframe selection', async () => {
      if (!scenario1Mounted) return;
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'NGL') as any : null;
        if (!v) return {ok: false};
        const df = tv.dataFrame;
        // Default: showSelectedRowsLigands=false.
        const initSelected = v.props.get('showSelectedRowsLigands');
        v.setOptions({showSelectedRowsLigands: true});
        const t0 = Date.now();
        while (Date.now() - t0 < 3000 && v.props.get('showSelectedRowsLigands') !== true)
          await new Promise((r) => setTimeout(r, 75));
        const selOn = v.props.get('showSelectedRowsLigands');
        // 1-row fixture: assert the property round-trip + selection driver (multi-row math is GROK-17967's).
        df.selection.init((i: number) => i === 0);
        const t1 = Date.now();
        while (Date.now() - t1 < 3000 && df.selection.trueCount < 1)
          await new Promise((r) => setTimeout(r, 75));
        const selectedCount = df.selection.trueCount;
        return {ok: true, initSelected, selOn, selectedCount, rowCount: df.rowCount};
      });
      expect(res.ok).toBe(true);
      expect(res.initSelected).toBe(false);
      expect(res.selOn).toBe(true);
      expect(res.selectedCount).toBeGreaterThanOrEqual(1);
      expect(res.selectedCount).toBeLessThanOrEqual(res.rowCount);
    });

    // SCENARIO 2 — importWithNgl is the registered handler for Mol*-incapable extensions (registry probe).
    await softStep('Scenario 2 — importWithNgl is the registered handler for mmtf/cns/prmtop/ccp4 (registry probe)', async () => {
      const res = await page.evaluate(() => {
        const importWithNglFns = DG.Func.find({name: 'importWithNgl', package: 'BiostructureViewer'});
        const importPdbFns = DG.Func.find({name: 'importPdb', package: 'BiostructureViewer'});
        const importPdbqtFns = DG.Func.find({name: 'importPdbqt', package: 'BiostructureViewer'});
        const importXYZFns = DG.Func.find({name: 'importXYZ', package: 'BiostructureViewer'});
        const importWithNgl = importWithNglFns && importWithNglFns[0];
        const importPdb = importPdbFns && importPdbFns[0];
        const importPdbqt = importPdbqtFns && importPdbqtFns[0];
        const importXYZ = importXYZFns && importXYZFns[0];
        // Normalize the ext option (comma-separated with optional spaces) to
        // a set of clean extensions.
        const extSet = (fn: any): Set<string> => {
          const raw = (fn?.options?.ext || '') as string;
          return new Set(
            raw.split(',').map((x: string) => x.trim().toLowerCase()).filter((x: string) => x.length > 0),
          );
        };
        const nglExt = extSet(importWithNgl);
        const pdbExt = extSet(importPdb);
        const pdbqtExt = extSet(importPdbqt);
        const xyzExt = extSet(importXYZ);
        // Each Mol*-incapable extension must be in importWithNgl's ext set and no other importer's.
        const nglOnlyExts = ['mmtf', 'cns', 'top', 'prmtop', 'ply', 'obj', 'ccp4'];
        const routingResults: Record<string, {inNgl: boolean; inPdb: boolean; inPdbqt: boolean; inXyz: boolean}> = {};
        for (const ext of nglOnlyExts) {
          routingResults[ext] = {
            inNgl: nglExt.has(ext),
            inPdb: pdbExt.has(ext),
            inPdbqt: pdbqtExt.has(ext),
            inXyz: xyzExt.has(ext),
          };
        }
        return {
          importWithNglRegistered: !!importWithNgl,
          importWithNglInputCount: importWithNgl?.inputs?.length ?? null,
          importWithNglFirstInputName: importWithNgl?.inputs?.[0]?.name ?? null,
          importWithNglFirstInputType: importWithNgl?.inputs?.[0]?.propertyType ?? null,
          nglExtList: Array.from(nglExt).sort(),
          routingResults,
        };
      });

      // importWithNgl must be registered with the canonical fileHandler signature.
      expect(res.importWithNglRegistered).toBe(true);
      expect(res.importWithNglInputCount).toBe(1);
      expect(res.importWithNglFirstInputName).toBe('fileContent');
      expect(res.importWithNglFirstInputType).toBe('string');

      // Each extension must route only to importWithNgl (no collision with importPdb/Pdbqt/XYZ).
      for (const ext of ['mmtf', 'cns', 'prmtop', 'ccp4']) {
        const r = res.routingResults[ext];
        expect(
          r.inNgl,
          `Expected importWithNgl to register extension '${ext}' but ext set is ${JSON.stringify(res.nglExtList)}`,
        ).toBe(true);
        expect(
          r.inPdb,
          `Routing collision: extension '${ext}' is also in importPdb's ext set (GROK-14442-shape regression).`,
        ).toBe(false);
        expect(
          r.inPdbqt,
          `Routing collision: extension '${ext}' is also in importPdbqt's ext set.`,
        ).toBe(false);
        expect(
          r.inXyz,
          `Routing collision: extension '${ext}' is also in importXYZ's ext set.`,
        ).toBe(false);
      }
    });

    // SCENARIO 3 — NGL preview file-viewers (structure/surface/density) registered with expected exts.
    await softStep('Scenario 3 — NGL preview file-viewers are registered with the expected ext sets (registry probe)', async () => {
      const res = await page.evaluate(() => {
        const previewNglStructureFns = DG.Func.find({name: 'previewNglStructure', package: 'BiostructureViewer'});
        const previewNglSurfaceFns = DG.Func.find({name: 'previewNglSurface', package: 'BiostructureViewer'});
        const previewNglDensityFns = DG.Func.find({name: 'previewNglDensity', package: 'BiostructureViewer'});
        const previewMolstarStructureFns = DG.Func.find({name: 'previewBiostructureStructure', package: 'BiostructureViewer'});
        const previewMolstarDensityFns = DG.Func.find({name: 'previewBiostructureDensity', package: 'BiostructureViewer'});

        const fn = (arr: any[]): any => arr && arr[0];
        // A fileViewer's match surface is options.fileViewer — probe it exactly (no ext fallback,
        // so a preview func registered without a fileViewer binding fails rather than passes on ext).
        const extSet = (f: any): Set<string> => {
          const raw = (f?.options?.fileViewer || '') as string;
          return new Set(
            raw.split(',').map((x: string) => x.trim().toLowerCase()).filter((x: string) => x.length > 0),
          );
        };

        const structureFn = fn(previewNglStructureFns);
        const surfaceFn = fn(previewNglSurfaceFns);
        const densityFn = fn(previewNglDensityFns);
        const molstarStructureFn = fn(previewMolstarStructureFns);
        const molstarDensityFn = fn(previewMolstarDensityFns);

        return {
          structureRegistered: !!structureFn,
          surfaceRegistered: !!surfaceFn,
          densityRegistered: !!densityFn,
          structureExt: Array.from(extSet(structureFn)).sort(),
          surfaceExt: Array.from(extSet(surfaceFn)).sort(),
          densityExt: Array.from(extSet(densityFn)).sort(),
          molstarStructureExt: Array.from(extSet(molstarStructureFn)).sort(),
          molstarDensityExt: Array.from(extSet(molstarDensityFn)).sort(),
        };
      });

      expect(res.structureRegistered).toBe(true);
      expect(res.surfaceRegistered).toBe(true);
      expect(res.densityRegistered).toBe(true);

      // Canonical NGL-only extensions: structure=mmtf, surface=ply, density=ccp4.
      expect(res.structureExt).toContain('mmtf');
      expect(res.surfaceExt).toContain('ply');
      expect(res.densityExt).toContain('ccp4');

      // NGL preview extensions must not collide with the Mol*-preview ext sets.
      for (const ext of ['mmtf', 'cns', 'prmtop']) {
        expect(
          res.molstarStructureExt.includes(ext),
          `Mol* preview-structure ext set collides with NGL-structure on '${ext}': ${JSON.stringify(res.molstarStructureExt)}`,
        ).toBe(false);
      }
      expect(
        res.molstarDensityExt.includes('ccp4'),
        `Mol* preview-density ext set collides with NGL-density on 'ccp4': ${JSON.stringify(res.molstarDensityExt)}`,
      ).toBe(false);
    });

    // SCENARIO 4 — Grid context menu Show -> NGL leaf on a Molecule3D cell (force-autostart + retry).
    let scenario4Mounted = false;
    let scenario4Result: any = null;

    await softStep('Scenario 4 setup — Stage Molecule3D table + force BSV autostart for deterministic context-menu wiring', async () => {
      const res = await page.evaluate(async (pdbPath) => {
        const w: any = window;
        const pollUntil = async (pred: () => boolean, timeoutMs = 3000, stepMs = 75) => {
          const t0 = Date.now();
          while (Date.now() - t0 < timeoutMs) {
            try { if (pred()) return true; } catch (_) { /* predicate not ready */ }
            await new Promise((r) => setTimeout(r, stepMs));
          }
          return false;
        };
        grok.shell.closeAll();
        const content = await grok.dapi.files.readAsText(pdbPath);
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['1bdq', '1bdq-clone']),
          DG.Column.fromStrings('structure', [content, content]),
        ]);
        const col = df.col('structure');
        col.semType = DG.SEMTYPE && DG.SEMTYPE.MOLECULE3D ? DG.SEMTYPE.MOLECULE3D : 'Molecule3D';
        try { col.setTag('cell.renderer', 'Molecule3D'); } catch (_) { /* tag set */ }
        try { col.meta.units = 'pdb'; } catch (_) { /* meta API variants */ }
        df.name = 'ngl-extension-context-menu-fixture';
        const tv = grok.shell.addTableView(df);
        await new Promise((resolve) => {
          try {
            const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
            setTimeout(resolve, 5000);
          } catch (_) { resolve(undefined); }
        });
        await pollUntil(() => !!document.querySelector('[name="viewer-Grid"]'));
        // Force-call autostart to wire the context-menu hook.
        let autostartCalled = false;
        try {
          const autoFns = DG.Func.find({package: 'BiostructureViewer', name: 'autostart'});
          if (autoFns && autoFns.length > 0) {
            await autoFns[0].apply({}, {processed: true});
            autostartCalled = true;
          }
        } catch (_) { /* best-effort */ }
        // Poll for the observable effect of autostart + grid render: overlay canvases present.
        await pollUntil(() => {
          const gr = tv && tv.grid ? tv.grid.root : null;
          return !!gr && gr.querySelectorAll('canvas').length >= 3;
        }, 10_000);
        if (!w.$biostructureViewer) w.$biostructureViewer = {};
        w.$biostructureViewer.contextMenuError = null;
        const gridRoot = tv && tv.grid ? tv.grid.root : null;
        const canvases = gridRoot ? Array.from(gridRoot.querySelectorAll('canvas')) : [];
        return {
          rowCount: df.rowCount,
          structureSemType: col.semType,
          hasGridDom: !!document.querySelector('[name="viewer-Grid"]'),
          canvasCount: canvases.length,
          autostartCalled,
        };
      }, sample1bdq);

      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 30_000});

      expect(res.hasGridDom).toBe(true);
      expect(res.rowCount).toBe(2);
      expect(res.structureSemType).toBe('Molecule3D');
      expect(res.canvasCount).toBeGreaterThanOrEqual(3);
      expect(res.autostartCalled).toBe(true);
      scenario4Mounted = true;
    });

    await softStep('Scenario 4 step 2-3 — Right-click populated Molecule3D cell; assert Show -> NGL leaf is injected', async () => {
      if (!scenario4Mounted) return;
      // Clear page errors before the load-bearing dispatch.
      pageErrors.length = 0;

      scenario4Result = await page.evaluate(async () => {
        const w: any = window;
        const tv = w.grok.shell.tv;
        if (!tv || !tv.grid) return {err: 'no grid'};
        const gridRoot = tv.grid.root;
        const gridRect = gridRoot.getBoundingClientRect();
        const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
        const overlay = canvases[2];
        if (!overlay) return {err: 'no overlay canvas'};

        const cell = tv.grid.cell('structure', 0);
        const sb = cell.bounds;
        const cx = gridRect.left + sb.x + Math.min(20, sb.width / 4);
        const cy = gridRect.top + sb.y + sb.height / 2;

        // Retry up to 3 times (cold-start race).
        let menuLabels: string[] = [];
        let hasShow = false, hasNgl = false, hasBio = false;
        let attemptCount = 0;
        for (let attempt = 0; attempt < 3; attempt++) {
          attemptCount = attempt + 1;
          document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
          document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
          await new Promise((r) => setTimeout(r, 400));
          const evt = new PointerEvent('contextmenu', {
            cancelable: true, bubbles: true, view: window, button: 2,
            clientX: cx, clientY: cy,
          });
          overlay.dispatchEvent(evt);
          // Poll for the popup to appear rather than waiting a fixed delay.
          const tPopup = Date.now();
          while (Date.now() - tPopup < 2500 &&
            document.querySelectorAll('.d4-menu-popup .d4-menu-item-label').length === 0)
            await new Promise((r) => setTimeout(r, 75));

          menuLabels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
            .map((el) => (el.textContent || '').trim());
          hasShow = menuLabels.includes('Show');
          hasNgl = menuLabels.includes('NGL');
          hasBio = menuLabels.includes('Biostructure');
          if (hasShow && hasNgl && hasBio) break;
        }

        // If a deferred submenu is used, hover the 'Show' group label to expand it.
        if (hasShow && (!hasNgl || !hasBio)) {
          const showLabelEl = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
            .find((el) => (el.textContent || '').trim() === 'Show');
          const showGroupEl = showLabelEl ? showLabelEl.closest('.d4-menu-item') : null;
          if (showGroupEl) {
            const r = (showGroupEl as HTMLElement).getBoundingClientRect();
            showGroupEl.dispatchEvent(new MouseEvent('mouseover', {
              bubbles: true, clientX: r.left + 10, clientY: r.top + 10,
            }));
            showGroupEl.dispatchEvent(new MouseEvent('mouseenter', {
              bubbles: true, clientX: r.left + 10, clientY: r.top + 10,
            }));
            // Poll for the submenu leaves to expand after hover.
            const tHover = Date.now();
            const hasLeaves = () => {
              const l = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
                .map((el) => (el.textContent || '').trim());
              return l.includes('NGL') && l.includes('Biostructure');
            };
            while (Date.now() - tHover < 2000 && !hasLeaves())
              await new Promise((rr) => setTimeout(rr, 75));
            const afterHoverLabels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
              .map((el) => (el.textContent || '').trim());
            hasNgl = hasNgl || afterHoverLabels.includes('NGL');
            hasBio = hasBio || afterHoverLabels.includes('Biostructure');
          }
        }

        const captured = w.$biostructureViewer && w.$biostructureViewer.contextMenuError;
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });

        return {
          contextMenuErrorIsNull: captured === null || captured === undefined,
          contextMenuErrorMessage: captured ? String(captured.message || captured) : null,
          hasShow, hasNgl, hasBio,
          menuItemCount: menuLabels.length,
          menuItemsSample: menuLabels.filter((t) => t.length > 0 && t.length < 50).slice(0, 30),
          attemptCount,
        };
      });

      // The 'NGL' leaf must be present under the 'Show' group on a populated Molecule3D cell.
      expect(
        scenario4Result.hasShow,
        `Show group missing from grid context menu. Menu items: ${JSON.stringify(scenario4Result.menuItemsSample)}.`,
      ).toBe(true);
      expect(
        scenario4Result.hasNgl,
        `Show -> NGL leaf missing from grid context menu — biostructure.grid-context-menu.show-ngl-viewer regression. ` +
        `Menu items: ${JSON.stringify(scenario4Result.menuItemsSample)}. ` +
        `Attempts: ${scenario4Result.attemptCount}.`,
      ).toBe(true);

      // Coexistence: Show -> Biostructure must appear alongside Show -> NGL.
      expect(
        scenario4Result.hasBio,
        `Coexistence invariant violated: Show -> Biostructure leaf missing alongside Show -> NGL. ` +
        `Menu items: ${JSON.stringify(scenario4Result.menuItemsSample)}.`,
      ).toBe(true);

      // The detectors.js hook must not throw on a populated cell (cross-check with GROK-14552).
      expect(
        scenario4Result.contextMenuErrorIsNull,
        `BSV context-menu hook threw on populated Molecule3D cell: ${JSON.stringify(scenario4Result.contextMenuErrorMessage)}.`,
      ).toBe(true);

      const errSig = pageErrors.filter((m) =>
        /TypeError|ReferenceError|Cannot read properties/i.test(m),
      );
      expect(
        errSig,
        `JS console error during context-menu engagement: ${JSON.stringify(errSig)}.`,
      ).toEqual([]);
    });

    // SCENARIO 5 — PDB id context-panel widget (NGL-wrapped): pane mount + direct pdbIdNglPanelWidget apply.
    let scenario5Mounted = false;

    await softStep('Scenario 5 step 1-3 — Stage PDB_ID table; set current cell; right pane PDB-id-viewer surfaces', async () => {
      pageErrors.length = 0;
      const res = await page.evaluate(async () => {
        const pollUntil = async (pred: () => boolean, timeoutMs = 3000, stepMs = 75) => {
          const t0 = Date.now();
          while (Date.now() - t0 < timeoutMs) {
            try { if (pred()) return true; } catch (_) { /* predicate not ready */ }
            await new Promise((r) => setTimeout(r, stepMs));
          }
          return false;
        };
        grok.shell.closeAll();
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('pdb_id', ['1QBS', '1BNA', '1CRN']),
        ]);
        df.col('pdb_id').semType = 'PDB_ID';
        df.name = 'ngl-extension-pdb-id-fixture';
        const tv = grok.shell.addTableView(df);
        await pollUntil(() => !!document.querySelector('[name="viewer-Grid"]'));
        df.currentRowIdx = 0;
        try { df.currentCell = df.cell(0, 'pdb_id'); } catch (_) { /* setter variants */ }
        // Force the context panel visible and set the current object explicitly to the PDB_ID
        // cell: the baseline sets simpleMode=true (hides the context panel), and the semantic
        // panes only surface into a shown panel driven by grok.shell.o.
        grok.shell.windows.showProperties = true;
        try { grok.shell.o = DG.SemanticValue.fromTableCell(df.cell(0, 'pdb_id')); } catch (_) { /* fall back to currentCell */ }
        await pollUntil(() => df.currentRowIdx === 0);
        return {
          rowCount: df.rowCount,
          pdbIdSemType: df.col('pdb_id').semType,
          currentRowIdx: df.currentRowIdx,
          paneNames: Array.from(document.querySelectorAll('[name^="pane-"]'))
            .map((el) => el.getAttribute('name')),
        };
      });
      expect(res.rowCount).toBe(3);
      expect(res.pdbIdSemType).toBe('PDB_ID');
      expect(res.currentRowIdx).toBe(0);
      await page.locator('[name="pane-PDB-id-viewer"]').waitFor({timeout: 30_000});
      expect(res.paneNames).toContain('pane-PDB-id-viewer');
      scenario5Mounted = true;
    });

    await softStep('Scenario 5 step 4 — pdbIdNglPanelWidget(pdbId) returns a renderable DG.Widget (DG.Func.find apply)', async () => {
      if (!scenario5Mounted) return;
      const res = await page.evaluate(async () => {
        const fns = DG.Func.find({name: 'pdbIdNglPanelWidget', package: 'BiostructureViewer'});
        const fn = fns && fns[0];
        if (!fn) return {ok: false, reason: 'pdbIdNglPanelWidget not registered'};
        let widget: any = null;
        let err: string | null = null;
        try { widget = await fn.apply({pdbId: '1QBS'}); }
        catch (e: any) { err = String(e?.message || e); }
        // Stash the root so step 5 can prove the re-invocation produced a fresh widget.
        const w: any = window;
        w.__nglStep4Root = widget && widget.root ? widget.root : null;
        return {
          ok: !!widget,
          err,
          widgetRootTagName: widget && widget.root ? widget.root.tagName : null,
          widgetRootNotNull: !!(widget && widget.root),
          // NGL Stage mounts its viewport into the widget host synchronously (RCSB structure
          // load is async/outbound and not asserted here).
          widgetRootChildCount: widget && widget.root ? widget.root.childElementCount : 0,
          widgetRootClass: widget && widget.root ? widget.root.className : null,
          panelRole: fn.options?.role,
          panelName: fn.options?.name,
          inputSemType: fn.inputs?.[0]?.semType,
        };
      });
      expect(
        res.ok,
        `pdbIdNglPanelWidget apply failed: ${JSON.stringify(res)}`,
      ).toBe(true);
      expect(res.widgetRootNotNull).toBe(true);
      expect(
        res.widgetRootChildCount,
        `pdbIdNglPanelWidget produced an empty root (no NGL viewport mounted): ${JSON.stringify(res)}`,
      ).toBeGreaterThan(0);
      expect(res.widgetRootClass).toContain('d4-ngl-viewer');
      expect(res.panelRole).toBe('panel');
      expect(res.inputSemType).toBe('PDB_ID');
    });

    await softStep('Scenario 5 step 5 — Switch current cell to a different PDB_ID; widget re-renders cleanly', async () => {
      if (!scenario5Mounted) return;
      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const df = tv?.dataFrame;
        if (!df) return {ok: false, reason: 'no df'};
        df.currentRowIdx = 1;
        try { df.currentCell = df.cell(1, 'pdb_id'); } catch (_) { /* setter variants */ }
        const t0 = Date.now();
        while (Date.now() - t0 < 3000 && df.currentRowIdx !== 1)
          await new Promise((r) => setTimeout(r, 75));
        // Re-invoke the widget func with the new PDB_ID.
        const fns = DG.Func.find({name: 'pdbIdNglPanelWidget', package: 'BiostructureViewer'});
        const fn = fns && fns[0];
        let widget: any = null;
        let err: string | null = null;
        try { widget = await fn.apply({pdbId: '1BNA'}); }
        catch (e: any) { err = String(e?.message || e); }
        const w: any = window;
        return {
          ok: !!widget,
          err,
          currentRowIdx: df.currentRowIdx,
          widgetRootNotNull: !!(widget && widget.root),
          widgetRootChildCount: widget && widget.root ? widget.root.childElementCount : 0,
          // A fresh widget instance for the new id is a distinct root element from step 4's.
          rootDiffersFromPrev: !!(widget && widget.root && widget.root !== w.__nglStep4Root),
          paneStillPresent: !!document.querySelector('[name="pane-PDB-id-viewer"]'),
        };
      });
      expect(res.ok).toBe(true);
      expect(res.currentRowIdx).toBe(1);
      expect(res.widgetRootNotNull).toBe(true);
      // The re-invocation for '1BNA' must produce a fresh, non-empty widget (not the stale '1QBS' one).
      // The id-specific 3D structure loads async from RCSB and is not asserted here.
      expect(res.rootDiffersFromPrev).toBe(true);
      expect(res.widgetRootChildCount).toBeGreaterThan(0);
      const errSig = pageErrors.filter((m) =>
        /TypeError|ReferenceError|Cannot read properties/i.test(m),
      );
      expect(
        errSig,
        `JS console error during PDB id viewer panel widget mount/re-mount: ${JSON.stringify(errSig)}.`,
      ).toEqual([]);
    });
  } finally {
    // Cleanup — no server-side state was created by this spec.
    try {
      await page.evaluate(() => {
        const w: any = window;
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
        document.querySelectorAll('.d4-dialog').forEach((d) => {
          const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
          if (cancel) cancel.click();
        });
        if (w.$biostructureViewer) w.$biostructureViewer.contextMenuError = null;
        try { w.grok.shell.closeAll(); } catch (_) { /* best-effort */ }
      });
    } catch (e) { /* best-effort */ }
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
