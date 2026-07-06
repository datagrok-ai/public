/* ---
sub_features_covered: [bio.manage.libraries-app, bio.manage.libraries-app.tree-browser, bio.manage.match-with-library, bio.manage.monomers-view, bio.manage.standardize-library]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
test('Bio Manage Monomer Libraries CRUD (app + tree browser + Monomers view + Match dispatch + standardize)', async ({page}) => {
  test.setTimeout(120_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
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
      for (let i = 0; i < 60; i++) {
        const c = document.querySelector('[name="viewer-Grid"] canvas') as HTMLElement | null;
        if (c && c.clientWidth > 0) break;
        await new Promise((r) => setTimeout(r, 200));
      }
    }
  }, 'System:AppData/Bio/tests/filter_HELM.csv');
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getMonomerLibHelper', 'Bio:getSeqHelper', 'Bio:getBioLib'];
    for (let i = 0; i < 15; i++) {
      for (const fn of probes) {
        try { await (grok as any).functions.call(fn, {}); return; } catch {  }
      }
      await new Promise((r) => setTimeout(r, 300));
    }
  });
  await page.evaluate(() => {
    const w: any = window as any;
    if (w.__balloonInstalled) return;
    w.__balloonInstalled = true;
    w.__balloonErrors = 0;
    w.__balloonWarnings = 0;
    const origErr = grok.shell.error?.bind(grok.shell);
    const origWarn = grok.shell.warning?.bind(grok.shell);
    if (origErr) {
      (grok.shell as any).error = (msg: any) => { w.__balloonErrors += 1; return origErr(msg); };
    }
    if (origWarn) {
      (grok.shell as any).warning = (msg: any) => { w.__balloonWarnings += 1; return origWarn(msg); };
    }
  });
  try {
    const balloonBefore1 = await page.evaluate(() => ({
      err: (window as any).__balloonErrors || 0,
      warn: (window as any).__balloonWarnings || 0,
    }));
    await softStep('S1.1-1.2: Bio:Manage Monomer Libraries app invocation opens the Manage Monomer Libraries view', async () => {
      // Direct app-entry invocation (bio.manage.libraries-app), disambiguated by the app tag
      // from the same-named func that opens a dialog.
      const invokeErr = await page.evaluate(async () => {
        try {
          // DG.Func.find matches the runtime Func.name (the JS export identifier), not the
          // decorator's friendly `name:` — the app registers as `manageMonomerLibrariesView`
          // (friendlyName 'Manage Monomer Libraries'). Querying the friendly name returns [].
          const fns = (window as any).DG.Func.find({package: 'Bio', name: 'manageMonomerLibrariesView', tags: ['app']});
          if (!fns || fns.length === 0) return 'Bio app "Manage Monomer Libraries" not registered';
          // The app entry returns its view with addView=false; when invoking the function
          // directly (not through the platform app-launcher) the returned view must be
          // attached to the shell ourselves so it becomes grok.shell.v.
          const v = await fns[0].apply();
          if (v && (window as any).grok.shell.v !== v) (window as any).grok.shell.addView(v);
          return null;
        } catch (e) { return String(e).slice(0, 250); }
      });
      expect(invokeErr, `Bio:Manage Monomer Libraries app invocation failed: ${invokeErr}`).toBeNull();
      await page.waitForFunction(() => {
        try {
          const n = (window as any).grok?.shell?.v?.name;
          return typeof n === 'string' && n === 'Manage Monomer Libraries';
        } catch (_) { return false; }
      }, null, {timeout: 45_000});
      const result = await page.evaluate(() => {
        const v = grok.shell.v;
        return {
          viewName: v?.name || null,
          viewType: (v as any)?.type || null,
          viewRootPresent: !!v?.root,
        };
      });
      expect(result.viewName).toBe('Manage Monomer Libraries');
      expect(result.viewType).toBe('view');
      expect(result.viewRootPresent).toBe(true);
    });
    await softStep('S1.3-1.4: Monomer Manager Tree Browser populates ≥1 library node + first node click is no-throw', async () => {
      const result = await page.evaluate(async () => {
        const treeRoot: any = (window as any).ui?.tree?.();
        if (!treeRoot)
          return {invokeErr: 'ui.tree() factory not available', nodeCount: 0, nodeNames: [] as string[], firstNodeClickErr: null};
        let invokeErr: string | null = null;
        try {
          // Call by the runtime nqName (JS export), not the space-containing friendly name —
          // grok.functions.call cannot parse 'Bio:Monomer Manager Tree Browser'.
          await (grok as any).functions.call('Bio:manageMonomerLibrariesViewTreeBrowser', {treeNode: treeRoot});
        } catch (e) {
          invokeErr = String(e).slice(0, 250);
        }
        const readItems = () => (treeRoot.items || treeRoot.children || []) as any[];
        for (let i = 0; i < 75; i++) {
          if (readItems().length > 0) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        const items = readItems();
        const nodeCount = items.length;
        const nodeNames = items.map((n: any) => {
          try { return String(n.text || n.value || n.name || ''); } catch (_) { return ''; }
        }).filter((s: string) => s.length > 0);
        let firstNodeClickErr: string | null = null;
        if (items.length > 0) {
          const node = items[0];
          try {
            if (node?.onSelected?.next) node.onSelected.next(node);
            else if (typeof node?.select === 'function') node.select();
            else if (typeof node?.click === 'function') node.click();
          } catch (e) { firstNodeClickErr = String(e).slice(0, 200); }
          // Bounded wait for the per-library monomer manager surface to mount (bounded posture — not asserted).
          for (let i = 0; i < 50; i++) {
            const root: any = (grok.shell.v as any)?.root;
            if (root?.querySelector?.('.d4-grid, .grok-grid, canvas')) break;
            await new Promise((r) => setTimeout(r, 200));
          }
        }
        return {invokeErr, nodeCount, nodeNames, firstNodeClickErr};
      });
      expect(result.invokeErr,
        `Bio:Monomer Manager Tree Browser did not resolve: ${result.invokeErr}`).toBeNull();
      expect(result.nodeCount,
        `expected ≥1 tree-browser library node; nodes=[${result.nodeNames.join(', ')}]`).toBeGreaterThanOrEqual(1);
      expect(result.firstNodeClickErr,
        `first tree-node click threw: ${result.firstNodeClickErr}`).toBeNull();
    });
    await softStep('S1.5: no error balloon raised during Scenario 1', async () => {
      const balloonAfter = await page.evaluate(() => ({
        err: (window as any).__balloonErrors || 0,
        warn: (window as any).__balloonWarnings || 0,
      }));
      expect(balloonAfter.err - balloonBefore1.err,
        `error balloon count increased by ${balloonAfter.err - balloonBefore1.err} during Scenario 1`).toBe(0);
    });
    // Scenario 2 — Bio | Manage | Monomers view open
    const balloonBefore2 = await page.evaluate(() => ({
      err: (window as any).__balloonErrors || 0,
      warn: (window as any).__balloonWarnings || 0,
    }));
    // Bring the HELM TableView forward — the Bio top-menu vanishes while the app view holds focus.
    await page.evaluate(async () => {
      const tvs: any[] = Array.from((grok.shell as any).tableViews || []);
      const helm: any = tvs.find((tv: any) => {
        const df = tv?.dataFrame;
        if (!df) return false;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.some((c: any) => c.semType === 'Macromolecule');
      });
      if (helm) {
        try { (grok.shell as any).v = helm; } catch (_) { /* read-only on some builds */ }
        for (let i = 0; i < 25; i++) {
          if ((grok.shell as any).v === helm) break;
          await new Promise((r) => setTimeout(r, 100));
        }
      }
    });
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    await softStep('S2.1-2.3: Bio | Manage | Monomers top-menu opens a view with the expected shape', async () => {
      await page.evaluate(async () => {
        const wait = async (sel: string) => {
          for (let i = 0; i < 50; i++) {
            const e = document.querySelector(sel) as HTMLElement | null;
            if (e) return e;
            await new Promise((r) => setTimeout(r, 100));
          }
          return null;
        };
        const bioMenu = document.querySelector('[name="div-Bio"]') as HTMLElement | null;
        if (!bioMenu) throw new Error('Bio top-menu [name="div-Bio"] not present after HELM TableView re-focus');
        bioMenu.click();
        const manage = await wait('[name="div-Bio---Manage"]');
        if (!manage) throw new Error('[name="div-Bio---Manage"] submenu not present');
        manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        const leaf = await wait('[name="div-Bio---Manage---Monomers"]');
        if (leaf) leaf.click();
      });
      await page.waitForFunction(() => {
        try {
          const n = (window as any).grok?.shell?.v?.name;
          return typeof n === 'string' && n.toLowerCase().includes('monomer');
        } catch (_) { return false; }
      }, null, {timeout: 30_000});
      await page.waitForFunction(() => {
        try {
          const root: any = (window as any).grok?.shell?.v?.root;
          return !!root && root.querySelectorAll('[name^="viewer-"], .d4-grid, .grok-grid, .d4-tree-view-root').length > 0;
        } catch (_) { return false; }
      }, null, {timeout: 30_000});
      const result = await page.evaluate(() => {
        const v = grok.shell.v;
        const root: any = v?.root;
        const rootIsElement = !!root && (root.nodeType === 1 || typeof root.querySelector === 'function');
        let firstChildTag: string | null = null;
        let hasChildElement = false;
        if (rootIsElement) {
          const candidates = root.querySelectorAll('[name^="viewer-"], .d4-grid, .grok-grid, .d4-tree-view-root');
          hasChildElement = candidates.length > 0;
          if (candidates.length > 0) firstChildTag = candidates[0].tagName.toLowerCase();
        }
        return {
          viewName: v?.name || null,
          viewType: v?.type || null,
          rootIsElement,
          hasChildElement,
          firstChildTag,
        };
      });
      expect((result.viewName || '').toLowerCase()).toContain('monomer');
      expect(result.viewType).toBe('view');
      expect(result.rootIsElement).toBe(true);
      expect(result.hasChildElement,
        `expected the Manage Monomers monomer-list surface under the view root; firstChildTag=${result.firstChildTag}`).toBe(true);
    });
    // Scenario 2 Expected: no error balloon raised.
    await softStep('S2.4: no error balloon raised during Scenario 2', async () => {
      const balloonAfter = await page.evaluate(() => ({
        err: (window as any).__balloonErrors || 0,
        warn: (window as any).__balloonWarnings || 0,
      }));
      expect(balloonAfter.err - balloonBefore2.err,
        `error balloon count increased by ${balloonAfter.err - balloonBefore2.err} during Scenario 2`).toBe(0);
    });
    // Scenario 3 — Match with Monomer Library dispatch + standardiseMonomerLibrary normalization
    const balloonBefore3 = await page.evaluate(() => ({
      err: (window as any).__balloonErrors || 0,
      warn: (window as any).__balloonWarnings || 0,
    }));
    await page.evaluate(async () => {
      const tvs: any[] = Array.from((grok.shell as any).tableViews || []);
      const helm: any = tvs.find((tv: any) => {
        const df = tv?.dataFrame;
        if (!df) return false;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.some((c: any) => c.semType === 'Macromolecule');
      });
      if (helm) {
        try { (grok.shell as any).v = helm; } catch (_) { /* read-only on some builds */ }
        for (let i = 0; i < 25; i++) {
          if ((grok.shell as any).v === helm) break;
          await new Promise((r) => setTimeout(r, 100));
        }
      }
    });
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 15_000});
    await softStep('S3.2-3.5: Match-with-Monomer-Library dialog opens with three host inputs + Polymer-Type carries PEPTIDE/RNA/CHEM', async () => {
      await page.evaluate(async () => {
        const wait = async (sel: string) => {
          for (let i = 0; i < 50; i++) {
            const e = document.querySelector(sel) as HTMLElement | null;
            if (e) return e;
            await new Promise((r) => setTimeout(r, 100));
          }
          return null;
        };
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        const manage = await wait('[name="div-Bio---Manage"]');
        if (!manage) throw new Error('[name="div-Bio---Manage"] submenu not present');
        manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        const leaf = await wait('[name="div-Bio---Manage---Match-with-Monomer-Library..."]');
        if (leaf) leaf.click();
      });
      await page.locator('[name="dialog-matchWithMonomerLibrary"]').waitFor({state: 'visible', timeout: 30_000});
      const result = await page.evaluate(() => {
        const dialog = document.querySelector('[name="dialog-matchWithMonomerLibrary"]');
        const hostTable = dialog?.querySelector('[name="input-host-Table"]') ?? null;
        const hostMolecules = dialog?.querySelector('[name="input-host-Molecules"]') ?? null;
        const hostPolymer = dialog?.querySelector('[name="input-host-Polymer-Type"]') ?? null;
        let polymerOptions: string[] = [];
        if (hostPolymer) {
          const sel = hostPolymer.querySelector('select') as HTMLSelectElement | null;
          if (sel) {
            polymerOptions = Array.from(sel.options).map((o) => o.value || o.textContent || '').filter((s) => s.length > 0);
          } else {
            const opts = hostPolymer.querySelectorAll('option, [role="option"], .d4-combo-popup-item');
            polymerOptions = Array.from(opts).map((o: any) => (o.value || o.textContent || '').trim()).filter((s: string) => s.length > 0);
          }
        }
        return {
          dialogPresent: !!dialog,
          hostTablePresent: !!hostTable,
          hostMoleculesPresent: !!hostMolecules,
          hostPolymerPresent: !!hostPolymer,
          polymerOptions,
        };
      });
      expect(result.dialogPresent).toBe(true);
      expect(result.hostTablePresent).toBe(true);
      expect(result.hostMoleculesPresent).toBe(true);
      expect(result.hostPolymerPresent).toBe(true);
      expect(result.polymerOptions.length,
        `Polymer-Type select produced no options; observed: [${result.polymerOptions.join(', ')}]`).toBeGreaterThanOrEqual(3);
      const upper = result.polymerOptions.map((s) => s.toUpperCase());
      for (const opt of ['PEPTIDE', 'RNA', 'CHEM']) {
        expect(upper,
          `expected Polymer-Type to include '${opt}'; observed: [${result.polymerOptions.join(', ')}]`).toContain(opt);
      }
    });
    await softStep('S3.6: standardiseMonomerLibrary resolves; normalized payload parses as a non-null object', async () => {
      const result = await page.evaluate(async () => {
        const stamp = Date.now();
        const sym = `XYZ_STD_${stamp}`;
        const lib: any[] = [{
          symbol: sym,
          name: sym,
          molfile: '\n     RDKit          2D\n\n  7  6  0  0  0  0  0  0  0  0999 V2000\n    1.6702    1.3929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3712    0.6429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9279    1.3929    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.2269    0.6429    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    0.3712   -0.8571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9279   -1.6071    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6702   -1.6071    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1\n  2  3  1  0\n  3  4  1  0\n  2  5  1  0\n  5  6  2  0\n  5  7  1  0\nM  RGP  2   4   1   7   2\nM  END\n',
          smiles: 'C[C@H](N[*:1])C(=O)[*:2]',
          polymerType: 'PEPTIDE',
          monomerType: 'Backbone',
          naturalAnalog: 'X',
          id: 0,
          rgroups: [
            {alternateId: 'R1-H', capGroupName: 'H', capGroupSMILES: '[*:1][H]', label: 'R1'},
            {alternateId: 'R2-OH', capGroupName: 'OH', capGroupSMILES: 'O[*:2]', label: 'R2'},
          ],
        }];
        const payload = JSON.stringify(lib);
        let stdErr: string | null = null;
        let returned: string | null = null;
        try {
          returned = await (grok as any).functions.call('Bio:standardiseMonomerLibrary', {library: payload});
        } catch (e) {
          stdErr = String(e).slice(0, 250);
        }
        let parsedIsArray = false;
        let parsedLength = -1;
        let firstPolymerType: string | null = null;
        let firstSymbol: string | null = null;
        let parseErr: string | null = null;
        if (typeof returned === 'string' && returned.length > 0) {
          try {
            const reparsed: any = JSON.parse(returned);
            parsedIsArray = Array.isArray(reparsed);
            if (parsedIsArray) {
              parsedLength = reparsed.length;
              if (reparsed[0]) {
                firstPolymerType = reparsed[0].polymerType ?? null;
                firstSymbol = reparsed[0].symbol ?? null;
              }
            }
          } catch (e) {
            parseErr = String(e).slice(0, 200);
          }
        }
        return {
          stdErr,
          returnedType: typeof returned,
          returnedNonEmpty: typeof returned === 'string' && returned.length > 0,
          parsedIsArray,
          parsedLength,
          firstPolymerType,
          firstSymbol,
          parseErr,
        };
      });
      expect(result.stdErr,
        `Bio:standardiseMonomerLibrary threw: ${result.stdErr}`).toBeNull();
      expect(result.returnedType).toBe('string');
      expect(result.returnedNonEmpty).toBe(true);
      expect(result.parseErr,
        `normalized library did not parse: ${result.parseErr}`).toBeNull();
      expect(result.parsedIsArray,
        'normalized library should JSON-parse to an array of monomers').toBe(true);
      expect(result.parsedLength,
        'normalized library should contain the standardized monomer').toBeGreaterThanOrEqual(1);
      expect(result.firstPolymerType,
        `standardized monomer should round-trip polymerType PEPTIDE; got ${result.firstPolymerType}`).toBe('PEPTIDE');
      expect((result.firstSymbol || '').length,
        'standardized monomer should carry a non-empty symbol').toBeGreaterThan(0);
    });
    await softStep('S3.7: no error balloon raised during Scenario 3', async () => {
      const balloonAfter = await page.evaluate(() => ({
        err: (window as any).__balloonErrors || 0,
        warn: (window as any).__balloonWarnings || 0,
      }));
      expect(balloonAfter.err - balloonBefore3.err,
        `error balloon count increased by ${balloonAfter.err - balloonBefore3.err} during Scenario 3`).toBe(0);
    });
  } finally {
    // Cleanup — close any open dialogs / manage views (best-effort).
    await page.evaluate(async () => {
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      for (const d of dialogs)
        d.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      try {
        const views = Array.from(grok.shell.views || []);
        for (const v of views) {
          const n: string = (v as any)?.name || '';
          const lower = n.toLowerCase();
          if ((lower.includes('manage') && lower.includes('monomer')) ||
              lower === 'manage monomers') {
            try { (v as any).close?.(); } catch (_) { /* best effort */ }
          }
        }
      } catch (_) { /* best effort */ }
    }).catch(() => {});
  }
  finishSpec();
});
