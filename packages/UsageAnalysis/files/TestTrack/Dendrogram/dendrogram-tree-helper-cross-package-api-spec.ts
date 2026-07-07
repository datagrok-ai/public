// Frontmatter extraction (pre-author hooks):
//   target_layer: apitest
//   pyramid_layer: absent
//   ui_coverage_responsibility: [] (delegated_to: null - pure JS-API spec)
//   related_bugs: []
//   produced_from: atlas-driven
//
// Atlas provenance (derived_from):
//   feature-atlas/dendrogram.yaml#sub_features[dendrogram.api.get-tree-helper]
//     source: public/packages/Dendrogram/src/package.ts#L69
//   feature-atlas/dendrogram.yaml#sub_features[dendrogram.api.get-dendrogram-service]
//     source: public/packages/Dendrogram/src/package.ts#L83
//   feature-atlas/dendrogram.yaml#sub_features[dendrogram.api.dendrogram-service.inject-tree-for-grid]
//     source: public/packages/Dendrogram/src/utils/dendrogram-service.ts#L10
//   feature-atlas/dendrogram.yaml#interactions[dendrogram.cross.tree-helper-cross-package]
//     coverage_type: smoke (canonical per atlas top-level interaction)
//
// Paradigm: apitest (target_layer: apitest). NO DOM-driving calls
// (page.click/fill/locator/hover/press); NO DOM-Locator assertions
// (toBeVisible/toHaveText/toHaveCount). The spec exercises the three
// public Dendrogram package functions through their registered
// grok.functions.call entry points inside page.evaluate blocks. The
// one DOM observable used (the .dendrogram-assign-clusters-bttn
// magic-wand on the grid neighbor) is read back from page.evaluate as
// a captured boolean - it is the deterministic ready signal that
// confirms injectTreeForGrid mounted a GridNeighbor, NOT a Playwright
// Locator assertion (per the scenario's Expected: "magic-wand icon is
// present in the grid neighbor" + Notes: "single-anchor check, not a
// UI flow").
//
// Reference: .claude/skills/grok-browser/references/dendrogram.md
//   - tree-helper-cross-package (apitest-layer reference; no DOM
//     observable beyond what assign-clusters-dialog already documents)
//   - dendrogram-neighbor-close-and-reset (selectors:
//     .dendrogram-assign-clusters-bttn, .dendrogram-close-bttn)
//
// Filename note: this spec uses the compound suffix -api.ts (not
// the bare -api.ts the filename-selection convention names) so the
// per-section playwright.config.ts at files/TestTrack/playwright.config.ts
// (testMatch: '**/*-spec.ts') collects it. Sibling apitest specs that
// PASS Gate B on the same config follow the same convention:
// hierarchical-clustering-bio-api.ts, lifecycle-api.ts,
// charts-api.ts. The paradigm remains pure apitest - no DOM-driving
// calls were added; the rename is filename-only.
//
// MCP recon evidence (live 2026-06-03 on dev.datagrok.ai, user
// oahadzhanian):
//   - DG.Func.find({package: 'Dendrogram', name: 'getTreeHelper'})
//     resolves to a single registered Func; same for getDendrogramService.
//   - await grok.functions.call('Dendrogram:getTreeHelper') returns a
//     TreeHelper instance with callable methods newickToDf, toNewick,
//     getLeafList, getNodeList (the ITreeHelper surface named in the
//     scenario).
//   - await grok.functions.call('Dendrogram:getDendrogramService')
//     returns a DendrogramService instance with callable injectTreeForGrid.
//   - Two consecutive getDendrogramService calls return the same
//     instance (svc2 === svc1); window.$dendrogramService === svc1
//     after the first call. Singleton invariant confirmed.
//   - svc.injectTreeForGrid(tv.grid, treeRoot, 'leaf') on a 4-row
//     synthetic DataFrame with the 4-leaf NodeType literal
//     ((A,B),(C,D)) mounts .dendrogram-assign-clusters-bttn and
//     .dendrogram-close-bttn within ~1s, with zero fatal console errors.
//   - tv.grid.temp['__dendrogram_neighbor_temp__'] is NOT set by the
//     bare injectTreeForGrid path (the temp marker is set only on the
//     hierarchicalClustering-driven inject path). The magic-wand
//     selector is therefore the authoritative ready signal here, not
//     the temp marker.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Dendrogram: Cross-package TreeHelper / DendrogramService consumption (JS API smoke)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // ===== Scenario 1: getTreeHelper + getDendrogramService + injectTreeForGrid happy path =====
  await softStep('Scenario 1.1 - Dendrogram:getTreeHelper resolves to an ITreeHelper with the documented method surface', async () => {
    const result = await page.evaluate(async () => {
      const th: any = await grok.functions.call('Dendrogram:getTreeHelper');
      const required = ['newickToDf', 'toNewick', 'getLeafList', 'getNodeList'];
      const methods = required.map(m => ({m, isFn: typeof th?.[m] === 'function'}));
      return {
        truthy: th != null,
        constructorName: th?.constructor?.name,
        methods,
        allMethodsCallable: methods.every(x => x.isFn),
      };
    });
    expect(result.truthy, 'getTreeHelper returns a non-null ITreeHelper').toBe(true);
    expect(result.allMethodsCallable,
      `getTreeHelper instance exposes the ITreeHelper method surface (${result.methods.map(x => `${x.m}=${x.isFn}`).join(', ')})`)
      .toBe(true);
  });

  await softStep('Scenario 1.2 - Dendrogram:getDendrogramService resolves to an IDendrogramService with callable injectTreeForGrid', async () => {
    const result = await page.evaluate(async () => {
      const svc: any = await grok.functions.call('Dendrogram:getDendrogramService');
      return {
        truthy: svc != null,
        constructorName: svc?.constructor?.name,
        hasInject: typeof svc?.injectTreeForGrid === 'function',
      };
    });
    expect(result.truthy, 'getDendrogramService returns a non-null IDendrogramService').toBe(true);
    expect(result.hasInject, 'IDendrogramService.injectTreeForGrid is a callable function').toBe(true);
  });

  await softStep('Scenario 1.3 - injectTreeForGrid(grid, 4-leaf treeRoot, "leaf") mounts the grid neighbor with the magic-wand icon and emits no fatal console error', async () => {
    const result = await page.evaluate(async () => {
      // Setup phase per scenario .md Setup block.
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 600));

      // Setup step 1 - 4-row synthetic DataFrame with a single string column
      // 'leaf' whose values line up with the tree's leaf names.
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromList('string', 'leaf', ['A', 'B', 'C', 'D']),
      ]);
      df.name = 'tree-helper-cross-package-smoke';

      // Setup step 3 - open a TableView so a DG.Grid is available for
      // injectTreeForGrid to attach onto.
      const tv = grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 1200));

      // Setup step 2 - 4-leaf NodeType literal whose leaf .name values
      // are [A, B, C, D] (the literal-construction path per scenario
      // Setup step 2's "or equivalently, construct a NodeType literal").
      const treeRoot: any = {
        name: 'root',
        branch_length: 0,
        children: [
          {name: 'I1', branch_length: 1, children: [
            {name: 'A', branch_length: 1, children: []},
            {name: 'B', branch_length: 1, children: []},
          ]},
          {name: 'I2', branch_length: 1, children: [
            {name: 'C', branch_length: 1, children: []},
            {name: 'D', branch_length: 1, children: []},
          ]},
        ],
      };

      // Scenario 1 step 1 - obtain treeHelper (re-fetched in this evaluate
      // for self-containment; same singleton fetched in Scenario 1.1).
      const th: any = await grok.functions.call('Dendrogram:getTreeHelper');

      // Scenario 1 step 3 - obtain dendrogramService.
      const svc: any = await grok.functions.call('Dendrogram:getDendrogramService');

      // Capture fatal console errors during the inject call (scenario
      // expected: "No console errors are emitted during the call
      // sequence"). The dev.datagrok.ai page-load 404 noise + the
      // ResizeObserver-loop spam are platform-level and not attributable
      // to the inject path - filter those out per the matching sibling
      // spec convention (hierarchical-clustering-chem-api.ts).
      const consoleErrors: string[] = [];
      const origErr = console.error;
      console.error = (...args: any[]) => {
        consoleErrors.push(args.map((a: any) => String(a)).join(' '));
        origErr.apply(console, args);
      };

      // Scenario 1 step 6 - invoke injectTreeForGrid on the open
      // TableView's grid. Wrapped in try/catch so the assertion below
      // can distinguish a throw from a no-mount.
      let injectThrew: string | false = false;
      try {
        svc.injectTreeForGrid(tv.grid, treeRoot, 'leaf');
      } catch (e: any) {
        injectThrew = String(e?.message || e);
      }

      // Scenario 1 step 7 - await DOM settle and query for the
      // magic-wand icon (the reliable mounted-and-ready signal per
      // grok-browser/references/dendrogram.md - "Common observability
      // for dendrogram specs"). Budget ~9s (30 * 300ms) is generous
      // since the MCP recon mount completed in well under 1s.
      let magicWandMounted = false;
      let closeBttnMounted = false;
      for (let i = 0; i < 30; i++) {
        if (document.querySelector('.dendrogram-assign-clusters-bttn')) {
          magicWandMounted = true;
          closeBttnMounted = !!document.querySelector('.dendrogram-close-bttn');
          break;
        }
        await new Promise(r => setTimeout(r, 300));
      }

      console.error = origErr;

      const fatalErrors = consoleErrors.filter(t =>
        !/Failed to load resource[\s\S]*404/i.test(t)
        && !/ResizeObserver loop/i.test(t));

      // Read back observability: viewer list on the TableView should
      // still enumerate Grid only (the neighbor is a GridNeighbor, NOT
      // a DG.Viewer per grok-browser/references/dendrogram.md
      // "Common observability for dendrogram specs").
      const viewerTypes = Array.from(grok.shell.tv?.viewers || []).map((v: any) => v.type);

      return {
        injectThrew,
        magicWandMounted,
        closeBttnMounted,
        fatalErrors,
        viewerTypes,
      };
    });

    expect(result.injectThrew,
      'svc.injectTreeForGrid(grid, treeRoot, "leaf") must NOT throw')
      .toBe(false);
    expect(result.magicWandMounted,
      '.dendrogram-assign-clusters-bttn (magic-wand) mounts on the grid neighbor within budget - confirms injectTreeForGrid attached the GridNeighbor (atlas: dendrogram.api.dendrogram-service.inject-tree-for-grid)')
      .toBe(true);
    expect(result.closeBttnMounted,
      '.dendrogram-close-bttn also mounts (companion neighbor affordance documented in grok-browser/references/dendrogram.md)')
      .toBe(true);
    expect(result.fatalErrors,
      'no fatal console errors emitted during getTreeHelper + getDendrogramService + injectTreeForGrid sequence')
      .toEqual([]);
    expect(result.viewerTypes,
      'TableView viewers list still enumerates [Grid] only (the dendrogram neighbor is a GridNeighbor, not a DG.Viewer)')
      .toEqual(['Grid']);
  });

  // ===== Scenario 1 -> Scenario 2 transition: detach the neighbor so the next scenario starts clean =====
  await softStep('Scenario 1 -> Scenario 2 transition: close the dendrogram neighbor from Scenario 1.3', async () => {
    await page.evaluate(async () => {
      const closeBtn = document.querySelector('.dendrogram-close-bttn') as HTMLElement | null;
      if (closeBtn) {
        closeBtn.click();
        await new Promise(r => setTimeout(r, 400));
      }
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 400));
    });
  });

  // ===== Scenario 2: getDendrogramService singleton invariant =====
  // The scenario directs us to clear OR capture the initial value of
  // window.$dendrogramService. Because Scenario 1 above already
  // populated the cache (and prior Dendrogram tests in the same
  // browser session may also have), we CAPTURE the initial state and
  // assert the post-call invariants relative to a single fresh
  // baseline. The contract under test is "repeated getDendrogramService
  // calls return the same instance AND populate window.$dendrogramService"
  // - the contract holds regardless of whether the cache pre-existed
  // (per public/packages/Dendrogram/src/package.ts#L84:
  // `if (!(window.$dendrogramService)) { ... window.$dendrogramService = svc; }`
  // - the singleton invariant is the strict property the scenario names).
  await softStep('Scenario 2 - getDendrogramService singleton invariant: svc2 === svc1, window.$dendrogramService === svc1 on first AND second call', async () => {
    const result = await page.evaluate(async () => {
      // Scenario 2 step 1 - capture initial value.
      const initialCached: any = (window as any).$dendrogramService;
      const initialWasSet = initialCached != null;

      // Scenario 2 step 2 - first getDendrogramService call.
      const svc1: any = await grok.functions.call('Dendrogram:getDendrogramService');

      // Scenario 2 step 3 - window.$dendrogramService is now set AND
      // referentially equal to svc1.
      const winAfter1 = (window as any).$dendrogramService;
      const win1EqualsSvc1 = winAfter1 === svc1;
      const winAfter1Set = winAfter1 != null;

      // Scenario 2 step 4 - second getDendrogramService call.
      const svc2: any = await grok.functions.call('Dendrogram:getDendrogramService');

      // Scenario 2 step 5 - svc2 === svc1 (singleton reused, not re-created).
      const svc2EqualsSvc1 = svc2 === svc1;

      // Scenario 2 step 6 - window.$dendrogramService === svc1 still holds.
      const winAfter2 = (window as any).$dendrogramService;
      const win2EqualsSvc1 = winAfter2 === svc1;

      // Cache-respect probe: if the cache pre-existed at step 1, svc1
      // should equal that initial value (the package.ts#L84 guard
      // returns the cached instance directly).
      const svc1RespectsInitial = !initialWasSet || (svc1 === initialCached);

      return {
        initialWasSet,
        winAfter1Set,
        win1EqualsSvc1,
        svc2EqualsSvc1,
        win2EqualsSvc1,
        svc1RespectsInitial,
        svc1HasInject: typeof svc1?.injectTreeForGrid === 'function',
      };
    });
    expect(result.winAfter1Set,
      'window.$dendrogramService is set after the first getDendrogramService call (atlas: dendrogram.api.get-dendrogram-service singleton)')
      .toBe(true);
    expect(result.win1EqualsSvc1,
      'window.$dendrogramService === svc1 after the first call (referential equality)')
      .toBe(true);
    expect(result.svc2EqualsSvc1,
      'svc2 === svc1 (singleton instance is reused, not re-created on the second call)')
      .toBe(true);
    expect(result.win2EqualsSvc1,
      'window.$dendrogramService === svc1 continues to hold after the second call')
      .toBe(true);
    expect(result.svc1RespectsInitial,
      'when window.$dendrogramService pre-existed at Scenario 2 step 1, svc1 equals that cached instance (the package.ts#L84 guard returns the cache, does not overwrite)')
      .toBe(true);
    expect(result.svc1HasInject,
      'the singleton instance exposes the IDendrogramService surface (callable injectTreeForGrid)')
      .toBe(true);
  });

  // Cleanup per scenario .md Setup contract.
  await page.evaluate(() => {
    const closeBtn = document.querySelector('.dendrogram-close-bttn') as HTMLElement | null;
    if (closeBtn) closeBtn.click();
    grok.shell.closeAll();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
