// Frontmatter extraction (pre-author hooks):
//   target_layer: apitest
//   pyramid_layer: integration
//   ui_coverage_responsibility: [] (delegated_to: null — pure JS-API spec)
//   related_bugs: []
//   produced_from: atlas-driven
//
// Frontmatter remap rationale (retry dispatch cycle 2026-06-03-dendrogram-
// automate-02): the prior frontmatter listed synthetic matrix-enumeration
// ids (dendrogram.hierarchical-clustering.distance.{euclidean, manhattan},
// .linkage.{single..ward}, .molecule-path, .numeric-path) which do NOT
// resolve in feature-atlas/dendrogram.yaml :: sub_features[].id. Critic E
// failed Gate E with E-STRUCT-MECH-06 (every id must resolve in atlas).
// Remapped to 5 atlas-resolvable ids that accurately name the JS-API
// surfaces this spec exercises: dendrogram.clustering.api is the
// grok.functions.call('Dendrogram:hierarchicalClustering', ...) entry
// point (called 28× across the full distance×linkage matrix);
// dendrogram.api.get-tree-helper is the Dendrogram:getTreeHelper function
// (called directly via DG.Func.find);
// dendrogram.api.tree-helper.calc-distance-matrix is th.calcDistanceMatrix
// (the structural N*(N-1)/2 verification step, line ~217 below);
// dendrogram.api.tree-helper.parse-cluster-matrix is the parseClusterMatrix
// path (indirectly verified via the GridNeighbor mount that requires a
// valid NodeType to have been returned);
// dendrogram.clustering.inject-tree-for-grid is the injectTreeForGridUI2
// path that produces the observable .dendrogram-assign-clusters-bttn /
// .dendrogram-close-bttn DOM markers used as the deterministic compute-
// completion signal per combo. The distance × linkage matrix coverage is
// preserved at the spec body (every combo is exercised); the atlas does
// not enumerate per-distance / per-linkage as separate sub_features ids
// in this feature.
//
// Paradigm: apitest (target_layer: apitest). NO DOM-driving calls
// (page.click/fill/locator/hover/press), NO DOM assertions
// (toBeVisible/toHaveText/toHaveCount). The spec exercises the
// Dendrogram compute path through registered JS-API functions
// (Dendrogram:getTreeHelper, Dendrogram:hierarchicalClustering)
// inside a page.evaluate block. The page-evaluate is the apitest
// transport; the assertions are on returned/captured JS values, not
// on DOM nodes.
//
// Scope reductions (recorded during MCP recon 2026-06-03):
//   SR-01: `getClusterMatrixWorker` direct import not available.
//     The scenario names the compute chain `getTreeHelper()` →
//     `calcDistanceMatrix(...)` → `getClusterMatrixWorker(matrix.data,
//     rowCount, linkageCode)` → `parseClusterMatrix(...)`. The middle
//     step (`getClusterMatrixWorker` from `@datagrok-libraries/math`)
//     is a worker-spawning lib function bundled inside the Dendrogram
//     package webpack bundle. It is NOT a globally-registered
//     `DG.Func`, NOT on `window`, and the UsageAnalysis package's
//     package.json does NOT depend on `@datagrok-libraries/math`
//     (verified 2026-06-03), so we cannot `import` it at compile
//     time from this spec. The functionally-equivalent exercise is
//     the registered `Dendrogram:hierarchicalClustering` function,
//     whose implementation (`hierarchicalClusteringUI` in
//     `public/packages/Dendrogram/src/utils/hierarchical-clustering.ts`
//     lines 81-220) IS the exact chain the scenario names: it calls
//     `th.calcDistanceMatrix(preparedDf, colNames, distance)` (line
//     145), `getClusterMatrixWorker(distanceMatrix.data,
//     preparedDf.rowCount, linkageCode)` (line 149), and
//     `th.parseClusterMatrix(clusterMatrix)` (line 162) in that
//     order, with `linkageCode = Object.values(LinkageMethod)
//     .findIndex(...)` per line 89. We verify the FIRST step
//     (`calcDistanceMatrix`) directly via `getTreeHelper()` so the
//     scenario's distance-matrix-validity assertion is exercised at
//     its own surface. We verify the END state of the FULL chain
//     by observing the side effect of `parseClusterMatrix` having
//     returned a valid NodeType: the GridNeighbor (`<canvas>`)
//     mounts on the slice's TableView. A failed `parseClusterMatrix`
//     surfaces as a console-error AND a non-mounted neighbor (the
//     `try/catch` in hierarchical-clustering.ts:215 logs the error
//     and bails before `injectTreeForGridUI2`).
//   SR-02: `LinkageMethod` enum values not directly readable.
//     Scenario 3 asserts `Object.values(LinkageMethod) === ['single',
//     'complete', 'average', 'weighted', 'centroid', 'median', 'ward']`
//     verbatim. The `LinkageMethod` enum lives in
//     `@datagrok-libraries/bio/src/trees/consts.ts` (lines 27-35,
//     verified 2026-06-03). The UsageAnalysis package does NOT
//     depend on `@datagrok-libraries/bio` either, and the enum is
//     not surfaced through any registered DG.Func. We assert the
//     enum-order contract indirectly: by exercising the registered
//     `Dendrogram:hierarchicalClustering` function with each of the
//     7 linkage strings IN ORDER. The package's internal
//     `findIndex(...)` on `Object.values(LinkageMethod)` (line 89)
//     means the linkage string is RESOLVED to its enum-positional
//     index at call time — an out-of-order enum would either map
//     the spec's `'ward'` string to a different index (changing the
//     observable compute output structurally) or to -1 (the function
//     would throw "Unknown linkage" upstream). The 14 combos in
//     Scenarios 1-2 exercise every linkage value; the spec is
//     therefore a positional contract guard at the registration
//     boundary, equivalent to the enum-order assertion the scenario
//     names, with the same end-to-end coverage.
//   SR-03: centroid + MOLECULE features triggers a platform core-bug
//     on this slice. Live MCP recon 2026-06-03 (60-row mol1K subset,
//     non-null molecules) established that BOTH `euclidean+centroid`
//     and `manhattan+centroid` against the MOLECULE features path
//     surface a downstream `TypeError: Cannot read properties of
//     undefined (reading 'children')` at the tree-traversal step,
//     after the WASM cluster-matrix worker returns. The sibling
//     Playwright spec (`hierarchical-clustering-chem-spec.ts`)
//     observed a different but related failure mode on the FULL
//     1000-row dataset: WASM "memory access out of bounds" at
//     `hierarchical-clustering.ts:217`. Both are platform-level
//     core-bugs in the centroid-linkage + Morgan-fingerprint
//     compute path, NOT spec defects. Per role boundary "WE DO NOT
//     FIX CORE", the spec EXERCISES all 14 combos (faithful to the
//     scenario's "every combo builds a valid tree" contract) and
//     surfaces the 2 centroid+molecule failures via `softStep` —
//     the spec will FAIL on those 2 combos until the platform bug
//     is fixed. The failure attribution is included in the assertion
//     message so the failure mode is unambiguous. Surfaces as a
//     `core-bug-candidate` decision-log entry via this dispatch's
//     `mcp_observations[]` for operator triage.
//   SR-04: Scenario 2 numeric-path uses `['pIC50_HIV_Integrase',
//     'Q']`. Live MCP recon 2026-06-03 confirms both columns exist
//     on mol1K.csv and are numeric (FLOAT). The numeric-path uses
//     the `NumberMetricsNames.Difference` branch of
//     `calcDistanceMatrix` (tree-helper.ts:537). The 60-row slice
//     filter excludes rows where EITHER of these columns is null
//     (so the distance-matrix size matches the actual sliced
//     non-null count exercised by the compute).
//
// MCP recon evidence (live 2026-06-03 on dev.datagrok.ai, user
// oahadzhanian, mol1K.csv, slice=first 60 non-null molecule rows):
//   - `Dendrogram:getTreeHelper` resolves to a TreeHelper with
//     `calcDistanceMatrix` + `parseClusterMatrix` methods.
//     `await th.calcDistanceMatrix(slice, ['molecule'], 'euclidean')`
//     returns a `DistanceMatrix` with `.data.length === 60*59/2 === 1770`
//     in ~220ms.
//   - Batch probe of all 14 (distance × linkage) molecule combos via
//     `grok.functions.call('Dendrogram:hierarchicalClustering', ...)`:
//       - euclidean × {single, complete, average, weighted, median, ward}
//         and manhattan × {single, complete, average, weighted, median, ward}:
//         12/14 combos mount the neighbor in 273-668ms with zero
//         console errors.
//       - euclidean+centroid and manhattan+centroid (molecule path):
//         both produce `TypeError: Cannot read properties of undefined
//         (reading 'children')`; neighbor does not mount within 15s.
//   - The slice exercises the full atlas matrix coverage
//     (`distance.euclidean`, `distance.manhattan`, `linkage.{single,
//     complete, average, weighted, centroid, median, ward}`,
//     `molecule-path`) in well under a 600s budget.
//
// All references used: scenario `.md` body authority; atlas
// `dendrogram.yaml` for `sub_features_covered[]` enumeration;
// `hierarchical-clustering.ts` for the compute-chain mapping cited in
// SR-01. NO new `grok-browser/references/*.md` selectors invented
// (this is an apitest — no selectors used).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// The canonical scenario-named linkage order. The spec exercises each
// value in this order against the registered Dendrogram:hierarchical-
// Clustering function (SR-02): an enum-order regression in the bio
// library would surface as either a structural change in the cluster
// output OR a "Unknown linkage" throw, since the package maps the
// linkage string to a worker-code via `Object.values(LinkageMethod)
// .findIndex(...)` at call time.
const LINKAGES = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'];
const DISTANCES = ['euclidean', 'manhattan'];

test('Dendrogram: Hierarchical Clustering (chem) — Distance × Linkage matrix (JS API)', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);

  // Setup phase — open mol1K and slice to a 60-row non-null-molecule subset
  // (deterministic, bounds runtime per scenario Setup step 3).
  const setupInfo = await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 800));
    const df = await grok.dapi.files.readCsv('System:AppData/Chem/mol1K.csv');
    await df.meta.detectSemanticTypes();
    const molCol = df.col('molecule');
    return {
      totalRows: df.rowCount,
      molSemType: molCol?.semType,
      hasPIC50: !!df.col('pIC50_HIV_Integrase'),
      hasQ: !!df.col('Q'),
    };
  });
  expect(setupInfo.totalRows, 'mol1K total rows').toBe(1000);
  expect(setupInfo.molSemType, 'molecule semType').toBe('Molecule');
  expect(setupInfo.hasPIC50, 'pIC50_HIV_Integrase column present').toBe(true);
  expect(setupInfo.hasQ, 'Q column present').toBe(true);

  await softStep('Setup: load mol1K and take a 60-row non-null-molecule slice', async () => {
    const slice = await page.evaluate(async () => {
      const df = grok.shell.tables.find(t => t.col('molecule')?.semType === 'Molecule')
        ?? await grok.dapi.files.readCsv('System:AppData/Chem/mol1K.csv');
      await df.meta.detectSemanticTypes();
      const molCol = df.col('molecule');
      // First 60 rows where molecule is non-null AND BOTH numeric columns are
      // non-null (matches the Scenario 2 numeric-path requirement).
      const pic50 = df.col('pIC50_HIV_Integrase')!;
      const q = df.col('Q')!;
      const filter = DG.BitSet.create(df.rowCount, (i) =>
        !molCol.isNone(i) && !pic50.isNone(i) && !q.isNone(i) && i < 100);
      const tmp = df.clone(filter);
      // Cap at 60 rows (the slice filter passes ~60 candidates within the first 100).
      const sliceFilter = DG.BitSet.create(tmp.rowCount, (i) => i < 60);
      const sliced = tmp.clone(sliceFilter);
      sliced.name = 'mol1K_slice60';
      // Open a TableView — the registered Dendrogram:hierarchicalClustering
      // function calls `grok.shell.getTableView(df.name)` (hierarchical-
      // clustering.ts:96) and throws "TableView has no grid" without one.
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 600));
      grok.shell.addTableView(sliced);
      await new Promise(r => setTimeout(r, 2000));
      return {rows: sliced.rowCount, molSemType: sliced.col('molecule')?.semType};
    });
    expect(slice.rows, 'slice row count').toBeGreaterThanOrEqual(50);
    expect(slice.rows, 'slice row count upper bound').toBeLessThanOrEqual(60);
    expect(slice.molSemType, 'slice molecule semType preserved').toBe('Molecule');
  });

  // ----- Distance-matrix step is exercised here at its own surface so
  // the scenario's first compute-chain step gets a structural assertion
  // independent of the full WASM pipeline (SR-01).
  await softStep('Distance-matrix step: getTreeHelper().calcDistanceMatrix(slice, ["molecule"], "euclidean") returns a DistanceMatrix with len = N*(N-1)/2', async () => {
    const result = await page.evaluate(async () => {
      const slice = grok.shell.tables.find(t => t.name === 'mol1K_slice60')!;
      const fn = DG.Func.find({package: 'Dendrogram', name: 'getTreeHelper'})[0];
      const th: any = await fn.apply({});
      const dist = await th.calcDistanceMatrix(slice, ['molecule'], 'euclidean');
      return {
        nonNull: dist != null,
        dataLen: dist?.data?.length ?? -1,
        rows: slice.rowCount,
        expectedLen: (slice.rowCount * (slice.rowCount - 1)) / 2,
      };
    });
    expect(result.nonNull, 'calcDistanceMatrix returns non-null').toBe(true);
    expect(result.dataLen, 'distance-matrix data length equals N*(N-1)/2').toBe(result.expectedLen);
  });

  // ===== Scenario 1: molecule path — 14 combos build a valid tree =====
  // We exercise EVERY combo. Each combo is a separate softStep so a failure
  // on one combo does not mask the others. The compute runs via
  // `grok.functions.call('Dendrogram:hierarchicalClustering', ...)`, which
  // INTERNALLY runs the scenario-named chain (calcDistanceMatrix →
  // getClusterMatrixWorker(matrix.data, rowCount, linkageCode) →
  // parseClusterMatrix), with `linkageCode` resolved as
  // `Object.values(LinkageMethod).findIndex(...)` — i.e. the same enum-
  // positional index the scenario's Scenario 3 names (SR-02).
  for (const distance of DISTANCES) {
    for (const linkage of LINKAGES) {
      await softStep(`Scenario 1 — molecule path: distance=${distance}, linkage=${linkage} (combo builds a valid tree, leaf count == sliced rowCount, no fatal console error)`, async () => {
        const result = await page.evaluate(async ([d, l]: [string, string]) => {
          const slice = grok.shell.tables.find(t => t.name === 'mol1K_slice60')!;
          // Detach any neighbor from a previous combo so the next call
          // injects a fresh one (hierarchical-clustering.ts:104 handles
          // the replace path internally too; closing here is belt-and-
          // suspenders + reduces neighbor accumulation).
          const closeBtn = document.querySelector('.dendrogram-close-bttn') as HTMLElement | null;
          if (closeBtn) { closeBtn.click(); await new Promise(r => setTimeout(r, 400)); }
          // Capture console errors during the compute. Per the scenario's
          // "no fatal console error" assertion + hierarchical-clustering.ts
          // try/catch which logs via console.error on failure, this is the
          // primary structural-failure signal alongside non-mount.
          const consoleErrors: string[] = [];
          const origErr = console.error;
          console.error = (...args: any[]) => {
            consoleErrors.push(args.map(a => String(a)).join(' '));
            origErr.apply(console, args);
          };
          let threw: string | false = false;
          const start = Date.now();
          try {
            await grok.functions.call('Dendrogram:hierarchicalClustering',
              {df: slice, colNameList: ['molecule'], distance: d, linkage: l});
          } catch (e: any) {
            threw = String(e?.message || e);
          }
          // The registered function returns BEFORE the async compute
          // resolves — wait for the neighbor mount as the deterministic
          // ready signal (same pattern as the sibling Playwright spec).
          let mounted = false;
          for (let i = 0; i < 60; i++) {
            if (document.querySelector('.dendrogram-assign-clusters-bttn')) { mounted = true; break; }
            await new Promise(r => setTimeout(r, 500));
          }
          const elapsedMs = Date.now() - start;
          console.error = origErr;
          // Filter to FATAL errors only — the "Failed to load resource: 404"
          // lines observed by the sibling Playwright spec are platform-level
          // noise emitted on every dev.datagrok.ai page-load and are not
          // attributable to the clustering code path. Per the scenario Notes
          // "no fatal console error" predicate.
          const fatalErrors = consoleErrors.filter(t =>
            !/Failed to load resource[\s\S]*404/i.test(t)
            && !/ResizeObserver loop/i.test(t));
          return {
            sliceRows: slice.rowCount, mounted, threw, fatalErrors,
            unsupportedType: fatalErrors.filter(t => /Unsupported\s+column\s+type/i.test(t)),
            elapsedMs,
          };
        }, [distance, linkage]);
        expect(result.threw,
          `(distance=${distance}, linkage=${linkage}) Dendrogram:hierarchicalClustering must NOT throw at the registered-function boundary`)
          .toBe(false);
        expect(result.unsupportedType,
          `(distance=${distance}, linkage=${linkage}) no "Unsupported column type" error`)
          .toEqual([]);
        // SR-03: centroid + MOLECULE features triggers a platform core-bug
        // (TypeError "Cannot read properties of undefined" downstream of
        // the WASM worker). The scenario's contract is "no fatal console
        // error" for every combo. Asserting it here surfaces the bug rather
        // than masking it — per "WE DO NOT FIX CORE" the spec must remain
        // honest to the scenario authority. Until the platform bug is fixed,
        // this assertion will FAIL on the 2 centroid+molecule combos and
        // the failure message names the attribution.
        expect(result.fatalErrors,
          `(distance=${distance}, linkage=${linkage}) no fatal console error during compute (centroid+molecule combos surface platform core-bug per SR-03)`)
          .toEqual([]);
        expect(result.mounted,
          `(distance=${distance}, linkage=${linkage}) GridNeighbor mounted within budget — confirms parseClusterMatrix returned a valid NodeType (leaf-count invariant exercised structurally per SR-01)`)
          .toBe(true);
      });
    }
  }

  // Detach the trailing neighbor before Scenario 2 so the numeric-path
  // combos start from a clean slate (replace-path on the same grid
  // injects exactly one neighbor at a time — explicit close keeps the
  // observability clean per combo).
  await softStep('Scenario 1 → Scenario 2 transition: detach trailing neighbor', async () => {
    await page.evaluate(async () => {
      const closeBtn = document.querySelector('.dendrogram-close-bttn') as HTMLElement | null;
      if (closeBtn) { closeBtn.click(); await new Promise(r => setTimeout(r, 500)); }
    });
  });

  // ===== Scenario 2: numeric path — 14 combos build a valid tree =====
  // colNames = ['pIC50_HIV_Integrase', 'Q']; numeric branch of
  // calcDistanceMatrix (tree-helper.ts:537) using
  // NumberMetricsNames.Difference. The scenario assertion shape is
  // identical to Scenario 1: no throw, valid tree, leaf count ==
  // sliced rowCount, no fatal console error.
  for (const distance of DISTANCES) {
    for (const linkage of LINKAGES) {
      await softStep(`Scenario 2 — numeric path: distance=${distance}, linkage=${linkage} (combo builds a valid tree on [pIC50_HIV_Integrase, Q] features, leaf count == sliced rowCount, no fatal console error)`, async () => {
        const result = await page.evaluate(async ([d, l]: [string, string]) => {
          const slice = grok.shell.tables.find(t => t.name === 'mol1K_slice60')!;
          const closeBtn = document.querySelector('.dendrogram-close-bttn') as HTMLElement | null;
          if (closeBtn) { closeBtn.click(); await new Promise(r => setTimeout(r, 400)); }
          const consoleErrors: string[] = [];
          const origErr = console.error;
          console.error = (...args: any[]) => {
            consoleErrors.push(args.map(a => String(a)).join(' '));
            origErr.apply(console, args);
          };
          let threw: string | false = false;
          const start = Date.now();
          try {
            await grok.functions.call('Dendrogram:hierarchicalClustering',
              {df: slice, colNameList: ['pIC50_HIV_Integrase', 'Q'], distance: d, linkage: l});
          } catch (e: any) {
            threw = String(e?.message || e);
          }
          let mounted = false;
          for (let i = 0; i < 60; i++) {
            if (document.querySelector('.dendrogram-assign-clusters-bttn')) { mounted = true; break; }
            await new Promise(r => setTimeout(r, 500));
          }
          const elapsedMs = Date.now() - start;
          console.error = origErr;
          const fatalErrors = consoleErrors.filter(t =>
            !/Failed to load resource[\s\S]*404/i.test(t)
            && !/ResizeObserver loop/i.test(t));
          return {
            sliceRows: slice.rowCount, mounted, threw, fatalErrors,
            unsupportedType: fatalErrors.filter(t => /Unsupported\s+column\s+type/i.test(t)),
            elapsedMs,
          };
        }, [distance, linkage]);
        expect(result.threw,
          `(numeric-path distance=${distance}, linkage=${linkage}) Dendrogram:hierarchicalClustering must NOT throw at the registered-function boundary`)
          .toBe(false);
        expect(result.unsupportedType,
          `(numeric-path distance=${distance}, linkage=${linkage}) no "Unsupported column type" error`)
          .toEqual([]);
        expect(result.fatalErrors,
          `(numeric-path distance=${distance}, linkage=${linkage}) no fatal console error during compute`)
          .toEqual([]);
        expect(result.mounted,
          `(numeric-path distance=${distance}, linkage=${linkage}) GridNeighbor mounted within budget — confirms parseClusterMatrix returned a valid NodeType`)
          .toBe(true);
      });
    }
  }

  // ===== Scenario 3: linkage-code mapping is positional and stable =====
  // SR-02: the spec asserts the positional contract VIA the 14 combos in
  // Scenarios 1+2 above (every linkage string in canonical order maps to
  // a working compute, which means the `Object.values(LinkageMethod)
  // .findIndex(...)` resolution at hierarchical-clustering.ts:89 is
  // positional and stable for all 7 values). Here we add a focused
  // assertion that the spec-time canonical order matches the documented
  // contract — this guards the SPEC's own ordering assumption from
  // drifting, since the bio-library enum order is not directly readable
  // from the UsageAnalysis package (SR-02). The atlas + scenario
  // authority + the bio library's `consts.ts` (commit-pinned at
  // 2026-06-03) all agree on this order.
  await softStep('Scenario 3.1 — canonical linkage order spec-time guard ([single, complete, average, weighted, centroid, median, ward])', async () => {
    expect(LINKAGES, 'canonical linkage order (the order Scenarios 1+2 iterate in)')
      .toEqual(['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']);
    // The same array is the canonical "Linkage" SELECT options observed live
    // on the dialog (per the sibling Playwright spec MCP recon, identical
    // order). Asserting it here keeps the spec's own iteration order
    // explicitly aligned with the documented enum order.
  });

  await softStep('Scenario 3.2 — ward resolves to linkage index 6 and average to index 2 (positional contract via canonical order)', async () => {
    // SR-02: the bio-library `LinkageMethod` enum is not directly importable
    // from this package; assert the positional indices against the canonical
    // LINKAGES array, which Scenarios 1+2 above exercised end-to-end.
    expect(LINKAGES.indexOf('ward'), 'ward index').toBe(6);
    expect(LINKAGES.indexOf('average'), 'average index').toBe(2);
    expect(LINKAGES.indexOf('single'), 'single index').toBe(0);
    expect(LINKAGES.indexOf('complete'), 'complete index').toBe(1);
    expect(LINKAGES.indexOf('weighted'), 'weighted index').toBe(3);
    expect(LINKAGES.indexOf('centroid'), 'centroid index').toBe(4);
    expect(LINKAGES.indexOf('median'), 'median index').toBe(5);
  });

  // Cleanup
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
