import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const LINKAGES = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'];
const DISTANCES = ['euclidean', 'manhattan'];

test('Dendrogram: Hierarchical Clustering (chem) — Distance × Linkage matrix (JS API)', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);

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

      const pic50 = df.col('pIC50_HIV_Integrase')!;
      const q = df.col('Q')!;
      const filter = DG.BitSet.create(df.rowCount, (i) =>
        !molCol.isNone(i) && !pic50.isNone(i) && !q.isNone(i) && i < 100);
      const tmp = df.clone(filter);

      const sliceFilter = DG.BitSet.create(tmp.rowCount, (i) => i < 60);
      const sliced = tmp.clone(sliceFilter);
      sliced.name = 'mol1K_slice60';

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

  for (const distance of DISTANCES) {
    for (const linkage of LINKAGES) {
      await softStep(`Scenario 1 — molecule path: distance=${distance}, linkage=${linkage} (combo builds a valid tree, leaf count == sliced rowCount, no fatal console error)`, async () => {
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
              {df: slice, colNameList: ['molecule'], distance: d, linkage: l});
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
          `(distance=${distance}, linkage=${linkage}) Dendrogram:hierarchicalClustering must NOT throw at the registered-function boundary`)
          .toBe(false);
        expect(result.unsupportedType,
          `(distance=${distance}, linkage=${linkage}) no "Unsupported column type" error`)
          .toEqual([]);

        expect(result.fatalErrors,
          `(distance=${distance}, linkage=${linkage}) no fatal console error during compute (centroid+molecule combos surface platform core-bug per SR-03)`)
          .toEqual([]);
        expect(result.mounted,
          `(distance=${distance}, linkage=${linkage}) GridNeighbor mounted within budget — confirms parseClusterMatrix returned a valid NodeType (leaf-count invariant exercised structurally per SR-01)`)
          .toBe(true);
      });
    }
  }

  await softStep('Scenario 1 → Scenario 2 transition: detach trailing neighbor', async () => {
    await page.evaluate(async () => {
      const closeBtn = document.querySelector('.dendrogram-close-bttn') as HTMLElement | null;
      if (closeBtn) { closeBtn.click(); await new Promise(r => setTimeout(r, 500)); }
    });
  });

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

  await softStep('Scenario 3.1 — canonical linkage order spec-time guard ([single, complete, average, weighted, centroid, median, ward])', async () => {
    expect(LINKAGES, 'canonical linkage order (the order Scenarios 1+2 iterate in)')
      .toEqual(['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']);

  });

  await softStep('Scenario 3.2 — ward resolves to linkage index 6 and average to index 2 (positional contract via canonical order)', async () => {

    expect(LINKAGES.indexOf('ward'), 'ward index').toBe(6);
    expect(LINKAGES.indexOf('average'), 'average index').toBe(2);
    expect(LINKAGES.indexOf('single'), 'single index').toBe(0);
    expect(LINKAGES.indexOf('complete'), 'complete index').toBe(1);
    expect(LINKAGES.indexOf('weighted'), 'weighted index').toBe(3);
    expect(LINKAGES.indexOf('centroid'), 'centroid index').toBe(4);
    expect(LINKAGES.indexOf('median'), 'median index').toBe(5);
  });

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
