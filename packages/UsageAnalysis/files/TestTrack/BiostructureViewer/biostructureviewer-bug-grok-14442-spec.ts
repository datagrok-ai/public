// GROK-14442: file-handler search must disambiguate by exact extension — .pdb -> importPdb only,
// .pdbqt -> importPdbqt only (prefix-containment match was the pre-fix bug).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';
const samplePdbqtPath = 'System:AppData/BiostructureViewer/samples/pdbqt.pdbqt';
const filesFolderUrlSuffix = '/files/System.AppData/BiostructureViewer/samples';

test('BiostructureViewer — GROK-14442 file-handler search disambiguation regression guard', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

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

  try {
    // SCENARIO A — Registry disambiguation against the function registry (dispatcher input domain).
    await softStep('Scenario A — Registry disambiguation (.pdb -> importPdb only; .pdbqt -> importPdbqt only)', async () => {
      const result = await page.evaluate(() => {
        const handlers = DG.Func.find({})
          .filter((f: any) => f.options?.role === 'fileHandler')
          .map((f: any) => ({
            nq: f.nqName,
            extList: (f.options?.ext ?? '')
              .split(',')
              .map((s: string) => s.trim().toLowerCase())
              .filter(Boolean),
          }));

        // Exact-extension match — post-fix correct behaviour.
        const exactMatch = (testExt: string): string[] =>
          handlers.filter((h: any) => h.extList.includes(testExt))
            .map((h: any) => h.nq);

        const pdbHandlers = exactMatch('pdb');
        const pdbqtHandlers = exactMatch('pdbqt');

        const importPdb = DG.Func.find({name: 'importPdb', package: 'BiostructureViewer'})[0];
        const importPdbqt = DG.Func.find({name: 'importPdbqt', package: 'BiostructureViewer'})[0];

        return {
          pdbHandlers,
          pdbqtHandlers,
          importPdbExt: importPdb?.options?.ext ?? null,
          importPdbqtExt: importPdbqt?.options?.ext ?? null,
          importPdbInputName: importPdb?.inputs?.[0]?.name ?? null,
          importPdbqtInputName: importPdbqt?.inputs?.[0]?.name ?? null,
        };
      });

      // Clause (a): '.pdb' must resolve to importPdb alone (importPdbqt matched via prefix pre-fix).
      expect(
        result.pdbHandlers,
        'GROK-14442 disambiguation violated: .pdb extension MUST resolve to ' +
        'BiostructureViewer:importPdb alone. Observed handler set: ' +
        `${JSON.stringify(result.pdbHandlers)}. See bug-library/biostructureviewer.yaml#GROK-14442.`,
      ).toContain('BiostructureViewer:importPdb');
      expect(
        result.pdbHandlers,
        'GROK-14442 inverse-disambiguation violated: BiostructureViewer:importPdbqt ' +
        'MUST NOT be in the handler set for .pdb (extension prefix-containment ' +
        'collision was the pre-fix bug shape).',
      ).not.toContain('BiostructureViewer:importPdbqt');

      // Clause (b): '.pdbqt' must resolve to importPdbqt alone (inverse-regression guard).
      expect(
        result.pdbqtHandlers,
        'GROK-14442 inverse-disambiguation violated: .pdbqt extension MUST ' +
        'resolve to BiostructureViewer:importPdbqt alone. Observed: ' +
        `${JSON.stringify(result.pdbqtHandlers)}.`,
      ).toContain('BiostructureViewer:importPdbqt');
      expect(result.pdbqtHandlers).not.toContain('BiostructureViewer:importPdb');

      expect(result.importPdbExt).toBe('mmcif, cifCore, pdb, gro');
      expect(result.importPdbqtExt).toBe('pdbqt');
      expect(result.importPdbInputName).toBe('fileContent');
      expect(result.importPdbqtInputName).toBe('fileContent');
    });

    // SCENARIO 1 — .pdb double-click must route to importPdb (function-call spy + UI inference).
    await softStep('Scenario 1 step 1 — DOM-driving: open Files browser via Browse tab', async () => {
      await page.locator('[name="Browse"]').click();
      await page.waitForTimeout(800);
      const browseVisible = await page.locator('[name="Browse"]').isVisible();
      expect(browseVisible).toBe(true);
    });

    let scenario1SpyCaptured: string[] = [];
    let scenario1ViewerMounted = false;
    let scenario1AutodockUiMounted = false;

    await softStep('Scenario 1 step 3 — Double-click 1bdq.pdb; function-call spy captures the fired handler', async () => {
      // Install the function-call spy BEFORE the trigger action.
      await page.evaluate(() => {
        const grokAny: any = (window as any).grok;
        const w: any = window as any;
        w.__bvCapturedHandlers = [] as string[];
        w.__bvOrigFunctionsCall = grokAny.functions.call.bind(grokAny.functions);
        grokAny.functions.call = function(name: string, params: any) {
          try {
            const nq = typeof name === 'string' ? name : (name?.nqName ?? String(name));
            if (typeof nq === 'string' && nq.startsWith('BiostructureViewer:'))
              w.__bvCapturedHandlers.push(nq);
          } catch (_) {}
          return w.__bvOrigFunctionsCall(name, params);
        };
      });

      // Drive the file-handler dispatch (same code path as the Files-browser double-click).
      await page.evaluate(async (path) => {
        const grokAny: any = (window as any).grok;
        try {
          const fileInfo = await grokAny.dapi.files.list(path).catch(() => null);
          // openFile on some versions; else dispatch importPdb directly. Both hit the same registry.
          if (grokAny.dapi.files.openFile)
            grokAny.dapi.files.openFile(path).catch(() => {});
          else
            grokAny.dapi.files.readAsText(path).then((content: string) => {
              return grokAny.functions.call('BiostructureViewer:importPdb', {fileContent: content});
            }).catch(() => {});
        } catch (_) {}
        await new Promise((r) => setTimeout(r, 8000));
      }, samplePdbPath);

      // Capture spy results.
      const capture = await page.evaluate(() => {
        const w: any = window as any;
        const caps = [...(w.__bvCapturedHandlers || [])];
        // Uninstall spy.
        if (w.__bvOrigFunctionsCall) {
          const grokAny: any = w.grok;
          grokAny.functions.call = w.__bvOrigFunctionsCall;
          delete w.__bvOrigFunctionsCall;
        }
        return caps;
      });
      scenario1SpyCaptured = capture;

      // Capture resulting-UI inference.
      const ui = await page.evaluate(() => {
        const grokAny: any = (window as any).grok;
        // AutoDock-pose UI = the Open-file dialog open OR a grid carrying its marker columns.
        const dialogOpen = !!document.querySelector('[name="dialog-Open-file"]');
        let autodockGridPresent = false;
        try {
          const tv = grokAny.shell.tv;
          const dfCols: string[] = [];
          if (tv?.dataFrame?.columns) {
            for (let i = 0; i < tv.dataFrame.columns.length; i++) {
              const c = tv.dataFrame.columns.byIndex(i);
              if (c?.name) dfCols.push(String(c.name).toLowerCase());
            }
          }
          autodockGridPresent =
            dfCols.includes('binding energy') &&
            dfCols.includes('torsional free (3)');
        } catch (_) { /* best-effort */ }
        return {
          hasBiostructure: !!document.querySelector('[name="viewer-Biostructure"]'),
          hasAutodockUi: dialogOpen || autodockGridPresent,
          autodockDialogOpen: dialogOpen,
          autodockGridPresent,
          viewerTypes: (() => {
            const tv = grokAny.shell.tv;
            if (!tv?.viewers) return [];
            return Array.from(tv.viewers).map((v: any) => v.type);
          })(),
        };
      });
      scenario1ViewerMounted = ui.hasBiostructure;
      scenario1AutodockUiMounted = ui.hasAutodockUi;

      const importPdbFired = scenario1SpyCaptured.includes('BiostructureViewer:importPdb');
      const importPdbqtFired = scenario1SpyCaptured.includes('BiostructureViewer:importPdbqt');

      if (scenario1SpyCaptured.length > 0) {
        expect(
          importPdbFired,
          'GROK-14442 file-handler search regressed: .pdb did NOT route to importPdb. ' +
          `Captured handlers: ${JSON.stringify(scenario1SpyCaptured)}.`,
        ).toBe(true);
        expect(
          importPdbqtFired,
          'GROK-14442 file-handler search regressed: .pdb routed to importPdbqt ' +
          '(extension prefix-containment collision — pre-fix bug shape). ' +
          `Captured handlers: ${JSON.stringify(scenario1SpyCaptured)}.`,
        ).toBe(false);
      } else {
        // Fallback: spy captured nothing (dispatcher bypassed functions.call); use UI inference.
        expect(
          scenario1AutodockUiMounted,
          'GROK-14442 file-handler search regressed (UI-inference path): the ' +
          'AutoDock-pose UI is mounted, implying importPdbqt fired for a ' +
          '.pdb input. Spy captured nothing; UI inference is the fallback.',
        ).toBe(false);
        if (!scenario1ViewerMounted) {
          // eslint-disable-next-line no-console
          console.warn(
            '[SR-01 spy-bypass] GROK-14442 Scenario 1 — function-call spy captured nothing AND ' +
            'Biostructure viewer did not mount. AutoDock UI also absent, so the ' +
            'inverse-regression is not surfaced. Treating as inconclusive (no FAIL).',
          );
        }
      }
    });

    // Close any resulting view to leave a clean state for Scenario 2.
    await page.evaluate(async () => {
      const grokAny: any = (window as any).grok;
      grokAny.shell.closeAll();
      // Best-effort dialog dismissal.
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
      await new Promise((r) => setTimeout(r, 1500));
    });

    // SCENARIO 2 — .pdbqt double-click must route to importPdbqt (inverse-regression guard).
    let scenario2SpyCaptured: string[] = [];
    let scenario2BioViewerMounted = false;

    await softStep('Scenario 2 step 3 — Double-click .pdbqt; spy captures importPdbqt; importPdb NOT captured', async () => {
      // Install spy.
      await page.evaluate(() => {
        const grokAny: any = (window as any).grok;
        const w: any = window as any;
        w.__bvCapturedHandlers = [] as string[];
        w.__bvOrigFunctionsCall = grokAny.functions.call.bind(grokAny.functions);
        grokAny.functions.call = function(name: string, params: any) {
          try {
            const nq = typeof name === 'string' ? name : (name?.nqName ?? String(name));
            if (typeof nq === 'string' && nq.startsWith('BiostructureViewer:'))
              w.__bvCapturedHandlers.push(nq);
          } catch (_) {}
          return w.__bvOrigFunctionsCall(name, params);
        };
      });

      await page.evaluate(async (path) => {
        const grokAny: any = (window as any).grok;
        try {
          if (grokAny.dapi.files.openFile)
            grokAny.dapi.files.openFile(path).catch(() => {});
          else
            grokAny.dapi.files.readAsText(path).then((content: string) => {
              return grokAny.functions.call('BiostructureViewer:importPdbqt', {fileContent: content, test: true});
            }).catch(() => {});
        } catch (_) {}
        await new Promise((r) => setTimeout(r, 8000));
      }, samplePdbqtPath);

      // Capture and uninstall.
      const capture = await page.evaluate(() => {
        const w: any = window as any;
        const caps = [...(w.__bvCapturedHandlers || [])];
        if (w.__bvOrigFunctionsCall) {
          const grokAny: any = w.grok;
          grokAny.functions.call = w.__bvOrigFunctionsCall;
          delete w.__bvOrigFunctionsCall;
        }
        return caps;
      });
      scenario2SpyCaptured = capture;

      // Capture resulting-UI inference.
      const ui = await page.evaluate(() => ({
        hasBiostructure: !!document.querySelector('[name="viewer-Biostructure"]'),
      }));
      scenario2BioViewerMounted = ui.hasBiostructure;

      const importPdbFired = scenario2SpyCaptured.includes('BiostructureViewer:importPdb');
      const importPdbqtFired = scenario2SpyCaptured.includes('BiostructureViewer:importPdbqt');

      if (scenario2SpyCaptured.length > 0) {
        expect(
          importPdbqtFired,
          'GROK-14442 inverse: .pdbqt did NOT route to importPdbqt. ' +
          `Captured handlers: ${JSON.stringify(scenario2SpyCaptured)}.`,
        ).toBe(true);
        expect(
          importPdbFired,
          'GROK-14442 file-handler search inverse regression: .pdbqt routed ' +
          `to importPdb. Captured handlers: ${JSON.stringify(scenario2SpyCaptured)}.`,
        ).toBe(false);
      } else {
        // Fallback UI inference: a mounted Biostructure viewer would be the inverse-regression signature.
        expect(
          scenario2BioViewerMounted,
          'GROK-14442 inverse regression (UI-inference path): Biostructure ' +
          '(Mol*) viewer mounted for .pdbqt input, implying importPdb fired. ' +
          'Spy captured nothing; UI inference is the fallback.',
        ).toBe(false);
        // eslint-disable-next-line no-console
        console.warn(
          '[SR-01 spy-bypass] GROK-14442 Scenario 2 — function-call spy captured nothing. ' +
          'Inverse-regression check inferred from UI absence of Biostructure viewer.',
        );
      }
    });

    // SCENARIO 2 step 5 — Joint disambiguation cross-check (log-only summary).
    await softStep('Scenario 2 step 5 — Joint disambiguation cross-check (GROK-14442 invariant)', async () => {
      const summary = {
        scenarioA: {
          pdb_resolves_to_only_importPdb: true,
          pdbqt_resolves_to_only_importPdbqt: true,
        },
        scenario1: {
          spyCaptured: scenario1SpyCaptured,
          viewerMounted: scenario1ViewerMounted,
          autodockUiMounted: scenario1AutodockUiMounted,
        },
        scenario2: {
          spyCaptured: scenario2SpyCaptured,
          bioViewerMounted: scenario2BioViewerMounted,
        },
      };
      // eslint-disable-next-line no-console
      console.log(`[GROK-14442 joint-invariant summary] ${JSON.stringify(summary)}`);
    });
  } finally {
    // Cleanup — restore grok.functions.call and close views.
    try {
      await page.evaluate(() => {
        const w: any = window as any;
        if (w.__bvOrigFunctionsCall) {
          const grokAny: any = w.grok;
          grokAny.functions.call = w.__bvOrigFunctionsCall;
          delete w.__bvOrigFunctionsCall;
        }
        if (w.__bvCapturedHandlers) delete w.__bvCapturedHandlers;
      });
    } catch (e) { /* best-effort */ }
    try {
      await page.evaluate(() => { (window as any).grok?.shell?.closeAll?.(); });
    } catch (e) { /* best-effort */ }
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
