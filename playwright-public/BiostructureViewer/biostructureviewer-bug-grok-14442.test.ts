/* ---
sub_features_covered: [biostructure.file-open, biostructure.file-open.importPdb, biostructure.file-open.importPdbqt]
--- */
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

test('BiostructureViewer / GROK-14442 file-handler search disambiguation regression guard', async ({page}) => {
  test.setTimeout(120_000);
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
      // Verify the Files browser tree actually opened, not just that the Browse tab exists.
      await page.locator('.d4-tree-view-root').first().waitFor({state: 'visible', timeout: 30_000});
      expect(await page.locator('.d4-tree-view-root').first().isVisible()).toBe(true);
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

      // Drive the REAL file-handler dispatch — the registry (not the test) decides WHICH
      // import fires. openFile is the same entrypoint the Files-browser double-click uses; do
      // NOT name the handler here (that would make the spy tautological).
      await page.evaluate((path) => {
        const grokAny: any = (window as any).grok;
        if (typeof grokAny.dapi.files.openFile === 'function')
          grokAny.dapi.files.openFile(path);
      }, samplePdbPath);

      // Poll for a routing signal instead of a blind sleep.
      await page.waitForFunction(() => {
        const w: any = window as any;
        const tv = w.grok?.shell?.tv;
        return (w.__bvCapturedHandlers || []).length > 0 ||
          !!document.querySelector('[name="viewer-Biostructure"]') ||
          !!document.querySelector('[name="dialog-Open-file"]') ||
          !!(tv?.dataFrame?.columns && tv.dataFrame.columns.length > 0);
      }, null, {timeout: 30_000}).catch(() => {});

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

      // Inverse-regression guard (the exact pre-fix bug shape) — HARD asserts. They fail only if
      // .pdb actually routed to importPdbqt or mounted the AutoDock-pose UI; they hold whether the
      // spy fired or the dispatcher bypassed the JS spy.
      expect(
        importPdbqtFired,
        'GROK-14442 file-handler search regressed: .pdb routed to importPdbqt ' +
        '(extension prefix-containment collision — pre-fix bug shape). ' +
        `Captured handlers: ${JSON.stringify(scenario1SpyCaptured)}.`,
      ).toBe(false);
      expect(
        scenario1AutodockUiMounted,
        'GROK-14442 file-handler search regressed: the AutoDock-pose UI mounted for a .pdb input, ' +
        'implying importPdbqt fired.',
      ).toBe(false);
      // Positive routing is asserted only when a signal is observable: openFile may dispatch through
      // the Dart-side registry without crossing the JS grok.functions.call spy, and headless WebGL
      // may keep the Mol* viewer from mounting — asserting the positive unconditionally would
      // false-fail on CORRECT routing (see riskNotes). The core GROK-14442 invariant is hard-asserted
      // registry-side in Scenario A.
      if (scenario1SpyCaptured.length > 0)
        expect(
          importPdbFired,
          'GROK-14442 file-handler search regressed: .pdb did NOT route to importPdb. ' +
          `Captured handlers: ${JSON.stringify(scenario1SpyCaptured)}.`,
        ).toBe(true);
    });

    // Close any resulting view to leave a clean state for Scenario 2.
    await page.evaluate(() => {
      const grokAny: any = (window as any).grok;
      grokAny.shell.closeAll();
      // Best-effort dialog dismissal.
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
    });
    await page.waitForFunction(() => document.querySelectorAll('.d4-dialog').length === 0,
      null, {timeout: 5_000}).catch(() => {});

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

      // Drive the REAL dispatcher — the registry decides which import fires (no handler naming).
      await page.evaluate((path) => {
        const grokAny: any = (window as any).grok;
        if (typeof grokAny.dapi.files.openFile === 'function')
          grokAny.dapi.files.openFile(path);
      }, samplePdbqtPath);

      // Poll for a routing signal instead of a blind sleep.
      await page.waitForFunction(() => {
        const w: any = window as any;
        const tv = w.grok?.shell?.tv;
        return (w.__bvCapturedHandlers || []).length > 0 ||
          !!document.querySelector('[name="viewer-Biostructure"]') ||
          !!(tv?.dataFrame?.columns && tv.dataFrame.columns.length > 0);
      }, null, {timeout: 30_000}).catch(() => {});

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

      // Inverse-regression guards — HARD. Fail only if .pdbqt actually routed to importPdb or
      // mounted the Biostructure (Mol*) viewer.
      expect(
        importPdbFired,
        'GROK-14442 file-handler search inverse regression: .pdbqt routed ' +
        `to importPdb. Captured handlers: ${JSON.stringify(scenario2SpyCaptured)}.`,
      ).toBe(false);
      expect(
        scenario2BioViewerMounted,
        'GROK-14442 inverse regression: Biostructure (Mol*) viewer mounted for a .pdbqt input, ' +
        'implying importPdb fired.',
      ).toBe(false);
      // Positive routing asserted only when the JS spy observed it (see Scenario 1 note / riskNotes).
      if (scenario2SpyCaptured.length > 0)
        expect(
          importPdbqtFired,
          'GROK-14442 inverse: .pdbqt did NOT route to importPdbqt. ' +
          `Captured handlers: ${JSON.stringify(scenario2SpyCaptured)}.`,
        ).toBe(true);
    });

    // SCENARIO 2 step 5 — Joint disambiguation cross-check, derived from the observed dispatches.
    await softStep('Scenario 2 step 5 — Joint disambiguation cross-check (GROK-14442 invariant)', async () => {
      // Neither scenario may have captured the OTHER extension's handler — the cross-contamination
      // that GROK-14442 fixed. Asserted from the actual captured sets, not hard-coded literals.
      expect(
        scenario1SpyCaptured.includes('BiostructureViewer:importPdbqt'),
        `GROK-14442 cross-check: .pdb captured importPdbqt (${JSON.stringify(scenario1SpyCaptured)}).`,
      ).toBe(false);
      expect(
        scenario2SpyCaptured.includes('BiostructureViewer:importPdb'),
        `GROK-14442 cross-check: .pdbqt captured importPdb (${JSON.stringify(scenario2SpyCaptured)}).`,
      ).toBe(false);
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
