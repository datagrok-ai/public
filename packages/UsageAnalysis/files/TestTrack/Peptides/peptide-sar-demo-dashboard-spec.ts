// Peptide SAR demo dashboard + Peptides app landing — entry-point smokes.
// Scenario 1 invokes the FASTA SAR demo (Peptides:macromoleculeSarFastaDemo); Scenario 2 invokes
// the landing View (Peptides:Peptides) and clicks its three demo buttons via real DOM.
// Notes: Peptides registers no #app, so the View is reached via grok.functions.call. Complex/HELM
// demos ship Macromolecule columns named "MSA"/"HELM" (not "AlignedSequence"), so the shape check
// is semType-by-name-agnostic. grok.shell.tableViews is Iterable (Array.from before .length);
// grok.shell.lastError is a Promise (must be awaited).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

test('Peptide SAR demo dashboard + Peptides app landing — entry-point smokes',
  async ({page}) => {
    test.setTimeout(360_000);

    await loginToDatagrok(page);

    // Clean shell, Windows mode, pre-warm Peptides:initPeptides (GROK-17557 cold-init race).
    await softStep('Setup: clean shell + pre-warm Peptides @init', async () => {
      const result = await page.evaluate(async () => {
        document.querySelectorAll('.d4-dialog').forEach((d) => {
          const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
          if (cancel) cancel.click();
        });
        grok.shell.closeAll();
        document.body.classList.add('selenium');
        grok.shell.settings.showFiltersIconsConstantly = true;
        grok.shell.windows.simpleMode = false;
        let initError: string | null = null;
        try { await grok.functions.call('Peptides:initPeptides'); }
        catch (e) { initError = String(e); }
        return {
          initError,
          tableViewsBefore: Array.from(grok.shell.tableViews).length,
        };
      });
      // initPeptides is idempotent; a non-null error is non-fatal (startAnalysis inits lazily).
      if (result.initError)
        console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', result.initError);
      expect(result.tableViewsBefore, 'shell should be clean after closeAll').toBe(0);
    });

    // Scenario 1 — Peptide SAR demo (dispatched via grok.functions.call; gallery card has no stable selector).
    await softStep('Scenario 1 (steps 1-3): activate Bioinformatics | Peptide SAR demo',
      async () => {
        const launched = await page.evaluate(async () => {
          const demoFunc = DG.Func.find({package: 'Peptides', name: 'macromoleculeSarFastaDemo'})[0];
          const demoPath = demoFunc?.options?.['demoPath'] ?? null;
          const isDemoDashboard = demoFunc?.options?.['isDemoDashboard'] ?? null;
          let demoError: string | null = null;
          try { await grok.functions.call('Peptides:macromoleculeSarFastaDemo'); }
          catch (e) { demoError = String(e); }
          return {
            demoFuncPresent: !!demoFunc,
            demoPath,
            isDemoDashboard,
            demoError,
          };
        });
        expect(launched.demoFuncPresent,
          'macromoleculeSarFastaDemo must be registered (demo-gallery discoverability)').toBe(true);
        expect(launched.demoPath,
          'demoPath must match the atlas-stated gallery position').toBe('Bioinformatics | Peptide SAR');
        expect(launched.isDemoDashboard,
          'isDemoDashboard must be "true" (renders as dashboard, not tutorial walkthrough)')
          .toBe('true');
        expect(launched.demoError,
          'Peptide SAR demo invocation must not throw').toBeNull();
      });

    // Poll for the full viewer set (model attaches before viewers, so gate on the viewers).
    await softStep('Scenario 1 (step 3 settle): wait for SAR dashboard viewers to attach',
      async () => {
        await page.waitForFunction(() => {
          const tv = Array.from(grok.shell.tableViews)
            .find((v) => v.dataFrame.temp['peptidesModel']);
          if (!tv) return false;
          const types = new Set(Array.from(tv.viewers).map((v) => v.type));
          return types.has('MCL') &&
            types.has('Logo Summary Table') &&
            types.has('Sequence Variability Map') &&
            types.has('Most Potent Residues');
        }, {timeout: 90_000});
        await page.waitForTimeout(2000);
      });

    await softStep('Scenario 1 (step 3 verify): "Simple peptides" TableView opened',
      async () => {
        const state = await page.evaluate(() => {
          const tv = Array.from(grok.shell.tableViews)
            .find((v) => v.dataFrame.temp['peptidesModel']);
          const df = tv?.dataFrame;
          const seq = df?.col('AlignedSequence');
          return {
            tvFound: !!tv,
            tableName: df?.name ?? null,
            rows: df?.rowCount ?? null,
            hasAlignedSeq: !!seq,
            semType: seq?.semType ?? null,
          };
        });
        expect(state.tvFound, 'SAR TableView (PeptidesModel-bearing) must exist').toBe(true);
        expect(state.tableName,
          'demo loads a DataFrame named "Simple peptides"').toBe('Simple peptides');
        expect(state.hasAlignedSeq,
          'AlignedSequence column must be present on the demo dataset').toBe(true);
        expect(state.semType,
          'AlignedSequence must be a Macromolecule column').toBe('Macromolecule');
      });

    await softStep('Scenario 1 (step 4): FASTA / aligned tags on AlignedSequence',
      async () => {
        const tags = await page.evaluate(() => {
          const tv = Array.from(grok.shell.tableViews)
            .find((v) => v.dataFrame.temp['peptidesModel']);
          const seq = tv?.dataFrame?.col('AlignedSequence');
          return {
            alphabet: seq?.getTag('alphabet') ?? null,
            units: seq?.meta?.units ?? null,
            aligned: seq?.getTag('aligned') ?? null,
          };
        });
        expect(tags.alphabet, 'alphabet tag must be PT').toBe('PT');
        expect(tags.units, 'units must be fasta').toBe('fasta');
        expect(tags.aligned, 'aligned tag must be SEQ.MSA').toBe('SEQ.MSA');
      });

    await softStep('Scenario 1 (step 5): activity scaling = -lg on the model', async () => {
      const scaling = await page.evaluate(() => {
        const tv = Array.from(grok.shell.tableViews)
          .find((v) => v.dataFrame.temp['peptidesModel']);
        const model = tv?.dataFrame?.temp['peptidesModel'] as any;
        return model?._settings?.activityScaling ?? null;
      });
      expect(scaling, 'model activityScaling must be the SCALING_METHODS.MINUS_LG canonical "-lg"')
        .toBe('-lg');
    });

    await softStep('Scenario 1 (step 6): MCL clustering on with threshold = 94',
      async () => {
        const mcl = await page.evaluate(() => {
          const tv = Array.from(grok.shell.tableViews)
            .find((v) => v.dataFrame.temp['peptidesModel']);
          const model = tv?.dataFrame?.temp['peptidesModel'] as any;
          const viewerTypes = tv ? Array.from(tv.viewers).map((v) => v.type) : [];
          return {
            threshold: model?._settings?.mclSettings?.threshold ?? null,
            hasMclViewer: viewerTypes.includes('MCL'),
            viewerTypes,
          };
        });
        expect(mcl.threshold, 'mclSettings.threshold must be 94').toBe(94);
        expect(mcl.hasMclViewer,
          'MCL viewer must attach when addMCL:true was passed to startAnalysis').toBe(true);
      });

    await softStep('Scenario 1 (step 7): SAR dashboard layout rendered, no fatal crash',
      async () => {
        const state = await page.evaluate(async () => {
          const tv = Array.from(grok.shell.tableViews)
            .find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
          const viewerTypes = Array.from(tv.viewers).map((v) => v.type);
          // grok.shell.lastError is a Promise — await it (raw String() yields "[object Promise]").
          let lastError = '';
          try { lastError = (await grok.shell.lastError) ?? ''; } catch (_) {}
          const hasNullCrash = /setTrue|null receiver|NoSuchMethodError|Cannot read prop|undefined is not/i
            .test(lastError);
          return {viewerTypes, lastError, hasNullCrash};
        });
        expect(state.viewerTypes,
          'Sequence Variability Map must attach on the demo dashboard').toContain('Sequence Variability Map');
        expect(state.viewerTypes,
          'Most Potent Residues must attach on the demo dashboard').toContain('Most Potent Residues');
        expect(state.viewerTypes,
          'Logo Summary Table must attach on the demo dashboard (MCL clustering is on)')
          .toContain('Logo Summary Table');
        expect(state.hasNullCrash,
          'no null-receiver / NoSuchMethodError crash during the demo load')
          .toBe(false);
      });

    // Scenario 2 — Peptides landing View invoked via grok.functions.call (no #app on this build).
    await softStep('Scenario 2 (steps 1-3): open Peptides landing view + verify shape',
      async () => {
        const result = await page.evaluate(async () => {
          grok.shell.closeAll();
          await new Promise((r) => setTimeout(r, 3000));

          const appsList = DG.Func.find({tags: ['app']})
            .map((f) => `${f.package?.name}:${f.name}`);
          const peptidesInAppsList = appsList.some((s) => /^Peptides:/i.test(s));

          const view = await grok.functions.call('Peptides:Peptides');
          grok.shell.addView(view);
          await new Promise((r) => setTimeout(r, 1500));

          const root = view?.root as HTMLElement | undefined;
          const buttonLabels = root
            ? Array.from(root.querySelectorAll('button'))
              .map((b) => (b as HTMLElement).textContent?.trim() ?? '')
            : [];

          return {
            peptidesInAppsList,
            viewName: view?.name ?? null,
            buttonLabels,
            // The three window-flag assertions (package.ts L107-L110).
            flags: {
              showToolbox: grok.shell.windows.showToolbox,
              showHelp: grok.shell.windows.showHelp,
              showProperties: grok.shell.windows.showProperties,
            },
          };
        });

        if (result.peptidesInAppsList) {
          console.log('[note] Peptides:Peptides IS in apps-browser app-list on this build.');
        } else {
          console.log('[note] Peptides:Peptides is NOT in apps-browser app-list on this build ' +
            '(scenario framing drift; spec drives the View-function path instead).');
        }
        expect(result.viewName, 'landing view must be named "Peptides"').toBe('Peptides');
        expect(result.buttonLabels, 'landing view must show exactly three demo buttons')
          .toEqual(['Simple demo', 'Complex demo', 'HELM demo']);
        expect(result.flags.showToolbox, 'showToolbox must be false on Peptides view').toBe(false);
        expect(result.flags.showHelp, 'showHelp must be false on Peptides view').toBe(false);
        expect(result.flags.showProperties, 'showProperties must be false on Peptides view').toBe(false);
      });

    await softStep('Scenario 2 (step 4): click "Simple demo" button (DOM-driving)',
      async () => {
        const tvsBefore = await page.evaluate(() => Array.from(grok.shell.tableViews).length);
        await page.locator('[name="button-Simple-demo"]').click();
        await page.waitForFunction((before) => Array.from(grok.shell.tableViews).length > before,
          tvsBefore, {timeout: 30_000});
        await page.waitForTimeout(1500);
        const newTv = await page.evaluate(() => {
          const tvs = Array.from(grok.shell.tableViews);
          const tv = tvs[tvs.length - 1];
          const df = tv?.dataFrame;
          const seq = df?.col('AlignedSequence');
          return {
            viewName: tv?.name ?? null,
            tableName: df?.name ?? null,
            rows: df?.rowCount ?? null,
            hasAlignedSeq: !!seq,
            macromoleculeColCount: df ? Array.from({length: df.columns.length}, (_, i) =>
              df.columns.byIndex(i)).filter((c) => c.semType === 'Macromolecule').length : 0,
          };
        });
        expect(newTv.viewName,
          'new TableView from Simple demo must be named "PeptidesView"').toBe('PeptidesView');
        expect(newTv.hasAlignedSeq,
          'Simple demo (aligned.csv) ships an AlignedSequence column').toBe(true);
        expect(newTv.macromoleculeColCount,
          'Simple demo dataset must carry >=1 Macromolecule-semType column').toBeGreaterThan(0);
      });

    // Complex demo (aligned_2.csv) ships a Macromolecule column named "MSA", not "AlignedSequence".
    await softStep('Scenario 2 (step 5): return to landing view + click "Complex demo"',
      async () => {
        await page.evaluate(() => {
          const v = Array.from(grok.shell.views).find((view) => view.name === 'Peptides');
          if (v) grok.shell.v = v;
        });
        await page.waitForTimeout(800);
        const tvsBefore = await page.evaluate(() => Array.from(grok.shell.tableViews).length);
        await page.locator('[name="button-Complex-demo"]').click();
        await page.waitForFunction((before) => Array.from(grok.shell.tableViews).length > before,
          tvsBefore, {timeout: 30_000});
        await page.waitForTimeout(1500);
        const newTv = await page.evaluate(() => {
          const tvs = Array.from(grok.shell.tableViews);
          const tv = tvs[tvs.length - 1];
          const df = tv?.dataFrame;
          const macromoleculeCols = df ? Array.from({length: df.columns.length}, (_, i) =>
            df.columns.byIndex(i)).filter((c) => c.semType === 'Macromolecule').map((c) => c.name) : [];
          return {
            viewName: tv?.name ?? null,
            tableName: df?.name ?? null,
            rows: df?.rowCount ?? null,
            macromoleculeCols,
          };
        });
        expect(newTv.viewName,
          'new TableView from Complex demo must be named "PeptidesView"').toBe('PeptidesView');
        expect(newTv.macromoleculeCols.length,
          'Complex demo dataset must carry >=1 Macromolecule-semType column ' +
          '(name is "MSA" on this build, not "AlignedSequence" — scenario column-name drift)')
          .toBeGreaterThan(0);
      });

    // HELM demo (aligned_3.csv) ships a Macromolecule column named "HELM" (units:"helm").
    await softStep('Scenario 2 (step 6): return to landing view + click "HELM demo"',
      async () => {
        await page.evaluate(() => {
          const v = Array.from(grok.shell.views).find((view) => view.name === 'Peptides');
          if (v) grok.shell.v = v;
        });
        await page.waitForTimeout(800);
        const tvsBefore = await page.evaluate(() => Array.from(grok.shell.tableViews).length);
        await page.locator('[name="button-HELM-demo"]').click();
        await page.waitForFunction((before) => Array.from(grok.shell.tableViews).length > before,
          tvsBefore, {timeout: 30_000});
        await page.waitForTimeout(1500);
        const newTv = await page.evaluate(() => {
          const tvs = Array.from(grok.shell.tableViews);
          const tv = tvs[tvs.length - 1];
          const df = tv?.dataFrame;
          const macromoleculeCols = df ? Array.from({length: df.columns.length}, (_, i) =>
            df.columns.byIndex(i)).filter((c) => c.semType === 'Macromolecule').map((c) => ({
              name: c.name, units: c.meta?.units ?? null,
            })) : [];
          return {
            viewName: tv?.name ?? null,
            tableName: df?.name ?? null,
            rows: df?.rowCount ?? null,
            macromoleculeCols,
          };
        });
        expect(newTv.viewName,
          'new TableView from HELM demo must be named "PeptidesView"').toBe('PeptidesView');
        expect(newTv.macromoleculeCols.length,
          'HELM demo dataset must carry >=1 Macromolecule-semType column ' +
          '(name is "HELM" with units:"helm" on this build, not "AlignedSequence" — ' +
          'scenario column-name drift)')
          .toBeGreaterThan(0);
        if (newTv.macromoleculeCols.some((c) => c.units === 'helm'))
          console.log('[note] HELM demo carries a Macromolecule column with units:"helm" ' +
            '(scenario\'s HELM-notation distinction holds).');
      });

    await softStep('Scenario 2 (step 7): no null-receiver / fatal crash across button clicks',
      async () => {
        const errState = await page.evaluate(async () => {
          let lastError = '';
          try { lastError = (await grok.shell.lastError) ?? ''; } catch (_) {}
          const hasNullCrash = /setTrue|null receiver|NoSuchMethodError|Cannot read prop|undefined is not/i
            .test(lastError);
          return {lastError, hasNullCrash};
        });
        expect(errState.hasNullCrash,
          'no null-receiver / NoSuchMethodError crash during the three demo button clicks')
          .toBe(false);
      });

    finishSpec();
  });
