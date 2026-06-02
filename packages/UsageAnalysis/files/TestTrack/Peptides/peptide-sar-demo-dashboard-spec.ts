/* ---
sub_features_covered: [peptides.demos.macromolecule-sar-fasta, peptides.app, peptides.lifecycle.init, peptides.workflow.start-analysis, peptides.model.add-sequence-space]
--- */
// target_layer: playwright; coverage_type: smoke. No pyramid_layer, so JS API
// substitution is permitted, but >=1 DOM-driving call is still required per
// E-LAYER-COMPLIANCE-01.
//   sub_features_covered: [peptides.demos.macromolecule-sar-fasta, peptides.app,
//                          peptides.lifecycle.init, peptides.workflow.start-analysis,
//                          peptides.model.add-sequence-space]
//
// Atlas provenance (derived_from):
//   feature-atlas/peptides.yaml#critical_paths[peptide-sar-demo-dashboard-from-gallery]
//     derived_from: public/packages/Peptides/src/package.ts#L270
//   feature-atlas/peptides.yaml#critical_paths[open-peptides-app-and-load-demo]
//     derived_from: public/packages/Peptides/src/package.ts#L93
//   feature-atlas/peptides.yaml#interactions[peptides-demo-dashboard-loads-from-gallery]
//
// Peptide SAR demo dashboard + Peptides app landing — entry-point smokes.
// Scenario 1 covers peptides.demos.macromolecule-sar-fasta (FASTA dashboard demo)
// and Scenario 2 covers peptides.app (three-button app landing view). Both retire
// coverage-gaps/peptides.yaml gap[1].proposal sub-clause (f).
//
// Empirical recon notes (live MCP recon 2026-05-29 on dev.datagrok.ai — these
// drive the deterministic assertions, contradicting two scenario claims that do
// NOT hold on this build):
//
//   - SCENARIO 2 FRAMING DRIFT (atlas-vs-build mismatch): the scenario describes
//     the Peptides app as "discoverable in the app browser", citing atlas
//     critical_path open-peptides-app-and-load-demo derived_from package.ts#L93.
//     Live recon: the Peptides function at package.ts#L93 is registered as a
//     plain func returning DG.View — NOT as an @app (no tags:['app'], no
//     meta.role:'app'). DG.Func.find({tags:['app']}).filter(f => f.package.name
//     === 'Peptides') == [] on this build (peptides v1.27.9). The
//     grok-browser/references/peptides.md "Peptide SAR demo & app entry" §
//     documents this finding explicitly ("No standalone Peptides app — Peptides
//     registers no #app"). The scenario's broader contract (the View function
//     EXISTS, returns the landing layout with three demo buttons + window-flag
//     toggling) DOES hold and IS driven here via grok.functions.call
//     ('Peptides:Peptides') + grok.shell.addView — the sanctioned API path that
//     bypasses the missing apps-browser surface. The View's three buttons are
//     then exercised through REAL DOM clicks (Locator API) for the smoke
//     contract (steps 4/5/6 of Scenario 2 — the only deterministic DOM gestures
//     on the landing view).
//
//   - SCENARIO 2 COLUMN-NAME DRIFT (steps 5 & 6): the scenario asserts presence
//     of an "AlignedSequence" column as the deterministic shape check for all
//     three demos. Live recon (2026-05-29) shows that this only holds for the
//     Simple demo (aligned.csv: AlignedSequence/Macromolecule/separator). The
//     Complex demo (aligned_2.csv) ships a Macromolecule column named "MSA"
//     (units:"separator"), and the HELM demo (aligned_3.csv) ships one named
//     "HELM" (units:"helm"). The semantically correct shape check is presence
//     of A Macromolecule-semType column, not the literal name "AlignedSequence"
//     — this spec asserts the semType-by-name-agnostic contract for Complex
//     and HELM (and asserts AlignedSequence specifically for Simple where it
//     does hold), preserving the smoke-class intent without papering over the
//     scenario's hand-coded name.
//
//   - SCENARIO 1 DEMO INVOCATION: the scenario describes navigating the demo
//     gallery (top-menu Help|Demo) and clicking the "Peptide SAR" entry under
//     "Bioinformatics". The demo function is registered with
//     options.demoPath:"Bioinformatics | Peptide SAR", isDemoDashboard:"true"
//     (verified live via DG.Func.find), and the canonical user-equivalent
//     invocation is grok.functions.call('Peptides:macromoleculeSarFastaDemo')
//     — the same call the gallery's dashboard-launcher button dispatches on
//     activation. The demo function returns Promise<void> and renders the SAR
//     dashboard end-to-end (Simple peptides TableView, FASTA/PT tags, -lg
//     scaling, MCL threshold 94, default attach set including Logo Summary
//     Table because MCL is on). Live recon confirms every assertion in
//     Scenario 1's Expected block: table name "Simple peptides", 647 rows,
//     AlignedSequence Macromolecule with alphabet:"PT", units:"fasta",
//     aligned:"SEQ.MSA", model._settings.activityScaling === "-lg",
//     model._settings.mclSettings.threshold === 94, viewers include Sequence
//     Variability Map / Most Potent Residues / MCL / Logo Summary Table.
//
//   - GROK-17557 init-prerequisite tolerance: Setup pre-warms the Peptides
//     @init (Peptides:initPeptides) so MonomerWorks + TreeHelper +
//     PeptideUtils.loadComponents() are resolved before the demo's startAnalysis
//     leg requests them. Removes the cold-package init race (same pattern as
//     peptide-space-spec.ts Setup). The demo function tolerates the warm path;
//     this just makes the spec's first-launch behavior deterministic.
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser
// reference). Only DOM selectors recorded here — JS-API surface use does not
// require a recon-note:
//   [name="button-Simple-demo"] / [name="button-Complex-demo"] /
//     [name="button-HELM-demo"] — the three ui.button(...) widgets on the
//     Peptides landing view (package.ts L117-L119). ui.button() DOES annotate
//     with name="button-<label>" (verified live 2026-05-29 via
//     chrome-devtools MCP evaluate_script — outerHTML shows
//     <button class="ui-btn ui-btn-ok" name="button-Simple-demo" role="button">
//     <span>Simple demo</span></button>). Not in peptides.md (its "Peptide SAR
//     demo & app entry" § documents the landing view exists but does not
//     enumerate per-button selectors). The earlier hand-coded conclusion in
//     this file that the buttons carry no name= was wrong; the canonical
//     Datagrok [name="button-<Label>"] pattern applies.
//
// Technique notes:
//
//   - grok.shell.tableViews is Iterable<TableView>, NOT Array (shell.ts L333)
//     — `.length` is undefined on it; wrap with Array.from(...) before reading
//     length. MCP-recon backed.
//
//   - The three demo-button clicks use page.locator('[name="button-Simple-demo"]')
//     (and Complex / HELM) rather than getByRole accessible-name. The
//     [name="button-X-demo"] attribute is set by ui.button via the platform's
//     annotate(); an attribute selector is more deterministic than accessible-name
//     computation, which on cold init can race with the accessibility-tree build.
//
//   - Scenario 1 readiness polls up to 90s for ALL of MCL + Logo Summary Table +
//     Sequence Variability Map + Most Potent Residues to attach. PeptidesModel
//     attaches at the START of startAnalysis (PeptidesModel.getInstance(df),
//     Peptides/src/model.ts L368-373), BEFORE addMonomerPosition /
//     addMostPotentResidues / addMCLClusters / addLogoSummaryTable run — so a
//     PeptidesModel-attach gate would fire well before the viewers materialize.
//     On a cold worker MCL clustering takes 30-60s; polling on the actual
//     viewer-set the assertions read is the deterministic readiness signal.
//
//   - grok.shell.lastError is Promise<string|undefined> (js-api shell.ts L107);
//     it must be awaited. Stringifying the Promise directly yields "[object
//     Promise]", which never matches the null-receiver regex (a silent no-op).
//     Awaiting turns it into a real crash signal without changing pass-fail
//     semantics on the warm path (lastError resolves to undefined or a benign
//     string).
//
//   - The Scenario-1 -> Scenario-2 transition waits 3000ms after closeAll(). On
//     cold workers the SAR dashboard cleanup (multiple viewers, MCL worker
//     shutdown) can take >800ms; 3s matches the conservative cleanup wait used in
//     sibling info-panels-spec.ts.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Peptide SAR demo dashboard + Peptides app landing — entry-point smokes',
  async ({page}) => {
    // Demo invocation runs startAnalysis end-to-end (MCL clustering + sequence
    // space embedding in Web Workers + Logo Summary Table population); plus the
    // app-landing-view flow opens three TableViews. Generous budget.
    test.setTimeout(360_000);

    await loginToDatagrok(page);

    // ---- Setup ----
    // Clean shell (no prior PeptidesModel / TableView), Windows mode (the SAR
    // dashboard docks several viewers — same mode as peptide-space-spec.ts),
    // and pre-warm Peptides:initPeptides so MonomerWorks/TreeHelper/PeptideUtils
    // singletons are resolved before the demo's startAnalysis leg requests them
    // (GROK-17557 family — cold init race).
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
        // Pre-warm @init — non-fatal if it has already run (cached singletons).
        let initError: string | null = null;
        try { await grok.functions.call('Peptides:initPeptides'); }
        catch (e) { initError = String(e); }
        return {
          initError,
          // grok.shell.tableViews is Iterable<TableView>, NOT Array — .length is undefined on it
          // Array.from(...) before reading length.
          tableViewsBefore: Array.from(grok.shell.tableViews).length,
        };
      });
      // initPeptides is idempotent; a non-null error here is non-fatal (the
      // demo's startAnalysis will trigger init lazily). Surface softly.
      if (result.initError)
        console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', result.initError);
      expect(result.tableViewsBefore, 'shell should be clean after closeAll').toBe(0);
    });

    // ====================================================================
    // Scenario 1 — Peptide SAR demo from the demo gallery loads the FASTA
    // dataset, applies -lg scaling and MCL clustering, and renders the SAR
    // dashboard.
    //
    // Atlas anchor: critical_paths[peptide-sar-demo-dashboard-from-gallery]
    //   + interactions[peptides-demo-dashboard-loads-from-gallery].
    //
    // Demo invocation: grok.functions.call('Peptides:macromoleculeSarFastaDemo')
    // — the same function the demo-gallery's "Peptide SAR" entry dispatches on
    // activation (registered with meta.demoPath="Bioinformatics | Peptide SAR",
    // isDemoDashboard="true"). The gallery's per-entry DOM card has no stable
    // platform-level [name=] anchor; the function-call surface IS the user
    // gesture for a demo-dashboard demo (the gallery is the discovery surface,
    // the function call is the dispatch surface).
    // ====================================================================

    // Steps 1-3: Open the gallery, navigate to Bioinformatics | Peptide SAR, and
    // activate the entry. Driven through grok.functions.call for the dispatch
    // (the gallery card has no stable selector); the demo function's
    // demoPath registration is the discoverability contract.
    await softStep('Scenario 1 (steps 1-3): activate Bioinformatics | Peptide SAR demo',
      async () => {
        const launched = await page.evaluate(async () => {
          // Discoverability contract: the demo function is registered with the
          // expected demoPath + isDemoDashboard flag. Read live from the
          // function registry (DG.Func.find with options.demoPath check).
          const demoFunc = DG.Func.find({package: 'Peptides', name: 'macromoleculeSarFastaDemo'})[0];
          const demoPath = demoFunc?.options?.['demoPath'] ?? null;
          const isDemoDashboard = demoFunc?.options?.['isDemoDashboard'] ?? null;
          // Dispatch (the gallery click equivalent).
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

    // Wait for the SAR analysis to fully settle. startAnalysis runs MCL
    // clustering + sequence-space embedding in Web Workers; PeptidesModel
    // attaches at the START of startAnalysis (Peptides/src/model.ts L368-373
    // — getInstance(df) sets df.temp.peptidesModel before any addViewer call),
    // so a PeptidesModel-attach gate fires WAY before the dashboard is
    // populated. The deterministic readiness signal is the full viewer set
    // the assertions below read: MCL + Logo Summary Table + Sequence
    // Variability Map + Most Potent Residues. Poll up to 90s for all four.
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
        // Brief settle for any final paint of the just-attached viewers.
        await page.waitForTimeout(2000);
      });

    // Step 3 (assertion leg): TableView named "Simple peptides" opens with an
    // AlignedSequence Macromolecule column present.
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

    // Step 4: FASTA / aligned tags applied — alphabet=PT, units=fasta,
    // aligned=SEQ.MSA. These are the demo's per-column setup (fasta.ts L37-L41).
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

    // Step 5: -lg scaling applied — model._settings.activityScaling === "-lg"
    // (the demo calls scaleActivity(activityCol, SCALING_METHODS.MINUS_LG) and
    // passes the same constant to startAnalysis as scalingMethod).
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

    // Step 6: MCL clustering enabled and mclSettings.threshold = 94. The demo
    // passes {addMCL: true, useEmbeddingsClusters: true, mclSettings: {threshold:94}}
    // (fasta.ts L43-L47). Verify via the model setting + the MCL viewer attaching.
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

    // Step 7: SAR dashboard layout rendered — the meta.isDemoDashboard:"true"
    // contract requires the full SAR layout to materialize. Assert the
    // deterministic-attach set: Grid + Sequence Variability Map + Most Potent
    // Residues + MCL + Logo Summary Table (LST attaches here because MCL is
    // on, per peptides.md "Logo Summary Table is cluster-dependent"). Also
    // assert no null-receiver / fatal crash signal (the only acceptable
    // grok.shell.lastError signal observed live is "[object Promise]", which
    // is innocuous wrapping by the platform — same as peptide-space-spec.ts).
    await softStep('Scenario 1 (step 7): SAR dashboard layout rendered, no fatal crash',
      async () => {
        const state = await page.evaluate(async () => {
          const tv = Array.from(grok.shell.tableViews)
            .find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
          const viewerTypes = Array.from(tv.viewers).map((v) => v.type);
          // grok.shell.lastError is Promise<string|undefined> (js-api shell.ts L107)
          // — must await it; the earlier `String(grok.shell.lastError)` returned
          // "[object Promise]" unconditionally (the null-receiver regex never
          // matched). On the warm path lastError resolves to undefined or a
          // benign Dart-side message; on a real null-receiver crash it carries
          // the stringified Dart error.
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

    // ====================================================================
    // Scenario 2 — Peptides app landing view opens with three demo buttons,
    // and each Simple / Complex / HELM button loads its sample into a new
    // PeptidesView TableView.
    //
    // Atlas anchor: critical_paths[open-peptides-app-and-load-demo].
    //
    // BUILD-DRIFT NOTE (see "SCENARIO 2 FRAMING DRIFT" in the header): the
    // scenario describes app-browser discovery, but the Peptides function at
    // package.ts#L93 lacks the @app tag (DG.Func.find({tags:['app']}) does
    // NOT list Peptides on this build). The function nonetheless EXISTS and
    // returns the landing layout; sanctioned API path is
    // grok.functions.call('Peptides:Peptides') + grok.shell.addView, then the
    // three buttons are exercised via REAL DOM clicks (Locator API) — the
    // smoke-class contract.
    // ====================================================================

    // Steps 1-3 of Scenario 2 (combined): close prior SAR view, invoke the
    // Peptides View function, verify landing-view shape (app header + 3
    // labeled buttons), and confirm the three window flags toggled to false.
    // The Locator-API click on the buttons (Step 4 below) is the load-bearing
    // DOM-driving call for this spec (satisfies E-LAYER-COMPLIANCE-01).
    await softStep('Scenario 2 (steps 1-3): open Peptides landing view + verify shape',
      async () => {
        const result = await page.evaluate(async () => {
          // Tear down the Scenario 1 SAR view so Scenario 2 starts clean.
          // On cold workers the SAR cleanup (multiple viewers + MCL worker
          // shutdown) can take >800ms; 3s matches the conservative cleanup
          // wait used in sibling info-panels-spec.ts.
          grok.shell.closeAll();
          await new Promise((r) => setTimeout(r, 3000));

          // BUILD-DRIFT CHECK: verify the documented finding that Peptides is
          // NOT in the apps-browser app-list (scenario says it should be — it
          // is not). Recorded as a soft signal, not a fail (the scenario's
          // broader contract about the View function holds).
          const appsList = DG.Func.find({tags: ['app']})
            .map((f) => `${f.package?.name}:${f.name}`);
          const peptidesInAppsList = appsList.some((s) => /^Peptides:/i.test(s));

          // Sanctioned API path: invoke the registered View function directly
          // (no apps-browser DOM affordance for it on this build).
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

        // Build-drift surfaced (soft): Peptides not in apps list. Not a fail
        // — the scenario's broader contract about landing-view shape holds.
        if (result.peptidesInAppsList) {
          // If the package starts registering @app in the future, this
          // becomes truthy — the scenario's framing then aligns with the
          // build. No-op note.
          console.log('[note] Peptides:Peptides IS in apps-browser app-list on this build.');
        } else {
          console.log('[note] Peptides:Peptides is NOT in apps-browser app-list on this build ' +
            '(scenario framing drift; spec drives the View-function path instead).');
        }
        expect(result.viewName, 'landing view must be named "Peptides"').toBe('Peptides');
        // Three buttons exactly: Simple demo, Complex demo, HELM demo —
        // no more, no fewer (scenario Expected block).
        expect(result.buttonLabels, 'landing view must show exactly three demo buttons')
          .toEqual(['Simple demo', 'Complex demo', 'HELM demo']);
        // Window flags (Scenario 2 step 3 expected): all three false.
        expect(result.flags.showToolbox, 'showToolbox must be false on Peptides view').toBe(false);
        expect(result.flags.showHelp, 'showHelp must be false on Peptides view').toBe(false);
        expect(result.flags.showProperties, 'showProperties must be false on Peptides view').toBe(false);
      });

    // Step 4: click the "Simple demo" button — REAL DOM gesture via Playwright
    // Locator API. The handler is openDemoData('aligned.csv'); result is a new
    // TableView named "PeptidesView" containing the Simple sample (which DOES
    // ship an AlignedSequence Macromolecule column).
    await softStep('Scenario 2 (step 4): click "Simple demo" button (DOM-driving)',
      async () => {
        // grok.shell.tableViews is Iterable<TableView>, not Array — wrap with
        // Array.from before reading .length.
        const tvsBefore = await page.evaluate(() => Array.from(grok.shell.tableViews).length);
        // Locator-API click on the button — the load-bearing DOM-driving call.
        // [name="button-Simple-demo"] is set by ui.button() via annotate(); the
        // attribute selector is more deterministic than getByRole accessible-name
        // on cold init.
        await page.locator('[name="button-Simple-demo"]').click();
        // openDemoData -> grok.data.loadTable -> addTableView (async).
        await page.waitForFunction((before) => Array.from(grok.shell.tableViews).length > before,
          tvsBefore, {timeout: 30_000});
        await page.waitForTimeout(1500);
        const newTv = await page.evaluate(() => {
          const tvs = Array.from(grok.shell.tableViews);
          const tv = tvs[tvs.length - 1];
          const df = tv?.dataFrame;
          const seq = df?.col('AlignedSequence');
          // For Simple, AlignedSequence IS the canonical column name.
          return {
            viewName: tv?.name ?? null,
            tableName: df?.name ?? null,
            rows: df?.rowCount ?? null,
            hasAlignedSeq: !!seq,
            // Build-drift-safe alternative: a Macromolecule-semType column
            // exists (semantically correct shape check regardless of name).
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

    // Step 5: return to Peptides app view, click "Complex demo" (Locator API
    // click). BUILD-DRIFT (per header note): the Complex demo dataset
    // (aligned_2.csv) carries a Macromolecule column named "MSA" — NOT
    // "AlignedSequence". The semantically correct shape check is Macromolecule
    // presence, not the literal name.
    await softStep('Scenario 2 (step 5): return to landing view + click "Complex demo"',
      async () => {
        // Switch back to the Peptides landing view.
        await page.evaluate(() => {
          const v = Array.from(grok.shell.views).find((view) => view.name === 'Peptides');
          if (v) grok.shell.v = v;
        });
        await page.waitForTimeout(800);
        // grok.shell.tableViews is Iterable<TableView>, not Array — wrap with
        // Array.from before reading .length.
        const tvsBefore = await page.evaluate(() => Array.from(grok.shell.tableViews).length);
        // Locator-API click on the [name="button-Complex-demo"] attribute selector.
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
            macromoleculeCols,  // empirical column name; "MSA" on this build
          };
        });
        expect(newTv.viewName,
          'new TableView from Complex demo must be named "PeptidesView"').toBe('PeptidesView');
        // Build-drift-tolerant assertion (NOT the scenario's literal
        // "AlignedSequence" name — Complex demo ships "MSA").
        expect(newTv.macromoleculeCols.length,
          'Complex demo dataset must carry >=1 Macromolecule-semType column ' +
          '(name is "MSA" on this build, not "AlignedSequence" — scenario column-name drift)')
          .toBeGreaterThan(0);
      });

    // Step 6: return to Peptides app view, click "HELM demo" (Locator API
    // click). BUILD-DRIFT (per header note): the HELM demo dataset
    // (aligned_3.csv) carries a Macromolecule column named "HELM" (with
    // units:"helm"). Same semantically correct shape check as Complex.
    await softStep('Scenario 2 (step 6): return to landing view + click "HELM demo"',
      async () => {
        await page.evaluate(() => {
          const v = Array.from(grok.shell.views).find((view) => view.name === 'Peptides');
          if (v) grok.shell.v = v;
        });
        await page.waitForTimeout(800);
        // grok.shell.tableViews is Iterable<TableView>, not Array — wrap with
        // Array.from before reading .length.
        const tvsBefore = await page.evaluate(() => Array.from(grok.shell.tableViews).length);
        // Locator-API click on the [name="button-HELM-demo"] attribute selector.
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
        // Soft note: the HELM-notation distinction is observable through the
        // units tag of the macromolecule column (per scenario step 6 hint
        // about HELM-notation peptides vs FASTA).
        if (newTv.macromoleculeCols.some((c) => c.units === 'helm'))
          console.log('[note] HELM demo carries a Macromolecule column with units:"helm" ' +
            '(scenario\'s HELM-notation distinction holds).');
      });

    // Step 7: no null-receiver / load-failure console errors across the three
    // demo-button clicks. Same fatal-crash detection as Scenario 1.
    await softStep('Scenario 2 (step 7): no null-receiver / fatal crash across button clicks',
      async () => {
        const errState = await page.evaluate(async () => {
          // Same lastError treatment as Scenario 1 step 7 — properly await
          // the Promise<string|undefined> instead of stringifying it raw.
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

    // Final: surface any softStep failure loudly so the run fails on a real
    // assertion miss (rather than silently passing on the page-level expect).
    if (stepErrors.length > 0) {
      throw new Error('soft-step failures:\n' +
        stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
    }
  });
