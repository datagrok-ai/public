---
feature: peptides
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [open-peptides-app-and-load-demo, peptide-sar-demo-dashboard-from-gallery]
realizes: [bio.analyze.sar, bio.analyze.sequence-space]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - peptide-sar-demo-dashboard-spec.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-31T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-31T01:10:00Z
    spec_runs:
      - spec: peptide-sar-demo-dashboard-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 117
        failure_keys: []
---

# Peptides — Demo gallery dashboard and app-browser entry points

Checks the two "getting started" entry points for Peptides: the **Bioinformatics | Peptide SAR** demo in the platform's demo gallery, which loads a FASTA sample and renders a full SAR dashboard; and the **Peptides** app itself, whose landing page offers Simple/Complex/HELM demo buttons that each open a sample dataset. Both paths are expected to work independently of any prior state, with no errors on first load.

## Setup

1. Start from a clean Datagrok shell — no Peptides analysis or PeptidesModel on any open DataFrame, no pre-existing `PeptidesView` TableView. The demo and app entry paths exercised by this scenario are entry-point smokes; they assume nothing pre-loaded and must succeed independently of prior state. Closing any pre-existing peptide-related TableViews and clearing `grok.shell.tableViews` from prior test runs is the operative precondition.
2. Confirm the Peptides package is loaded and its `@init` `initPeptides` function has completed at least once during the session (the package's lifecycle init wires `MonomerWorks`, `TreeHelper`, and `PeptideUtils.loadComponents()` — the demo's `startAnalysis` call and the app's `openDemoData(...)` button handlers both depend on these singletons being resolvable). When invoking the package for the first time in a session, allow up to a few seconds after the first `Bio | Analyze | SAR...` or app open for the init handler to complete; subsequent invocations short-circuit on the cached singletons (`monomerWorks ??=`, `treeHelper ??=` per `package.ts#L82`).

## Scenarios

### Scenario 1 — The Bioinformatics | Peptide SAR demo loads a FASTA dataset with -lg scaling and MCL clustering, and renders as a dashboard

Opens the demo gallery, navigates to Bioinformatics | Peptide SAR, and confirms the FASTA sample loads with -lg activity scaling and MCL clustering (threshold 94) applied automatically, rendering the full SAR layout as a demo dashboard rather than a step-by-step walkthrough.

1. Open the Datagrok demo gallery (top-menu `Help | Demo` or the platform's "Demo" affordance in the sidebar — the entry-point name surfaced by the platform is `Demo`, listing all functions registered with `meta.demoPath`). Confirm the gallery view opens listing the available demo categories.
2. Navigate the gallery to the `Bioinformatics` category and locate the `Peptide SAR` demo entry — this is the demo registered at `package.ts#L262` via `@grok.decorators.func({meta: {demoPath: 'Bioinformatics | Peptide SAR', isDemoDashboard: 'true'}, name: 'Peptide SAR', description: 'Peptide SAR Analysis demo on peptide sequences in FASTA format'})`. The `isDemoDashboard: 'true'` flag marks the demo as a dashboard-style demo (renders as a layout, not as a step-by-step tutorial walkthrough).
3. Activate the `Peptide SAR` demo entry — this invokes the registered `macromoleculeSarFastaDemo` function, which in turn calls `macromoleculeSarFastaDemoUI()` (`src/demo/fasta.ts#L12`). The function sets `grok.shell.windows.showContextPanel = true`, then reads the package's sample `aligned.csv` via `_package.files.readAsText('aligned.csv')` and constructs a DataFrame named `Simple peptides`. Confirm a new TableView opens with the FASTA peptides loaded — the active table's name is `Simple peptides` and an `AlignedSequence` column is present.
4. Confirm the FASTA / aligned tags applied to the sequence column — the demo sets `simpleAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE`, `simpleAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.PT)`, `simpleAlignedSeqCol.meta.units = NOTATION.FASTA`, and `simpleAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA)` (per `fasta.ts#L37-L41`). Read the `AlignedSequence` column tags back via the JS API and confirm `semType === 'Macromolecule'`, the alphabet tag equals `PT`, the units equal `fasta`, and the aligned tag equals `SEQ.MSA`. The demo's FASTA-vs-HELM distinction matters because the SAR pipeline's `SeqHandler.splitter()` resolution branches on `meta.units` — the demo intentionally exercises the FASTA path.
5. Confirm -lg scaling applied — the demo calls `scaleActivity(simpleActivityCol, C.SCALING_METHODS.MINUS_LG)` against the `IC50` column (`fasta.ts#L42`). After `startAnalysis` completes, the SAR setup carries a scaled activity column on the DataFrame (the model wires this as the working activity series for stats). Read the model's active scaling via `PeptidesModel.getInstance(table)?.settings?.activityScaling` and confirm it equals `'-lg'` (the canonical `SCALING_METHODS.MINUS_LG` value).
6. Confirm MCL clustering enabled and `mclSettings.threshold = 94` — `startAnalysis` is invoked with `{addMCL: true, useEmbeddingsClusters: true, mclSettings: mclSettings}` where `mclSettings.threshold = 94` (`fasta.ts#L43-L47`). Confirm the SAR layout includes the MCL viewer (rendered as a scatter plot of cluster size vs max activity / sequence-space scatter — `peptides.model.add-sequence-space` populates this) and the cluster column is present on the DataFrame (atlas critical_path entry explicitly names `peptides.model.add-sequence-space` as a covered sub_feature). Read `PeptidesModel.getInstance(table)?.settings?.mclSettings?.threshold` and confirm it equals `94`.
7. Confirm SAR dashboard layout rendered — the dashboard demo registers with `isDemoDashboard: 'true'`, so the entire SAR layout should appear: the peptides grid with per-position WebLogo column headers, the Sequence Variability Map (`MonomerPosition`) and/or Most Potent Residues viewer (per the demo's default viewer toggles), the Logo Summary Table viewer (cluster summary; populated because MCL is on), and the sequence-space scatter / MCL viewer. No null-receiver or `splitter()` console errors must surface during the load.

Expected (assertion summary):
- The `Peptide SAR` demo entry is discoverable under `Bioinformatics` in the demo gallery (registered via `meta.demoPath: 'Bioinformatics | Peptide SAR'`, `meta.isDemoDashboard: 'true'`).
- Activating the demo opens a TableView named `Simple peptides` with an `AlignedSequence` Macromolecule column tagged `alphabet: PT`, `units: fasta`, `aligned: SEQ.MSA`.
- The SAR model is created with `activityScaling: '-lg'` and `mclSettings.threshold: 94`.
- The SAR layout includes the Logo Summary Table viewer and the sequence-space / MCL viewer (MCL clustering ran on the dataset).
- No null-receiver console errors during the demo load (the dashboard renders cleanly end-to-end).

### Scenario 2 — The Peptides app's landing page opens with three demo buttons, each loading its sample into a new table view

Opens the Peptides app from the app browser and confirms its landing page shows exactly three buttons — Simple demo, Complex demo, HELM demo — each of which loads its corresponding sample CSV into a new table view named `PeptidesView`. These buttons are the only interactive surface on the landing page, so exercising all three confirms the app-browser entry path works across all three sample shapes.

1. Open the Datagrok app browser (the platform surfaces this via `Apps` in the sidebar or `Help | Apps`). The app browser lists all functions registered with role `app` — `Peptides` should appear as one of the registered apps (registered at `package.ts#L93` via `@grok.decorators.func()` returning `DG.View` with `view.name = 'Peptides'`).
2. Activate the `Peptides` app entry — this invokes the registered `Peptides` function which constructs a new `View` named `Peptides`. Confirm the new app view opens with the Peptides landing layout: an app header rendered via `u2.appHeader(...)` (containing the package icon, the `learnMoreUrl` link, and the marketing-copy description listing the SAR capabilities), and a horizontal button row (`ui.divH([...])`) with exactly three buttons labeled `Simple demo`, `Complex demo`, and `HELM demo`.
3. Confirm the shell window flags toggled — the app function sets `grok.shell.windows.showToolbox = false`, `showHelp = false`, `showProperties = false` (per `package.ts#L107-L110`); read these flags back via the JS API and confirm all three are `false`. These flags collapse the side panels so the landing view occupies the full canvas — the assertion guards against a regression that leaves the side panels open on the app entry.
4. Click the `Simple demo` button — the handler is `() => openDemoData('aligned.csv')` (`package.ts#L117`). The `openDemoData` helper loads the named CSV from the package's `files/` folder via `grok.data.loadTable(...)` and opens it as a new TableView; the resulting view's name is `PeptidesView` per the atlas sub_feature `peptides.app` description ("each calls `openDemoData(<csv>)` which loads a sample CSV via `grok.data.loadTable` into a new TableView named `PeptidesView`"). Confirm a new TableView appears, the active table contains the `aligned.csv` rows (presence of an `AlignedSequence` column is the deterministic shape check), and the view's `name` equals `PeptidesView`.
5. Return to the Peptides app view (re-open the app or switch back to the `Peptides` view in the views list), click the `Complex demo` button — the handler is `() => openDemoData('aligned_2.csv')` (`package.ts#L118`). Confirm a new TableView opens with the `aligned_2.csv` rows loaded; the resulting view's `name` is `PeptidesView`.
6. Return to the Peptides app view again, click the `HELM demo` button — the handler is `() => openDemoData('aligned_3.csv')` (`package.ts#L119`). Confirm a new TableView opens with the `aligned_3.csv` rows loaded; the resulting view's `name` is `PeptidesView`. The HELM sample's `AlignedSequence` column typically carries HELM-notation peptides (vs the FASTA notation of `aligned.csv`); presence of the `AlignedSequence` column is the deterministic shape check (semtype + units may vary between the three samples per the source CSVs).
7. Across all three demo button clicks, confirm no null-receiver or load-failure console errors fire (the package init must have completed at least once during the session per Setup step 2; otherwise the first button click is the de-facto init trigger and may show a brief delay while `initPeptides` resolves the singletons — but it must not error).

Expected:
- The `Peptides` app entry is discoverable in the app browser.
- The app's landing view shows the appHeader plus three labeled buttons: `Simple demo`, `Complex demo`, `HELM demo` — no more, no fewer.
- The shell window flags `showToolbox`, `showHelp`, `showProperties` are all `false` after the app opens.
- Each of the three buttons opens a new TableView named `PeptidesView` carrying the corresponding sample (`aligned.csv`, `aligned_2.csv`, `aligned_3.csv`).
- No null-receiver console errors across the three button clicks.

## Notes

- **Why bundled together.** Both scenarios are smoke-class, independent-of-prior-state entry points registered side by side in the package, so they're covered in one file rather than two.
- **No related bugs.** No curated bug currently anchors the demo or app-browser surfaces. The closest sibling bug class is GROK-17557 (an init-prerequisite race on the context-panel Launch SAR path) — the demo and app paths don't hit the same race, since both go through their own dedicated entry handlers rather than the context-panel button.
- **Setup composition.** Both scenarios are expected to succeed from a clean shell with no prior state. Package init is allowed to run lazily on first invocation — the assertion is that both paths complete cleanly whether the package was already warmed up or not.
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Table View`; `public/packages/Peptides/README.md#Detection and usage`.

## Original trailing metadata

```json
{
  "order": 11,
  "datasets": []
}
```
