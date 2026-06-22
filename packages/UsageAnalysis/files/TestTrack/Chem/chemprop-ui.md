---
feature: chem
sub_features_covered: [chem.chemprop, chem.chemprop.train, chem.chemprop.apply, chem.chemprop.is-applicable, chem.chemprop.is-interactive]
target_layer: manual-only
coverage_type: regression
produced_from: manual-only-rename
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/chemprop.md
migration_date: 2026-05-11
related_bugs: []
ui_coverage_responsibility:
  - ml-models-train-model-top-menu
  - ml-train-model-feature-column-selector
  - ml-train-model-predict-column-selector
  - ml-apply-model
  - scatter-plot-add-for-correlation
ui_coverage_delegated_to: null
manual_only_rename_date: 2026-05-12
manual_execution_notes: |
  Manual-only execution required because the `chem-chemprop` Docker container
  is needed for training (multi-minute compute) and reliable container
  startup state cannot be guaranteed in CI Playwright runs. Per user
  directive 2026-05-12T17:41Z ("нужен запуск контейнера — пропускай этот
  тест и переводи в manual-only"). 8 rounds of Playwright Gate B attempts
  exhausted hypothesis budget (Round 1-3: column-selector widget; Round 4-5:
  JS API substitution attempts; Round 6: MCP-recon-validated UI patterns
  for Predict + Features confirmed but training step hung >12 min on
  Container start; Round 7: pre-cleanup unguarded page.evaluate hang;
  Round 8: bounded pre-cleanup + recon-validated patterns, but still
  exceeded budget due to container provisioning latency). The recon
  patterns (Predict via mousedown+search-input+dispatched-click,
  Features via None+dispatched-click row 0) are preserved in the
  decision-log for future re-enablement when container provisioning
  becomes reliable.
---

# Chem | Chemprop ML model training walk

End-to-end ML > Models > Train Model walk: open `mol1K.sdf`, launch the Train Model view, select
the molecule column as features and `pIC50_HIV_Integrase` as the prediction target, train a
Chemprop model, apply it back to the same dataset, then verify the prediction column tracks the
ground-truth `pIC50_HIV_Integrase` column via a scatter plot.

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario (`depends_on: []`),
`classification: simple`, `pyramid_layer: integration`, `target_layer: playwright`, strategy
`simple`. The defining test focus is the cross-subsystem coupling between the Datagrok ML
predictive-models registry (Train Model dialog, Apply Model), the Chem trainer
(`chem.chemprop.*` — train / apply / isApplicable / isInteractive), and the `chem-chemprop`
Docker container backend. UI coverage owned (`ui_coverage_delegated_to: null`) over the ML >
Models > Train Model entry, feature column selector, predict column selector, Apply Model
control, and scatter-plot correlation check.

## Setup

1. **Provision linked dataset.** `System:AppData/Chem/mol1K.sdf` (per source JSON footer
   `order: 13`, `datasets: ["System:AppData/Chem/mol1K.sdf"]`). The file is bundled in the
   platform's System file share — no external provisioning required. The dataset includes a
   molecule column plus the numeric activity column `pIC50_HIV_Integrase` used as the prediction
   target.
2. **Confirm Chem package and ML predictive-models surface are loaded** so that the **ML > Models
   > Train Model...** top-menu entry is registered and the Chemprop trainer (`chem.chemprop`,
   `mlname: Chemprop`, atlas line 1188; `chem.chemprop.train` package source
   `Chem/src/package.ts#L2556`) surfaces in the Train Model view for molecule-column input. The
   `chem-chemprop` Docker container must be running and reachable from the server — the train
   step calls into the container per atlas `chem.chemprop` description ("Message-passing-NN
   training & inference via the `chem-chemprop` Docker container").
3. **No fixture consumed.** Per chain YAML `depends_on: []`. The dataset is opened fresh inside
   the scenario and the trained model is in-session only (`produces` field of chain entry:
   "in-session only"); the scenario does not persist the model.

## Scenarios

### Train Chemprop on mol1K.sdf, apply, and verify prediction–activity correlation

Single linear walk (5 steps mirroring the original 5 numbered bullets, with expected-result
verifications woven inline per D-STEP-02):

1. Open `System:AppData/Chem/mol1K.sdf` (close any previously open views first to start from a
   clean state). Wait for the table view to render with the molecule column populated and the
   RDKit cell renderer applied. Verify the `pIC50_HIV_Integrase` column is present and is
   detected as numeric.
2. From the top menu, run **ML > Models > Train Model...**. The Train Model view opens (its
   active view name becomes `Predictive model`, mirroring the sibling
   `chemprop-spec.ts:46-47` `expect(viewName).toBe('Predictive model')` invariant). The
   `mlname: Chemprop` trainer is registered (`chem.chemprop.is-applicable` /
   `isApplicableNN` returns true for the molecule column, atlas
   `chem.chemprop.is-applicable` `Chem/src/package.ts#L2664`).
3. In the Train Model view, select the molecule column as the **features** input and
   `pIC50_HIV_Integrase` as the **predict** target. Choose the Chemprop trainer (`mlname:
   Chemprop`) — the cross-subsystem gating point per atlas `chem.chemprop.is-applicable`
   ("`isApplicableNN` — gate for showing Chemprop in the model dialog"). Verify the Chemprop
   trainer is selectable in the model picker; no console errors.
4. Run training. The Chemprop trainer dispatches to the `chem-chemprop` Docker container
   (atlas `chem.chemprop.train` — `trainChemprop(df, predictColumn, ...hyperparams)`,
   `mlrole: train`). Wait for training to complete (Chemprop training can take minutes — the
   sibling `chemprop-spec.ts:7` sets `test.setTimeout(300_000)` for this reason). Once the
   model is built, apply it back to the same dataset (atlas `chem.chemprop.apply` —
   `applyChemprop(df, model)`, `mlrole: apply`). Verify a new prediction column is appended to
   the table (commonly named `outcome` or a predict-column-derived label per the ML
   registry's apply convention); no console errors.
5. Add a scatter plot to the active view comparing the new prediction column against
   `pIC50_HIV_Integrase`. Verify the points cluster around the diagonal — the prediction
   column is "nearly equal" to `pIC50_HIV_Integrase` per the original assertion (the
   prediction-vs-actual correlation is the cross-subsystem correctness check that ties the ML
   train + apply round-trip back to the source dataset). No console errors during scatter
   creation; the visual correlation is the verifiable invariant.

Implicit cross-step invariants:

- The training step is the long-running step; the Apply step is fast once the model exists.
- The prediction column produced by Apply is positionally aligned with the source rows — the
  scatter plot can compare row-by-row without any join logic.
- The visual correlation (points clustering around the diagonal `y = x`) is the only
  perceptual judgment in the scenario; the Automator may quantify it with a Pearson / Spearman
  correlation coefficient or a per-row absolute-difference threshold for the scripted spec.

## Notes

- **`coverage_type: regression`** — happy-path end-to-end ML train + apply + correlation walk
  coupling the Datagrok ML predictive-models registry, the Chem Chemprop trainer, and the
  `chem-chemprop` Docker container backend. Not `smoke` (the section's smoke is
  `Advanced/scaffold-tree-functions.md` per chain `ui_coverage_plan.smoke_scenario`); not
  `edge` / `perf` (no specific failure-mode invariant being asserted here — the Docker
  container unavailability surface is owned by atlas critical path
  `chem.cp.calculate-descriptors-docker` p1 with the note "Chemprop training shares this
  surface" (atlas line 1538), but is not exercised by this happy-path walk). The scenario
  walks the happy training round-trip on a single dataset.
- **No JS API substitution.** Every entry in chain `ui_coverage_responsibility`
  (`ml-models-train-model-top-menu`, `ml-train-model-feature-column-selector`,
  `ml-train-model-predict-column-selector`, `ml-apply-model`,
  `scatter-plot-add-for-correlation`) is exercised via UI driving — top-menu walk + Train
  Model view interaction (feature column picker, predict column picker, trainer selector) +
  Apply Model control + scatter-plot add. Direct invocation of `trainChemprop` /
  `applyChemprop` is NOT a substitute for the dialog / view-driven flow that this scenario
  asserts; the Datagrok ML predictive-models registry integration (`mlname: Chemprop` surfaces
  in the Train Model view) is the cross-subsystem invariant.
- **Existing sibling spec.** A test file already exists at
  `public/packages/UsageAnalysis/files/TestTrack/Chem/chemprop-spec.ts` (per
  `existing-test-index.yaml` line 32445). It currently covers steps 1–2 only (open mol1K.sdf
  + open Train Model view) and stops short of column selection, training, apply, and scatter
  correlation. The Automator will extend it to cover the full 5-step walk per this migrated
  scenario. A second `chemprop-spec.ts` exists under
  `public/packages/UsageAnalysis/files/TestTrack/Models/` (per `existing-test-index.yaml`
  line 33060) — sibling Models-section test in the same family; cross-reference only,
  Chem-section spec is the canonical home for this scenario.
- **Docker container backend dependency.** The `chem-chemprop` Docker container must be
  running for training to complete. Atlas `chem.cp.calculate-descriptors-docker` (p1) cites
  this surface as a critical infrastructure-dependency path ("Container unavailability
  produces a clear error within bounded time (NOT infinite timeout per GROK-17621). Critical
  infrastructure-dependency path — Chemprop training shares this surface", atlas lines
  1531-1544). The unavailable-container error-bounding invariant is NOT exercised by this
  happy-path walk; it is the natural extension for an edge-coverage spec (parallel to
  `chem.cp.calculate-descriptors-docker` happy-path vs. failure-mode partition).
- **Long-running training timeout.** Chemprop training is a multi-minute operation per the
  sibling spec's `test.setTimeout(300_000)`. The Automator's extension MUST preserve a
  generous timeout for the train step (and possibly extend to 600s if Apply + scatter add
  push the spec over 5 minutes total).
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` + `closeAllViews` from
  `helpers-registry.yaml` cover the section harness. No registered helper currently abstracts
  the ML > Models > Train Model flow or the Chemprop train + apply round-trip; candidate
  helpers surfaced in migration report Decisions.
- **Order in chain.** `order: 13` per source JSON footer.
- **No source-text fixes applied during migration.** The original 5 numbered bullets are
  semantically clean (no step-numbering chaos, no typos, no JSON footer / body-prose
  divergence). The original is preserved verbatim into the 5-step migrated body with
  expected-result verifications woven inline per D-STEP-02. The original assertion
  "Make sure that column with predictions is nearly equal to `pIC50_HIV_Integrase` (use the
  scatterplot)" is the scenario's central correctness invariant and is preserved verbatim in
  Scenarios step 5.
