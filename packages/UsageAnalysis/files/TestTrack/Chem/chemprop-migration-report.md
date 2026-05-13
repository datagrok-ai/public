# Migration Report — chemprop.md

First-cycle migration of `chemprop.md` (TestTrack/Chem section) per chain YAML
`scenario-chains/chem.yaml` rev 2 (`output_plan.chemprop.md`: `target_layer: playwright`,
strategy `simple`, classification `simple`, `pyramid_layer: integration`, `depends_on: []`).
End-to-end Chemprop ML model training walk on `System:AppData/Chem/mol1K.sdf` covering five
chain `ui_coverage_responsibility` entries spanning ML > Models > Train Model entry, feature
column selector, predict column selector, Apply Model control, and scatter-plot correlation
check. Cross-subsystem coupling: Datagrok ML predictive-models registry ↔ Chem Chemprop
trainer (`chem.chemprop.*`) ↔ `chem-chemprop` Docker container.

## Step mapping

The original is a 14-line scenario (5 numbered steps + JSON footer) with clean step numbering
and no typos / phrasing defects. The body parses semantically into a single linear 5-step
walk (open → Train Model view → column selection → train + apply → scatter correlation
check). No source-text fixes required; original is preserved verbatim with expected-result
verifications woven inline per D-STEP-02.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Line 1 — "Open the file (Chem -> mol1k.sdf)." | Setup step 1 (provisioning — `System:AppData/Chem/mol1K.sdf` per JSON footer `order: 13`) + Scenarios step 1 (open dataset as the scenario entry; verify table view renders + molecule column + `pIC50_HIV_Integrase` numeric column present) | preserved (split for clarity — provisioning in Setup; open action in scenario body; expected-result verification added per D-STEP-02 — table renders + molecule column populated + RDKit cell renderer + `pIC50_HIV_Integrase` numeric) |
| Line 2 — "ML > Models > Train Model..." | Scenarios step 2 (top-menu walk → Train Model view opens; verify view name = `Predictive model` mirroring sibling chemprop-spec.ts:46-47 invariant; verify `mlname: Chemprop` registered + `isApplicableNN` returns true for the molecule column) | preserved (expected-result verification added per D-STEP-02 — view name + Chemprop trainer registration gate `chem.chemprop.is-applicable`) |
| Line 3 — "Select 'molecule' column as features and 'pIC50_HIV_Integrase' for prediction." | Scenarios step 3 (Train Model view: select molecule column as features input + `pIC50_HIV_Integrase` as predict target + Chemprop trainer in model picker; verify trainer selectable + no console errors) | preserved (expected-result verification added per D-STEP-02 — Chemprop trainer selectable + atlas cross-reference `chem.chemprop.is-applicable` cited; the original phrasing "molecule column" is preserved verbatim as the column selector target) |
| Line 4 — "Run the training, when the model is built run prediction for the same dataset." | Scenarios step 4 (train via Chemprop → `chem-chemprop` Docker container; wait for completion; apply model back to same dataset; verify prediction column appended; no console errors) | preserved (split semantically into train + apply phases within one step per the original's single-bullet phrasing; expected-result verifications added per D-STEP-02 — prediction column appended + atlas `chem.chemprop.train` / `chem.chemprop.apply` cross-references; long-running training timeout note moved to migrated Notes section) |
| Line 5 — "Make sure that column with predictions is nearly equal to 'pIC50_HIV_Integrase' (use the scatterplot)." | Scenarios step 5 (add scatter plot prediction-column vs. `pIC50_HIV_Integrase`; verify points cluster around diagonal; no console errors) | preserved (expected-result verification preserved verbatim — "Make sure that column with predictions is nearly equal to `pIC50_HIV_Integrase` (use the scatterplot)" mapped to the prediction–actual diagonal-clustering invariant per D-STEP-02; visual correlation quantification options surfaced in migrated body's "Implicit cross-step invariants" note for Automator) |
| JSON footer `{ "order": 13, "datasets": ["System:AppData/Chem/mol1K.sdf"] }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` + migrated Setup step 1 dataset enumeration) |

No original step is silently dropped. All 5 numbered bullets in the original body land in
either Setup step 1 (Line 1 provisioning half) or Scenarios steps 1–5. No source-text defect
present; no silent fix applied per D-SAN-02.

## Decisions

- **Why this `target_layer`:** chose `playwright` per `scenario-chains/chem.yaml`
  `output_plan.chemprop.md.target_layer = playwright`. The scenario requires DOM-level UI
  driving: top-menu **ML > Models > Train Model...** entry click, Train Model view
  interaction (feature column picker, predict column picker, trainer selector), Apply Model
  control, and scatter-plot add for the prediction-vs-actual correlation visualization —
  all UI behaviors that Playwright drives natively. The cross-subsystem invariant (Datagrok
  ML predictive-models registry surfaces `mlname: Chemprop` in the Train Model view via
  `chem.chemprop.is-applicable`) specifically requires the dialog-driven flow; direct
  `trainChemprop` / `applyChemprop` invocation via JS API does not exercise the registry
  integration gate. Chain YAML reason field: "ML > Models > Train Model walk coupling ML
  predictive-models registry + Chem chemprop train/apply + chem-chemprop Docker container.
  Sibling chemprop-spec.ts already exists at playwright per existing-test-index. Docker-
  backed training flow with UI interaction routed to playwright. Single scenario; simple
  strategy."
- **Why this `coverage_type`:** chose `regression` per chain rev 2 (`pyramid_layer:
  integration` — cross-subsystem coupling between ML registry / Chem trainer / Docker
  container; `classification: simple` — single linear walk, no matrix). Not `smoke` (the
  section's smoke is `Advanced/scaffold-tree-functions.md` per chain
  `ui_coverage_plan.smoke_scenario`; chemprop's `ui_coverage_delegated_to: null` records
  "Scaffold Tree smoke does not witness ML flows. Scenario owns its own coverage."). Not
  `edge` / `perf` since the scenario walks the happy training round-trip only — the Docker
  container unavailability error-bounding invariant per `GROK-17621` (atlas
  `chem.cp.calculate-descriptors-docker` p1 with the note "Chemprop training shares this
  surface", atlas line 1538) is NOT exercised by this happy-path walk; the long-running
  training timeout invariant is a runtime concern (Automator handles via
  `test.setTimeout`), not a coverage-type signal. `regression` is the natural fit for a
  cross-subsystem happy-path walk that guards against regression-of-the-integration across
  the ML registry ↔ Chem trainer ↔ Docker container surface.
- **Why this `strategy`:** `simple` per chain YAML `output_plan.chemprop.md.strategy =
  simple`. Single scenario, single dataset, no cross-file fixture; chain analyzer
  classified `simple` (single linear 5-step walk, no matrix axes, no decomposition
  proposal). Pattern 1 (decomposition) is NOT triggered (chain classifier returned
  `simple`, not `complex-standalone`; no child scenarios with `produced_from: decomposed`).
  Pattern 2 (bug-focused slice) does NOT apply — `chemprop.md` is NOT a span in any of the
  chain's 10 `bug_focused_candidates[]` (verified via grep on chain YAML). Pattern 3
  (UI-smoke consolidation) is NOT triggered (no `ui_consolidation_proposals[]` from Gate
  F). Standard per-scenario migration per "Migrated scenario shape".
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/chemprop-spec.ts` — existing
    playwright-layer test ("Chem: ChemProp training", category Chem, `test.setTimeout(300_000)`).
    Currently covers steps 1–2 only (open mol1K.sdf via `grok.data.files.openTable` + open
    Train Model view; assertion `viewName === 'Predictive model'`). The spec stops short of
    column selection, training, apply, and scatter correlation — Automator at extension time
    will extend it to the full 5-step walk. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Models/chemprop-spec.ts` — sibling
    Models-section playwright test (per `existing-test-index.yaml` line 33060) in the same
    chemprop family; cross-reference only, Chem-section spec is the canonical home for this
    scenario. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/chemical-space-migration-report.md`
    and `chemical-space.md` (2026-05-11 predecessor in same migration cycle) — shape anchor
    for the migrated `.md` (single-scenario body + `coverage_type: regression` precedent
    for a cross-subsystem integration walk + SR-01 carryforward block). Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/filter-panel-migration-report.md`
    and `filter-panel.md` (2026-05-11 predecessor in same migration cycle) — shape anchor
    for the in-line expected-result verification weave + Implicit cross-step invariants
    note pattern. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/elemental-analysis-migration-report.md`
    (2026-05-11 predecessor in same migration cycle) — SR-01 A-STRUCT-02 carryforward
    template precedent (simple single-scenario `coverage_type: regression`). Read-only.
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`)
      — section-standard fixture login; Automator will use in spec `beforeEach`.
    - `softStep` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`) —
      per-step soft-fail wrapper; Automator wraps each of the 5 walk steps in `softStep`
      for granular failure reporting (Train step is the long-running step where soft-fail
      visibility matters most).
    - `closeAllViews` (`public/packages/UITests/playwright/scripts/helpers.ts:129`) —
      cleanup at the end of the scenario, mirroring sibling `chemprop-spec.ts:50`
      `await page.evaluate(() => grok.shell.closeAll())` pattern.
  - **Candidate helpers (NOT yet in registry — flagged for addition via helpers-registry
    curator; per migration-prompt Helpers discipline §):**
    - `helpers.playwright.ml.openTrainModelView(page)` — top-menu walk
      **ML > Models > Train Model...** and waits for the `Predictive model` view to
      render. Surfaced by Scenarios step 2. Sibling
      `chemprop-spec.ts:36-48` already inlines this pattern via the
      `[name="div-ML"]` selector + `.d4-menu-item-label` lookup; convenience abstraction
      candidate.
    - `helpers.playwright.ml.selectTrainModelInputs(page, {features, predict, trainer})`
      — select feature column, predict column, and trainer in the Train Model view.
      Surfaced by Scenarios step 3.
    - `helpers.playwright.ml.runTrainAndApply(page, options)` — run training (wait up to
      300s) then Apply Model back to the active table; returns the appended prediction
      column name. Surfaced by Scenarios step 4. Long-running operation — helper must
      preserve generous timeout.
    - `helpers.playwright.viewers.addScatterPlot(page, {x, y})` — add a scatter plot
      viewer comparing two columns on the active table view. Surfaced by Scenarios step 5.
- **Bug library consulted:** yes — `bug-library/chem.yaml` revision 1 (12 curated_bugs,
  generated 2026-05-05). Grepped for `chemprop` — **zero matches**. No curated bug
  intersects this scenario's `sub_features_covered` (`chem.chemprop`, `chem.chemprop.train`,
  `chem.chemprop.apply`, `chem.chemprop.is-applicable`, `chem.chemprop.is-interactive`).
  `related_bugs: []` in frontmatter.
- **Decision log queried:** yes — `decision-log.yaml` grepped for `feature: chem` — zero
  per-scenario `chemprop`-specific entries. The chem section is in its first migration
  cycle. No prior `migration_decisions`, `layer_decisions`, `manual_only`, or
  `failed_attempts` entries exist for chemprop. No "approaches off the table" constraints
  apply.
- **Cross-cutting bug citations (chain YAML `bug_focused_candidates[]`).** Per
  `migration-prompt.md` "Cross-cutting bug citations from chain YAML" §: `chemprop.md` is
  NOT a span in any of the chain's 10 `bug_focused_candidates[]` (verified via grep on the
  chain YAML — chemprop does not appear in any `spans:` list). No cross-cutting bug
  citation applicable. F-BUG-COVERAGE-01 at section-complete is authoritative for cross-
  cutting bug coverage; no migrator-side citation required here.
- **Atlas critical-path adjacency.** This scenario covers `chem.chemprop.*` sub-features but
  is NOT directly named in any atlas `critical_paths[]` entry. The infrastructure-dependency
  adjacency is via `chem.cp.calculate-descriptors-docker` (p1) which explicitly cross-
  references Chemprop training in its description ("Critical infrastructure-dependency path
  — Chemprop training shares this surface", atlas line 1538). The Docker container
  unavailability error-bounding invariant per `GROK-17621` is owned by that critical-path
  entry; this scenario walks the happy training round-trip only and intentionally does NOT
  exercise the container-down failure mode. Surfaced in migrated Notes (`Docker container
  backend dependency`) for awareness.
- **UI delegation status.** Per chain YAML `ui_coverage_plan.smoke_scenario:
  Advanced/scaffold-tree-functions.md` (NOT this scenario) and chain
  `ui_coverage_plan.delegated_scenarios[] → chemprop.md → delegated_to: null` ("UI
  coverage owned (Datagrok ML > Models > Train Model dialog + feature/predict column
  selectors + Apply Model + scatter plot correlation check). Cross-subsystem (ML + Chem +
  Docker). Scaffold Tree smoke does not witness ML flows. Scenario owns its own
  coverage."), this scenario's 5 `ui_coverage_responsibility` flows are owned in-scenario.
  No SCOPE_REDUCTION proposal substitutes JS API for any of the 5 UI flows; therefore the
  UI delegation discipline of migration-prompt § "UI delegation in SCOPE_REDUCTION
  proposals" does not apply (no SR proposal to delegate). Reaffirmed as a hard constraint
  in the migrated Notes ("No JS API substitution").
- **Scenario constraint extraction (per migration-prompt § Scenario constraint
  extraction).**
  - **(a) FORBIDDEN substitutions:** None explicit in the original body (no "must not",
    "forbidden", "not allowed" phrasing). HOWEVER `pyramid_layer: integration` per chain +
    `ui_coverage_delegated_to: null` (scenario owns its own coverage) + the cross-subsystem
    coupling invariant (ML registry ↔ Chem trainer ↔ Docker container) imply JS API
    substitution for the 5 `ui_coverage_responsibility` flows is NOT acceptable — the
    Datagrok ML predictive-models registry integration gate (`mlname: Chemprop` surfaces
    via `chem.chemprop.is-applicable`) is observable only through the Train Model view;
    direct `trainChemprop` / `applyChemprop` invocation bypasses the registry. Surfaced
    explicitly in the migrated Notes (`No JS API substitution`) as a structured
    behavioral guard.
  - **(b) REQUIRED actions:** All 5 entries of `ui_coverage_responsibility` are exercised
    via UI driving in the migrated body — ML > Models > Train Model top-menu walk
    (Scenarios step 2), feature column selector (Scenarios step 3), predict column
    selector (Scenarios step 3), Apply Model control (Scenarios step 4),
    scatter-plot-add-for-correlation (Scenarios step 5). No flow is deferred.
  - **(c) Missing-selector escalation:** Migrator does NOT have explicit selector
    definitions for the Train Model view's feature column picker, predict column picker,
    trainer selector, train run control, Apply Model control, or scatter-plot add-via-
    toolbar button. The sibling `chemprop-spec.ts:36-48` already drives the top-menu
    via `[name="div-ML"]` and `.d4-menu-item-label` selectors — Automator at extension
    time will inherit and extend; spec-time selector discovery for the unexplored steps
    (3–5) is the Automator's responsibility (per existing sibling spec patterns + `page.
    evaluate` introspection). A reference-file approval-required proposal to add ML-Train-
    Model-specific selectors to `grok-browser/references/ml-train-model.md` (or
    equivalent) is a Phase 2 deliverable — out of this migration's scope.
  - **(d) Reference templates:** Code-style anchors for the Automator (cited verbatim,
    not modified here): `public/packages/UsageAnalysis/files/TestTrack/Chem/chemprop-
    spec.ts` (current existing test — sibling shape; `test.setTimeout(300_000)` long-
    running pattern + `[name="div-ML"]` ML top-menu selector) +
    `public/packages/UsageAnalysis/files/TestTrack/Models/chemprop-spec.ts` (sibling
    Models-section playwright test in the same chemprop family) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/chemical-space-spec.ts`
    (section-mate playwright walk; `spec-login` + `page.evaluate` + `page.locator`
    patterns for a Chem-section dialog-driven analyzer flow).
  - **(e) Scenario authority clause:** No prompt-vs-scenario conflict surfaced. The
    original scenario is a clean 5-bullet linear walk; the migrated body preserves all
    original bullets + adds inline expected-result verifications per D-STEP-02. No
    `prompt_scenario_conflict:` field to flag for orchestrator.
- **Source-text fixes silently applied (per chain rev 2 directive footer note (c) — Olena
  2026-05-11).**
  - **None.** The original 5 numbered bullets are semantically clean — no step-numbering
    chaos, no typos, no `Drag'n'drop` / `sketcherr`-style canonical-phrasing fixes, no
    JSON footer / body-prose divergence. The original is preserved verbatim into the
    5-step migrated body with expected-result verifications woven inline per D-STEP-02.
- **No invented sub_features / helpers.** Every id in `sub_features_covered` is present in
  `feature-atlas/chem.yaml` rev 2 (verified per atlas lines 1188, 1195, 1202, 1209, 1216
  for `chem.chemprop` + `chem.chemprop.train` + `chem.chemprop.apply` +
  `chem.chemprop.is-applicable` + `chem.chemprop.is-interactive` — 5 ids total). Every
  reused helper is in `helpers-registry.yaml`; candidate helpers are flagged as
  candidates, NOT referenced by name in the migrated body.

## Opt-outs (SCOPE_REDUCTION proposals)

### SR-01: A-STRUCT-02 carryforward (chain-level edge/perf coverage)

**Cited technical dependency:** A-STRUCT-02 ("at least one scenario is edge case or
negative path — i.e. its frontmatter `coverage_type` is `edge` or `perf`") is a
**section-level structural invariant**. This scenario (`chemprop.md`) is a happy-path
cross-subsystem ML train + apply + correlation walk — `coverage_type: regression` is the
natural fit per chain rev 2 (`pyramid_layer: integration`, `classification: simple`, single
linear 5-step walk on a single dataset). Forcing `coverage_type: edge` would mis-classify
(the scenario asserts happy-path Chemprop training round-trip on `mol1K.sdf`, not a
specific failure-mode invariant such as Docker container unavailability error bounding per
`GROK-17621`, infeasible-prediction-column error, or training timeout). The container
unavailability failure mode is owned by atlas critical path
`chem.cp.calculate-descriptors-docker` p1 with the note "Chemprop training shares this
surface" (atlas line 1538) — natural extension for an edge-coverage spec, but distinct from
this happy-path migration.

**Resolution path:** A-STRUCT-02 satisfaction is chain-wide via 10 `bug_focused_candidates[]`
in `scenario-chains/chem.yaml` rev 2 (each chem-bug-focused spec naturally carries
`coverage_type: edge` once authored — none of the 10 directly target Chemprop, but the
section-level invariant is satisfied chain-wide). Additionally `r-group-analysis.md`
(same-cycle, same section, `coverage_type: edge`) already organically satisfies the
section-level invariant. Same SR-01 carried by sibling-cycle reports
`activity-cliffs-migration-report.md`, `calculate-migration-report.md`,
`scaffold-tree-functions-migration-report.md`, `scaffold-tree-migration-report.md`,
`info-panels-migration-report.md`, `Advanced/structure-filter-migration-report.md`,
`Advanced/similarity-search-migration-report.md`, `elemental-analysis-migration-report.md`,
`chemical-space-migration-report.md`, `sketcher-migration-report.md`,
`filter-panel-migration-report.md`.

**Scope of waiver:** A-STRUCT-02 verdict on `chemprop.md` is deferred to chain-level
evaluation. Critic A should return SR (carryforward), not FAIL.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through A-STRUCT-06, A-COVERAGE-*,
A-MERIT-*) are expected to PASS on `chemprop.md` without any other SR proposals. No step
is opted out for effort. The single entry above cites a real section-level structural
property whose satisfaction path is owned by the chain (and reinforced organically by
`r-group-analysis.md`'s `coverage_type: edge`), not by this per-scenario migration
(D-MERIT-01 compliant)._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist yet. The migrated body's
5 steps (Setup step 1 reused as Scenarios step 1 entry + 4 walk steps) are realizable
against the current playwright + helpers-registry + atlas state. The candidate helpers
surfaced in Decisions are convenience abstractions that the Automator may inline pending
registration — they do NOT block spec realization. The selector-reference write proposal
(Decisions § Scenario constraint extraction (c)) is a downstream Phase 2 enhancement, not
a blocker on this scenario's automation (existing `chemprop-spec.ts` already drives the
ML top-menu via `page.evaluate` + selectors per the sibling-test pattern; Automator can
extend that pattern without a reference update). The long-running training timeout
(multi-minute) is a runtime concern handled by `test.setTimeout` in the sibling spec —
not a Migrator-side deferral._

## Edge cases

The original lists no explicit "edge case" keyword. Implicit edge cases derivable from the
scenario content:

- **Long-running training timeout.** Chemprop training is a multi-minute Docker-container-
  backed operation. Sibling `chemprop-spec.ts:7` sets `test.setTimeout(300_000)` (5
  minutes); the Automator's extension MUST preserve a generous timeout for the train step
  (Step 4 in this migration) and possibly extend to 600s if Apply + scatter add push the
  spec over 5 minutes total. PRESERVED as scenario step + Notes-level caveat (migrated
  Notes "Long-running training timeout"). Not a coverage-type signal — runtime concern.
- **Docker container unavailability.** Implicit edge case surfaced by atlas
  `chem.cp.calculate-descriptors-docker` (p1) which cross-references Chemprop training:
  "Container unavailability produces a clear error within bounded time (NOT infinite
  timeout per GROK-17621)". The scenario does NOT exercise this failure mode directly —
  it walks the happy path only on the assumption that `chem-chemprop` is up and reachable.
  PRESERVED-AS-NOTE (no scenario step here; the error-bounding invariant is owned by the
  atlas critical path entry and would be the natural target for an edge-coverage spec).
  Atlas critical_path `chem.cp.calculate-descriptors-docker` already captures the
  cross-Chemprop adjacency — no atlas write required.
- **Prediction column naming convention.** The prediction column appended by Apply Model
  is auto-named by the ML registry per the trainer's convention (commonly `outcome` or a
  predict-column-derived label like `pIC50_HIV_Integrase (predicted)` depending on the
  trainer / registry version). PRESERVED as scenario step (Scenarios step 4 — explicit
  "prediction column appended (commonly named `outcome` or a predict-column-derived label
  per the ML registry's apply convention)"). Automator at spec time picks the actual
  column name by introspection. Implicit edge case: the column name is auto-determined
  and may vary across registry versions — discovery via column-difference between pre-
  and post-apply column lists is the robust pattern.
- **Visual correlation quantification.** Original Line 5 "Make sure that column with
  predictions is nearly equal to `pIC50_HIV_Integrase` (use the scatterplot)" is a
  perceptual judgment ("nearly equal" + "use the scatterplot"). PRESERVED as scenario
  step (Scenarios step 5 — "points cluster around the diagonal" verbatim assertion +
  Implicit cross-step invariants note "the Automator may quantify it with a Pearson /
  Spearman correlation coefficient or a per-row absolute-difference threshold"). Implicit
  edge case: the perceptual "nearly equal" assertion is non-trivial to automate without
  picking a quantification threshold; the Automator's choice of correlation coefficient +
  threshold (e.g. Pearson r > 0.7, or per-row max-abs-diff < 1.0 on the pIC50 scale) is
  the spec-side decision.
- **Single-dataset coverage.** Original is single-dataset (`mol1K.sdf` only). PRESERVED
  as single-dataset migration body (no multi-format extension applied — distinct from
  `activity-cliffs.md` / `chemical-space.md` / `calculate.md` / `elemental-analysis.md`
  multi-format extensions which were driven by source-body prose "should be tested on
  smiles, molV2000, molV3000 formats"). chemprop's original body has no such prose; the
  scenario's source defines the dataset surface as `mol1K.sdf` only. Multi-format
  Chemprop training is a natural extension for a future edge-coverage or multi-dataset
  spec, but is out of scope for this happy-path migration.
- **Console errors throughout.** Implicit across all 5 original bullets: each action
  (open dataset, top-menu walk, column selection, train + apply, scatter add) must
  complete without console errors. PRESERVED as per-step verification weave (Scenarios
  steps 3 + 4 + 5 explicitly include "no console errors" assertions; steps 1 + 2
  implicitly assume same).

No edge case is moved to atlas, manual_only, deferred, or a separate scenario. All
implicit edge cases are PRESERVED as scenario steps OR Notes-level caveats in the
migrated body, with the Docker container unavailability error-bounding invariant flagged
for future edge-coverage spec consideration (already owned by atlas
`chem.cp.calculate-descriptors-docker` p1).

## Unresolved ambiguities

- **"Run the training, when the model is built run prediction" — separate user actions or
  single fused action?** Original Line 4 conjoins train + apply in a single numbered
  bullet ("Run the training, when the model is built run prediction for the same
  dataset."). Migrated body resolves by treating them as a single semantic step (Scenarios
  step 4) with the train + apply phases described within the step body — the cross-
  subsystem train + apply round-trip is the unit of interest, and Apply Model is
  triggered by the user after train completes (not auto-triggered by training itself).
  Some predictive-model UIs auto-display an Apply control once training succeeds; others
  require a separate top-menu invocation **ML > Models > Apply Model...** Automator at
  spec time picks the actual UI invocation pattern by observation. Flag for QA pair
  review.
- **"Column with predictions" — exact column name produced by Apply.** Original Line 5
  references "column with predictions" without naming it. The ML registry's apply
  convention auto-names the column (commonly `outcome` or a predict-column-derived label
  like `pIC50_HIV_Integrase (predicted)`). Migrated body resolves by deferring to spec-
  time introspection (column-difference between pre- and post-apply column lists is the
  robust pattern). Automator at spec time picks the actual column name. Flag for QA pair
  review.
- **"Nearly equal" — quantification threshold.** Original Line 5 ("Make sure that column
  with predictions is nearly equal to `pIC50_HIV_Integrase` (use the scatterplot)") is a
  perceptual judgment. Migrated body resolves by preserving the verbatim assertion +
  surfacing Pearson / Spearman correlation or per-row absolute-difference threshold as
  Automator-side quantification options. Automator at spec time picks the specific
  threshold (e.g. Pearson r > 0.7, or median per-row |delta| < 1.0 on the pIC50 scale).
  Implicit assumption: the trained model is "good enough" to clear whatever threshold
  the Automator picks — `mol1K.sdf` with `pIC50_HIV_Integrase` is the canonical Chemprop
  training dataset, and the original scenario assumes the model will train
  successfully and produce a high-correlation prediction column; if training fails or
  produces a low-correlation prediction (model under-fit), the scenario's verification
  fails — this is the intended behavior. Flag for QA pair review.
- **"Add scatter plot" — toolbar button vs. menu entry vs. context-menu vs. shortcut.**
  Migrated Scenarios step 5 says "Add a scatter plot to the active view comparing the new
  prediction column against `pIC50_HIV_Integrase`" without specifying the invocation
  pattern. Datagrok has multiple paths: viewer-strip add button on the active view,
  **Add | Scatter Plot** top-menu, right-click context-menu, or a keyboard shortcut.
  Original is silent. Automator at spec time picks one deterministic invocation
  pattern (most likely the viewer toolbar add button on the active table view, mirroring
  the sibling Chem-section scatter-plot patterns). Flag for QA pair review.
