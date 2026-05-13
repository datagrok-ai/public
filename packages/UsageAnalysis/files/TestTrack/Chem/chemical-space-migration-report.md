# Migration Report — chemical-space.md

First-cycle migration of `chemical-space.md` (TestTrack/Chem section, original order 8) per chain YAML
`scenario-chains/chem.yaml` rev 2 (`output_plan.chemical-space.md`: `target_layer: playwright`, strategy
`simple`, classification `medium`, `pyramid_layer: integration`, `depends_on: []`). Multi-format Chemical
Space dialog + embedding-transform walk across SMILES / molV2000 / molV3000 dataset variants. Overlaps
atlas critical path `chem.cp.chemical-space` (p1); parallel-coverage with bug-focused candidate
`chem-grok-18407-spec.ts` on the same `chem.analyze.chemical-space.*` sub_features (different invariant
scope — happy path here, empty-column-input edge case there).

## Step mapping

The original is a 13-line scenario with 6 numbered steps plus a JSON footer with `order: 8` and 1 dataset
(`System:DemoFiles/chem/smiles.csv`). Body Step 1 prose says "should be tested on smiles, molV2000,
molV3000 formats" — JSON/body mismatch flagged in chain rev 2 `dependency_graph.chemical-space.md.notes`.
Per chain rev 2 directive footer note (c) "dataset list extension" (Olena 2026-05-11, chain footer note
(c) "source-text defects in scenario bodies — Migrator fixes during migration"), the dataset list is
extended to 3 entries here to match body prose, following sibling activity-cliffs.md (5-entry extension)
/ calculate.md / elemental-analysis.md (3-entry extension) precedent in the same cycle.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Step 1 — "Open linked datasets (should be tested on smiles, molV2000, molV3000 formats)" | Setup step 1 (provisioning — datasets D1 / D2 / D3 enumerated for the three notation formats) + Scenarios > "Chemical Space dialog walk per dataset" step 1 (dataset open as the per-cell loop entry) | preserved (split for clarity — provisioning declaration in Setup; the open action in scenario body; JSON dataset list extended from 1 entry to 3 entries to match body prose per chain rev 2 directive) |
| Step 2 — "On the menu ribbon, open **Chem > Analyze > Chemical space**. A dialog opens." | Scenarios > "Chemical Space dialog walk per dataset" step 2 (top-menu walk → `ChemSpaceEditor` dialog open + column-picker + method selector + Cluster MCS toggle surfaced) | preserved (with explicit dialog-shape verification per D-STEP-02 — naming the four editor surfaces per atlas `chem.analyze.chemical-space.editor`) |
| Step 3 — "Click OK to run a function with default parameters." | Scenarios > "Chemical Space dialog walk per dataset" step 3 (OK click with defaults; verify 2D scatter rendered with `chem-space-embedding-col` tagged columns appended; no console errors) | preserved (with expected-result verification made explicit per D-STEP-02 — embedding output column tagging is part of the atlas `chem.analyze.chemical-space.transform` contract) |
| Step 4 — "Run **Chemical space**." | Scenarios > "Chemical Space dialog walk per dataset" step 4 (re-run Chemical Space top-menu on the same active table) | preserved (interpreted as the second-pass entry per the body's implicit two-pass structure — default-pass then custom-pass) |
| Step 5 — "In the dialog, change the parameters arbitrarily." | Scenarios > "Chemical Space dialog walk per dataset" step 5 (change at least one editor parameter — method UMAP↔t-SNE, similarity metric, or Cluster MCS toggle) | preserved (arbitrary-change made concrete via the enumerated parameter list per atlas `chem.analyze.chemical-space.editor` — same approach as activity-cliffs.md step 10) |
| Step 6 — "Click OK to run a function with edited parameters." | Scenarios > "Chemical Space dialog walk per dataset" step 6 (OK click with edited params; verify fresh embedding scatter reflects edited parameters; Cluster-MCS-on adds `chem-space-cluster-col`; no console errors) | preserved (renumbered cleanly; expected-result verification added per D-STEP-02 — point-distribution difference and conditional cluster-col column append) |
| (implicit — close active view before next dataset) | Scenarios > "Chemical Space dialog walk per dataset" step 7 (close active view between cells) | preserved as cleanup (implicit per-cell isolation invariant made explicit) |
| JSON footer `{ "order": 8, "datasets": ["System:DemoFiles/chem/smiles.csv"] }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` + migrated Setup step 1 dataset enumeration — extended to 3 entries to match body Step 1 prose) |

No original step is silently dropped. All 6 original numbered steps are preserved (no source-text
defects in step numbering for this scenario, unlike sibling activity-cliffs.md). The migrated body adds
one cleanup step (close active view) to make the implicit per-cell isolation invariant explicit. The
dataset-list extension (1 → 3 entries) is the sole silent fix per chain rev 2 directive.

## Decisions

- **Why this `target_layer`:** chose `playwright` per `scenario-chains/chem.yaml`
  `output_plan.chemical-space.md.target_layer = playwright`. The scenario requires DOM-level UI
  driving: top-menu walk (**Chem | Analyze | Chemical Space...**), `ChemSpaceEditor` dialog
  interaction with column picker + method selector + Cluster MCS toggle, and re-run with edited
  parameters. Chain YAML reason field: "Chemical Space dialog walk — default then custom parameters.
  Sibling chemical-space-spec.ts already exists at playwright per existing-test-index. Dialog +
  method selector + Cluster MCS toggle routed to playwright. Single scenario; simple strategy."
- **Why this `coverage_type`:** chose `regression` per chain rev 2 (`pyramid_layer: integration` —
  multi-format × dialog-walk; not `smoke` since the section's smoke is
  `Advanced/scaffold-tree-functions.md` per chain `ui_coverage_plan.smoke_scenario`; not `edge` /
  `perf` since the scenario walks the happy path only — the empty-column-input silent-close
  failure-mode invariant (atlas `chem.cp.chemical-space` p1 edge case, `GROK-18407`) is owned by
  bug-focused candidate `chem-grok-18407-spec.ts` per chain `bug_focused_candidates[]`, not by this
  scenario). `regression` is the natural fit for a multi-format dialog walk that guards against
  regression-of-the-set across the Chemical Space dialog + transform surface. Sister-pattern to
  same-cycle predecessors `activity-cliffs.md`, `calculate.md`, `elemental-analysis.md` (all
  `coverage_type: regression` per chain rev 2).
- **Why this `strategy`:** `simple` per chain YAML `output_plan.chemical-space.md.strategy = simple`.
  Single scenario, no cross-file fixture; chain analyzer classified `medium` (multi-format prose
  promise + dialog interaction) but kept strategy `simple` — single migrated `.md`, single downstream
  spec. Pattern 1 (decomposition) is NOT triggered (chain classifier returned `medium`, not
  `complex-standalone`; no child scenarios with `produced_from: decomposed`). Pattern 2 (bug-focused
  slice) does NOT apply to this scenario directly — `chemical-space.md` IS a span in chain
  `bug_focused_candidates[chem-grok-18407-spec.ts].spans` (steps 3 and 5), but Pattern 2 governs the
  bug-focused spec authoring (separate, parallel-coverage), not this happy-path migration. Pattern 3
  (UI-smoke consolidation) is NOT triggered (no `ui_consolidation_proposals[]` from Gate F).
  Standard per-scenario migration per "Migrated scenario shape".
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/chemical-space-spec.ts` — existing
    playwright-layer test ("Chem: Chemical Space", `category: Chem`, `layer: playwright` per
    `existing-test-index.yaml`). The Automator will extend it to cover the full enumerated 3-format
    dataset variant matrix per this migrated scenario. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/activity-cliffs.md` +
    `activity-cliffs-migration-report.md` (2026-05-11 same-cycle predecessor — adjacent in chain
    order, both `chem.analyze.*` analyze-dialog walks). Shape anchor for the migrated `.md`
    (multi-format dataset matrix + per-cell dialog walk + `coverage_type: regression` precedent +
    JSON-footer-vs-body-prose extension precedent + SR-01 A-STRUCT-02 carryforward template). The
    two are sister scenarios — same analyze-dialog UX pattern (top-menu walk → editor dialog →
    OK with defaults → re-run with edited params), different sub_feature scope. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/calculate.md` +
    `calculate-migration-report.md` (2026-05-11 same-cycle predecessor — 3-format matrix precedent).
    Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/elemental-analysis.md` +
    `elemental-analysis-migration-report.md` (2026-05-11 same-cycle predecessor — 3-format matrix
    precedent + SR-01 carryforward template). Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts` — sibling Chem
    playwright test ("Chem: Info Panels on smiles.csv"). Same shape pattern; confirms the section-
    wide playwright convention. Read-only.
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`) —
      section-standard fixture login; Automator will use in spec `beforeEach`.
    - `softStep` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`) — per-step
      soft-fail wrapper; Automator wraps each cell verification in `softStep` for granular failure
      reporting (especially important for a 3-variant × 7-step dialog walk where individual cells
      may surface format-specific dialog branches).
    - `closeAllViews` (`public/packages/UITests/playwright/scripts/helpers.ts:129`) — cleanup at the
      end of each dataset cell + at the end of the scenario.
  - **Candidate helpers (NOT yet in registry — flagged for addition via helpers-registry curator;
    per migration-prompt Helpers discipline §):**
    - `helpers.playwright.chem.openChemicalSpaceDialog(page, dataset, params?)` — opens dataset +
      walks **Chem | Analyze | Chemical Space...** top-menu + accepts dialog with default or
      supplied params; returns a locator for the resulting embedding scatter plot. Surfaced by
      Scenarios steps 1-3. Convenience abstraction; spec Author can inline the top-menu walk until
      the helper lands.
    - `helpers.playwright.chem.changeChemicalSpaceMethod(editor, method)` — changes the method
      selector between UMAP and t-SNE in the open `ChemSpaceEditor` dialog. Surfaced by Scenarios
      step 5.
    - `helpers.playwright.chem.toggleClusterMCS(editor)` — toggles the Cluster MCS checkbox in the
      open `ChemSpaceEditor` dialog. Surfaced by Scenarios step 5.
- **Bug library consulted:** yes — `bug-library/chem.yaml` revision 1 (12 curated_bugs, generated
  2026-05-05). One curated bug intersects this scenario's `sub_features_covered`:
  - `GROK-18407` (Chemical Space: do not allow to run analysis with empty column data) —
    `affects: [chem.analyze.chemical-space, chem.analyze.chemical-space.top-menu,
    chem.analyze.chemical-space.editor, chem.analyze.chemical-space.transform]`,
    `status: regression-risk`, `test_coverage: needed`. The bug is OWNED by chain
    `bug_focused_candidates[chem-grok-18407-spec.ts]` (Pattern 2 bug-focused spec). Per chain rev 2
    `dependency_graph.chemical-space.md.notes` ("Body does NOT test empty-column-input edge case (the
    GROK-18407 silent-close invariant) — this gap is captured under bug_focused_candidates[]"), this
    scenario does NOT propagate the bug into `related_bugs: []`. The bug-focused spec is the
    authoritative coverage owner for that failure-mode invariant; this scenario provides parallel
    happy-path coverage on the same sub_features. `related_bugs: []` in the migrated frontmatter;
    the bug is cited in Notes as parallel-coverage awareness only, NOT propagated as `related_bugs`.
- **Decision log queried:** yes — `decision-log.yaml` grepped for `feature: chem`. No prior
  `migration_decisions`, `layer_decisions`, `manual_only`, or `failed_attempts` entries exist for
  this specific scenario file (chemical-space.md). The same-cycle predecessors
  (activity-cliffs-migration-report.md, calculate-migration-report.md, elemental-analysis-migration-
  report.md, scaffold-tree-functions-migration-report.md) establish the chem-section pattern at this
  first cycle: `coverage_type: regression` for multi-format happy-path walks + SR-01 A-STRUCT-02
  carryforward + dataset-list extension per chain rev 2 directive. No "approaches off the table"
  constraints apply.
- **Cross-cutting bug citations (chain YAML `bug_focused_candidates[]`).** Per `migration-prompt.md`
  "Cross-cutting bug citations from chain YAML" §: this scenario IS a span in chain
  `bug_focused_candidates[chem-grok-18407-spec.ts]` (`spans: [chemical-space.md:Step 3,
  chemical-space.md:Step 5]`, `proposed_spec: chem-grok-18407-spec.ts`, rationale: "Chemical Space
  dialog silent-close on empty column input ... matched_steps = [3, 5] non-adjacent (gap=2 > 1) →
  Trigger 2 strict threshold met"). The dedicated bug-focused spec (`chem-grok-18407-spec.ts`)
  carries `coverage_type: edge` and exercises the empty-column-input invariant; this scenario carries
  `coverage_type: regression` and exercises the happy-path multi-format walk. The two specs are
  parallel-coverage on the same `chem.analyze.chemical-space.*` sub_features but assert orthogonal
  invariants. Per Pattern 2 discipline: the bug is NOT propagated into this scenario's
  `related_bugs` (it is the bug-focused spec's responsibility); the bug IS cited in Notes for
  cross-reference awareness. F-BUG-COVERAGE-01 at section-complete is authoritative for cross-
  cutting bug coverage.
- **UI delegation status.** Per chain YAML `ui_coverage_plan.smoke_scenario:
  Advanced/scaffold-tree-functions.md` (NOT this scenario) and chain
  `ui_coverage_plan.delegated_scenarios[] → chemical-space.md → delegated_to: null` ("UI coverage
  owned (Chemical Space editor + method selector + Cluster MCS toggle). Disjoint from Scaffold Tree
  smoke. Scenario owns its own coverage."), this scenario's `ui_coverage_responsibility` flows are
  owned in-scenario. No SCOPE_REDUCTION proposal substitutes JS API for any of the 4 UI flows;
  therefore the UI delegation discipline of migration-prompt § "UI delegation in SCOPE_REDUCTION
  proposals" does not apply (no SR proposal to delegate). Reaffirmed as a hard constraint in the
  migrated Notes ("No JS API substitution").
- **Scenario constraint extraction (per migration-prompt § Scenario constraint extraction).**
  - **(a) FORBIDDEN substitutions:** None explicit in the original body (no "must not", "forbidden",
    "not allowed" phrasing). HOWEVER `pyramid_layer: integration` per chain + `ui_coverage_delegated_
    to: null` (scenario owns its own coverage) imply JS API substitution for the 4
    `ui_coverage_responsibility` flows is NOT acceptable — surfaced explicitly in the migrated Notes
    (`No JS API substitution`) as a structured behavioral guard. In particular, direct invocation of
    `chemSpaceTransform` or `getChemSpaceEmbeddings` (atlas `chem.analyze.chemical-space.transform` /
    `chem.analyze.chemical-space.embeddings`) is NOT a substitute for the dialog-driven flow.
  - **(b) REQUIRED actions:** All 4 entries of `ui_coverage_responsibility` are exercised via UI
    driving in the migrated body (`chem-add-chemical-space` via top-menu walk,
    `chem-chemical-space-editor-dialog` via dialog open + interaction,
    `chem-chemical-space-method-selector` via the editor's method selector interaction at step 5,
    `chem-chemical-space-cluster-mcs-toggle` via the Cluster MCS toggle at step 5). No flow is
    deferred.
  - **(c) Missing-selector escalation:** Migrator does NOT have explicit selector definitions for
    the `ChemSpaceEditor` dialog method selector, Cluster MCS checkbox, or the embedding scatter-
    plot container. Spec-time selector discovery is the Automator's responsibility (per existing
    sibling `chemical-space-spec.ts` patterns + `page.evaluate` introspection); a reference-file
    approval-required proposal to add Chemical-Space-specific selectors to
    `grok-browser/references/viewers/chemical-space.md` is a Phase 2 deliverable — out of this
    migration's scope.
  - **(d) Reference templates:** Code-style anchors for the Automator (cited verbatim, not modified
    here): `public/packages/UsageAnalysis/files/TestTrack/Chem/chemical-space-spec.ts` (current
    existing test — sibling shape) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/activity-cliffs-spec.ts` (adjacent
    same-cycle sister — same analyze-dialog walk pattern, same `spec-login` / `page.evaluate` /
    `page.locator` patterns) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/calculate-spec.ts` (multi-format walk
    precedent for the dataset matrix loop).
  - **(e) Scenario authority clause:** No prompt-vs-scenario conflict surfaced. The original
    scenario is a 6-step linear dialog walk; the migrated body preserves all original steps + makes
    the implicit cross-cell isolation invariant explicit. No `prompt_scenario_conflict:` field to
    flag for orchestrator.
- **Source-text fixes silently applied (per chain rev 2 directive footer note (c) — Olena
  2026-05-11).**
  - **Dataset list extended from 1 entry to 3 entries.** Original JSON footer's `datasets[]` lists
    only `System:DemoFiles/chem/smiles.csv`, while body Step 1 prose says "should be tested on
    smiles, molV2000, molV3000 formats". Per chain rev 2 directive footer note (c) "body/JSON
    dataset mismatches in calculate / elemental-analysis / activity-cliffs / chemical-space /
    filter-panel / r-group-analysis / Advanced/scaffold-tree — fix during migration. Migrator MAY
    edit scenario body for typo correction + step renumbering + dataset list extension", the
    dataset list is extended to three entries here:
    - `System:DemoFiles/chem/smiles.csv` (Variant A — smiles, the original JSON footer's sole
      entry)
    - `System:AppData/Chem/mol1K.sdf` (Variant B — molV2000, per same-cycle sibling
      activity-cliffs.md + calculate.md + elemental-analysis.md canonical molV2000 source)
    - `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` (Variant C — molV3000, per same-cycle
      sibling canonical molV3000 source)
- **No invented sub_features / helpers.** Every id in `sub_features_covered` is present in
  `feature-atlas/chem.yaml` rev 2 (verified per atlas lines 294-335 for
  `chem.analyze.chemical-space.*` ids — 6 ids total: parent + top-menu + transform + embeddings +
  editor + similarity-matrix; this scenario covers 5 of 6 — `chem.analyze.chemical-space.similarity-
  matrix` is a programmatic helper not exposed via dialog walk and is intentionally omitted). Every
  reused helper is in `helpers-registry.yaml`; candidate helpers are flagged as candidates, NOT
  referenced by name in the migrated body.

## Opt-outs (SCOPE_REDUCTION proposals)

### SR-01: A-STRUCT-02 carryforward (chain-level edge/perf coverage)

**Cited technical dependency:** A-STRUCT-02 ("at least one scenario is edge case or negative path —
i.e. its frontmatter `coverage_type` is `edge` or `perf`") is a **section-level structural
invariant**. This scenario (`chemical-space.md`) is a multi-format dialog + transform walk —
`coverage_type: regression` is the natural fit per chain `pyramid_layer: integration` + the
multi-format × dialog-walk shape. Forcing `coverage_type: edge` would mis-classify (no specific
failure-mode invariant being asserted here — the empty-column-input silent-close invariant is owned
by the dedicated bug-focused spec `chem-grok-18407-spec.ts`, not by this happy-path migration).

**Resolution path:** A-STRUCT-02 satisfaction is chain-wide via 10 `bug_focused_candidates[]` in
`scenario-chains/chem.yaml` rev 2 (each chem-bug-focused spec naturally carries `coverage_type: edge`
once authored — `chem-grok-18407-spec.ts` directly targets `chem.analyze.chemical-space.*` and will
organically satisfy the section-level invariant in the chem-section bug-focused authoring phase).
Additionally `r-group-analysis.md` (same-cycle, same section, `coverage_type: edge`) already
organically satisfies the section-level invariant. Same SR-01 carried by sibling-cycle reports
`activity-cliffs-migration-report.md`, `calculate-migration-report.md`,
`scaffold-tree-functions-migration-report.md`, `scaffold-tree-migration-report.md`,
`info-panels-migration-report.md`, `Advanced/structure-filter-migration-report.md`,
`Advanced/similarity-search-migration-report.md`, `elemental-analysis-migration-report.md`.

**Scope of waiver:** A-STRUCT-02 verdict on `chemical-space.md` is deferred to chain-level
evaluation. Critic A should return SR (carryforward), not FAIL.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through A-STRUCT-06, A-COVERAGE-*, A-MERIT-*) are
expected to PASS on `chemical-space.md` without any other SR proposals. No step is opted out for
effort. The single entry above cites a real section-level structural property whose satisfaction path
is owned by the chain (and reinforced by the directly-targeting bug-focused spec
`chem-grok-18407-spec.ts`), not by this per-scenario migration (D-MERIT-01 compliant)._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist yet. The migrated body's 7 steps per
dataset variant are realizable against the current playwright + helpers-registry + atlas state. The
candidate helpers surfaced in Decisions are convenience abstractions that the Automator may inline
pending registration — they do NOT block spec realization. The selector-reference write proposal
(Decisions § Scenario constraint extraction (c)) is a downstream Phase 2 enhancement, not a blocker
on this scenario's automation (existing `chemical-space-spec.ts` already drives the Chemical Space
dialog DOM via `page.evaluate` + `page.locator` per the sibling-test pattern; Automator can extend
that pattern without a reference update)._

## Edge cases

The original lists no explicit "edge case" keyword. Implicit edge cases derivable from the scenario
content:

- **Multi-format input variant (smiles / molV2000 / molV3000).** Original Step 1 enumerates three
  format variants in body prose; the JSON footer only carried one entry. PRESERVED as scenario steps
  (the D1-D3 dataset matrix — every variant is exercised as a separate cell in the per-dataset
  loop). The dataset-list extension is the silent fix per chain rev 2 directive.
- **Default-parameters vs custom-parameters pass.** Original Step 3 (default) and Steps 4-6 (custom)
  test the editor dialog under two distinct parameter configurations. PRESERVED as scenario steps
  (Scenarios steps 3 and 4-6 — explicit two-pass invariant). Implicit edge case: re-running Chemical
  Space on the same active table produces a fresh embedding scatter; the existing
  `chem-space-embedding-col` columns from the previous run may be overwritten or appended (per atlas
  `chem.analyze.chemical-space.transform` — `tags result columns`). The migrated body asserts
  "fresh embedding scatter plot is produced that reflects the edited parameters" verification at
  step 6.
- **Method-selector branch (UMAP vs t-SNE).** Atlas `chem.analyze.chemical-space.embeddings` exposes
  `methodName` as a parameter; the editor surfaces UMAP and t-SNE as the two selectable methods.
  Migrated body Step 5 enumerates "switch the method selector between UMAP and t-SNE" as one of the
  parameter-change options. Implicit: each method produces a different embedding for the same
  similarity matrix. Spec-time the Automator may pick a deterministic switch (e.g. default-pass
  UMAP, custom-pass t-SNE) for reproducible verification.
- **Cluster MCS toggle.** Atlas `chem.analyze.chemical-space.editor` notes the Cluster MCS toggle is
  exposed on top of the standard dim-reduction editor; toggling it ON adds DBSCAN clustering + an
  MCS overlay producing the additional `chem-space-cluster-col` tagged column. Migrated body Step 5
  enumerates the toggle as one of the parameter-change options; Step 6 verifies the conditional
  `chem-space-cluster-col` column append.
- **Empty-column-input silent close (GROK-18407, atlas `chem.cp.chemical-space` p1 edge case).**
  Original body does NOT exercise this edge case (per chain rev 2
  `dependency_graph.chemical-space.md.notes`: "Body does NOT test empty-column-input edge case (the
  GROK-18407 silent-close invariant) — this gap is captured under bug_focused_candidates[]"). NOT
  PRESERVED here — owned by bug-focused spec `chem-grok-18407-spec.ts` per chain
  `bug_focused_candidates[]`. Cited in Notes for cross-reference awareness; NOT propagated into
  `related_bugs: []` (Pattern 2 discipline — the bug-focused spec is the authoritative coverage
  owner).
- **Row-count / dataset-volume behavior.** Atlas `chem.analyze.chemical-space.top-menu` does NOT
  document an explicit row-count cap or long-running-confirm threshold for Chemical Space (unlike
  Activity Cliffs which has explicit 10,000-row cap + 5000/2000-row confirm). None of the three
  linked datasets exceed any plausible threshold (smiles.csv ≈ 700 rows; mol1K.sdf = 1000 rows;
  ApprovedDrugs2015.sdf ≈ 1900 rows). PRESERVED implicitly — the happy-path walk on these datasets
  exercises the typical analyst workload; any large-volume cap invariant would warrant a dedicated
  bug-library entry, not propagation here.
- **Console errors throughout.** Implicit across all 6 original steps: each action ("Open", "Run top
  menu", "Click OK", re-run, change parameters, OK) must complete without console errors. PRESERVED
  as per-step verification weave (Scenarios steps 3 and 6 explicitly include "no console errors"
  assertions). Sister with `calculate.md` and `activity-cliffs.md` patterns.

No edge case is moved to atlas, manual_only, deferred, or a separate scenario. All implicit edge
cases are PRESERVED as scenario steps OR Notes-level caveats in the migrated body, except the
empty-column-input invariant which is owned by the parallel bug-focused spec per chain Pattern 2.

## Unresolved ambiguities

- **"Change the parameters arbitrarily" — which specific parameter to change.** Original Step 5
  asks to "change the parameters arbitrarily" — silent on which parameter(s), or how many. The
  Chemical Space editor (`ChemSpaceEditor`, atlas `chem.analyze.chemical-space.editor`) exposes
  the molecule column picker + method selector (UMAP / t-SNE) + similarity metric + Cluster MCS
  toggle. Migrated body resolves by asserting "change at least one parameter from its default value
  (e.g. switch the method selector between UMAP and t-SNE, change the similarity metric, or toggle
  Cluster MCS on / off)" — any deterministic change is acceptable. Automator at spec time may pick
  a single high-signal parameter (e.g. method UMAP→t-SNE, or Cluster MCS off→on which adds a new
  column observably) for reproducible verification. Flag for QA pair review at spec realization.
- **Method selector default value.** Atlas `chem.analyze.chemical-space.embeddings` accepts
  `methodName` as a parameter but does not document the editor's default selection. The default may
  be UMAP or t-SNE depending on the editor's initialization logic. Migrated body does not depend on
  the default value — the custom-parameters pass at Step 5 picks an arbitrary other method.
  Implementation-defined; spec-time the Automator may read the default at dialog-open time and
  switch to the other.
- **"Fresh scatter plot" — overwrite vs append semantics.** Original Step 6 implies that re-running
  Chemical Space with edited parameters produces a fresh result. The `chemSpaceTransform` per atlas
  `chem.analyze.chemical-space.transform` "tags result columns with `chem-space-embedding-col` /
  `chem-space-cluster-col`"; whether the second run overwrites the existing tagged columns or
  appends new ones (e.g. `chem-space-embedding-col-2`) is implementation-defined. The migrated body
  asserts the scatter plot reflects edited parameters (visible distribution difference) without
  prescribing the column-tagging mechanism. Flag for QA pair review at spec realization.
- **Cluster MCS toggle visibility on the editor.** Atlas `chem.analyze.chemical-space.editor`
  describes the Cluster MCS toggle but does not document its DOM-selector format. Spec-time
  selector strategy is the Automator's responsibility (per existing sibling `chemical-space-spec.ts`
  patterns + dialog-element introspection). Flag for QA pair review.
