# Migration Report — activity-cliffs.md

First-cycle migration of `activity-cliffs.md` (TestTrack/Chem section) per chain YAML
`scenario-chains/chem.yaml` rev 2 (`output_plan.activity-cliffs.md`: `target_layer: playwright`,
strategy `simple`, classification `medium`, `pyramid_layer: integration`, `depends_on: []`).
Multi-format Activity Cliffs viewer + cliffs-grid + property-panel interaction walk across SMILES /
molV2000 / molV3000 / 2-column SMILES variants, plus the SPGI reference dataset declared in the
original's JSON footer.

## Step mapping

The original is a 23-line scenario with 9 numbered steps (gap between step 7 and step 9 — no step 8 in
the original body; chain rev 2 footer note (c) "source-text defects" — fix during migration) plus a
JSON footer with `order: 7` and 5 datasets. The body lists "smiles (smiles), molV2000 (mol1K) and
molV3000 (ApprovedDrugs2015)" as linked datasets at step 1 with a TODO bullet for
`smiles_2_columns.csv`, while the JSON footer lists those four plus `System:DemoFiles/SPGI.csv` as a
fifth dataset; the JSON footer's `smiles_2_columns.csv` path is at
`System:AppData/UsageAnalysis/test_datasets/smiles_2_columns.csv` while the file actually lives at
`System:AppData/Chem/tests/smiles_2_columns.csv` (verified per repo search at
`public/packages/Chem/files/tests/smiles_2_columns.csv`, per the calculate.md migration finding).
Defects fixed silently per chain rev 2 directive (Olena 2026-05-11, chain footer note (c) "source-text
defects in scenario bodies — Migrator fixes during migration").

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Step 1 — "Open: Linked datasets (should be tested on smiles (smiles), molV2000 (mol1K) and molV3000 (ApprovedDrugs2015) formats); smiles_2_columns.csv (TODO: add to linked datasets)" | Setup step 1 (provisioning — datasets D1 / D2 / D3 / D4 / D5 enumerated, with `smiles_2_columns.csv` TODO resolved at the verified path) + Scenarios > "Activity Cliffs viewer walk per dataset" step 1 (dataset open as the per-cell loop entry) | preserved (split for clarity — provisioning declaration in Setup; the open action in scenario body; TODO resolved + JSON dataset path corrected) |
| Step 2 — "Run **Chem > Analyze > Activity cliffs**. A dialog opens." | Scenarios > "Activity Cliffs viewer walk per dataset" step 2 (top-menu walk → editor dialog open + column-picker selection for Variant D) | preserved (split for clarity — Variant D column-selection branch added per dialog-shape variance) |
| Step 3 — "Click OK to run a function with default parameters" | Scenarios > "Activity Cliffs viewer walk per dataset" step 3 (OK click with defaults; verify scatter + cliffs overlay rendered; `activityCliffsParams` tag set; no console errors) | preserved (with expected-result verification made explicit per D-STEP-02) |
| Step 4 — "Check 'Show only cliffs' - only points with cliffs should be left on scatter plot" | Scenarios > "Activity Cliffs viewer walk per dataset" step 4 (Show-only-cliffs toggle + verify only cliff-participating points remain visible) | preserved (action + verification kept together; phrasing kept verbatim where possible) |
| Step 5 — "Click on the link with number of cliffs - grid with cliffs should appear on the bottom" | Scenarios > "Activity Cliffs viewer walk per dataset" step 5 (cliff-count link click + verify cliffs grid appears docked below + populated with cliff pair rows) | preserved (action + verification kept together) |
| Step 6 — "Click on the first row of the grid - the scatter plot is zoomed to the cliff, property panel is shown with the pair of molecules (highlighting uncommon parts), activity difference. When hovering over line the tooltip with the same information should be showed." | Scenarios > "Activity Cliffs viewer walk per dataset" step 6 (grid row click + scatter zoom + property panel pair + uncommon-parts highlight + activity difference) + step 7 (cliff-line hover tooltip with same information) | preserved (split for clarity — grid-row-click action / verification block and the line-hover tooltip verification block separated; both halves preserved verbatim) |
| Step 7 — "Double click to unzoom the scatter plot. And zoom any cluster to see any line. Click on the line - the scatter plot will be zoomed. And corresponding line will become current in the grid below." | Scenarios > "Activity Cliffs viewer walk per dataset" step 8 (double-click unzoom + cluster zoom + line click + verify scatter zoom + corresponding grid row becomes current) | preserved (composite step expanded into a single migrated step covering the entire unzoom → re-zoom → line-click → grid-current sequence) |
| (gap — original has no Step 8; Step 7 → Step 9 directly) | n/a | source-text defect (silently fixed — migrated body numbers steps cleanly without gap; chain rev 2 footer note (c) "source-text defects" — Migrator MAY edit for step renumbering) |
| Step 9 (original numbering) — "Run **Activity cliffs** one more time." | Scenarios > "Activity Cliffs viewer walk per dataset" step 9 (re-run Activity Cliffs top-menu) | preserved (renumbered cleanly) |
| Step 10 (original numbering) — "Change the parameters arbitrarily" | Scenarios > "Activity Cliffs viewer walk per dataset" step 10 (change at least one editor dialog parameter — similarity threshold / method / metric / preprocessing) | preserved (renumbered cleanly; arbitrary-change made concrete via the enumerated parameter list) |
| Step 11 (original numbering) — "Click OK to run a function with edited parameters." | Scenarios > "Activity Cliffs viewer walk per dataset" step 11 (OK click with edited params; verify fresh scatter + cliffs overlay reflects edited parameters; no console errors) | preserved (renumbered cleanly; expected-result verification added per D-STEP-02 — fresh scatter and / or differing cliff distribution from step 3) |
| (implicit — close active view before next dataset) | Scenarios > "Activity Cliffs viewer walk per dataset" step 12 (close active view between cells) | preserved as cleanup (implicit per-cell isolation invariant made explicit) |
| JSON footer `{ "order": 7, "datasets": ["System:DemoFiles/SPGI.csv", "System:DemoFiles/chem/smiles.csv", "System:AppData/Chem/mol1K.sdf", "System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf", "System:AppData/UsageAnalysis/test_datasets/smiles_2_columns.csv"] }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` + migrated Setup step 1 dataset enumeration — extended to 5 entries including the resolved `smiles_2_columns.csv` path at `System:AppData/Chem/tests/smiles_2_columns.csv` and the SPGI reference dataset) |

No original step is silently dropped. Step-numbering gap (original `7.` → `9.` with no `8.`) renumbered
cleanly into the migrated body's 12-step linear sequence.

## Decisions

- **Why this `target_layer`:** chose `playwright` per `scenario-chains/chem.yaml`
  `output_plan.activity-cliffs.md.target_layer = playwright`. The scenario requires DOM-level UI
  driving: top-menu walk (**Chem | Analyze | Activity Cliffs**), Activity Cliffs editor dialog
  interaction, Show-only-cliffs toggle, cliff-count link click, cliffs-grid row click, scatter-plot
  double-click unzoom + cluster zoom + line click + line-hover tooltip, and property-panel pair-of-
  molecules display with uncommon-parts highlighting — all UI behaviors that Playwright drives
  natively and the apitest layer cannot exercise as an interaction walk. Chain YAML reason field:
  "Activity Cliffs viewer + grid interaction walk across 4 format variants. Sibling
  activity-cliffs-spec.ts already exists at playwright per existing-test-index. Scatter zoom +
  cliff-count link + grid-row click + property-panel pair-display require real DOM interaction.
  Single scenario; simple strategy."
- **Why this `coverage_type`:** chose `regression` per chain rev 2 (`pyramid_layer: integration` —
  multi-format × interaction-walk; not `smoke` since the section's smoke is
  `Advanced/scaffold-tree-functions.md` per chain `ui_coverage_plan.smoke_scenario`; not `edge` /
  `perf` since the scenario walks the happy path only — no specific failure-mode invariant being
  asserted, no row-count cap or long-running-confirm threshold exercised at the given dataset volumes).
  The dialog-level empty-input / invalid-input failure-mode invariant is owned by sibling
  `chemical-space.md` + bug-focused spec `chem-grok-18407-spec.ts` (`GROK-18407` targets
  `chem.analyze.chemical-space.*`, not `chem.analyze.activity-cliffs.*` — parallel-coverage, different
  sub_feature scope). `regression` is the natural fit for a multi-format × multi-interaction walk that
  guards against regression-of-the-set across the Activity Cliffs viewer surface.
- **Why this `strategy`:** `simple` per chain YAML `output_plan.activity-cliffs.md.strategy = simple`.
  Single scenario, no cross-file fixture; chain analyzer classified `medium` (multi-format /
  multi-interaction) but kept strategy `simple` — single migrated `.md`, single downstream spec.
  Pattern 1 (decomposition) is NOT triggered (chain classifier returned `medium`, not
  `complex-standalone`; no child scenarios with `produced_from: decomposed`). Pattern 2 (bug-focused
  slice) does NOT apply — `activity-cliffs.md` is NOT a span in any of the chain's 10
  `bug_focused_candidates[]`; no GROK ticket invariant intersects `chem.analyze.activity-cliffs.*`.
  Pattern 3 (UI-smoke consolidation) is NOT triggered (no `ui_consolidation_proposals[]` from Gate F).
  Standard per-scenario migration per "Migrated scenario shape".
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/activity-cliffs-spec.ts` — existing
    playwright-layer test ("Chem: Activity Cliffs", category Chem, `features_covered:
    [chem.activity-cliffs]`, helpers `spec-login`, patterns `uses-fixture-runners` / `uses-grok.dapi`
    / `uses-grok.shell` / `uses-page.evaluate` / `uses-page.locator` / `uses-playwright-test`). Per
    `existing-test-index.yaml` line 32321-32338. The Automator will extend it to cover the full
    enumerated dataset variant matrix + interaction walk per this migrated scenario. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/calculate.md` +
    `calculate-migration-report.md` (2026-05-11 predecessor in same migration cycle) — shape anchor
    for the migrated `.md` (multi-format dataset matrix + per-cell interaction loop + `coverage_type:
    regression` precedent + JSON dataset path correction + step-numbering fix precedent for
    `smiles_2_columns.csv` TODO resolution). Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts` — sibling Chem
    playwright test ("Chem: Info Panels on smiles.csv"). Same shape pattern; confirms the section-wide
    playwright convention. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/scaffold-tree-functions.md` +
    `scaffold-tree-functions-migration-report.md` (2026-05-11 predecessor in same migration cycle) —
    SR-01 A-STRUCT-02 carryforward template adopted here. Read-only.
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`) —
      section-standard fixture login; Automator will use in spec `beforeEach`.
    - `softStep` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`) — per-step
      soft-fail wrapper; Automator wraps each cell verification in `softStep` for granular failure
      reporting (especially important for a 5-variant × 12-step interaction walk where individual
      cells may surface format-specific dialog branches).
    - `closeAllViews` (`public/packages/UITests/playwright/scripts/helpers.ts:129`) — cleanup at the
      end of each dataset cell + at the end of the scenario.
    - `expandAccordionPane` (`public/packages/UITests/playwright/scripts/helpers.ts:141`) —
      applicable in step 6 (expand property-panel accordion groups when the pair-of-molecules +
      activity-difference info is grouped under collapsed accordion sections).
  - **Candidate helpers (NOT yet in registry — flagged for addition via helpers-registry curator;
    per migration-prompt Helpers discipline §):**
    - `helpers.playwright.chem.addActivityCliffsViewer(page, dataset, params?)` — opens dataset +
      walks **Chem | Analyze | Activity Cliffs** top-menu + accepts dialog with default or supplied
      params; returns a locator for the Activity Cliffs scatter plot. Surfaced by Setup step 1 +
      Scenarios steps 2-3. Convenience abstraction; spec Author can inline the top-menu walk until
      the helper lands.
    - `helpers.playwright.chem.toggleActivityCliffsShowOnlyCliffs(viewer)` — toggles the
      Show-only-cliffs property on the Activity Cliffs scatter plot. Surfaced by Scenarios step 4.
    - `helpers.playwright.chem.clickActivityCliffsCliffCountLink(viewer)` — clicks the cliff-count
      link to open the cliffs grid below. Surfaced by Scenarios step 5.
    - `helpers.playwright.chem.clickActivityCliffsGridRow(grid, rowIndex)` — clicks the Nth row in
      the cliffs grid. Surfaced by Scenarios step 6.
- **Bug library consulted:** yes — `bug-library/chem.yaml` revision 1 (12 curated_bugs, generated
  2026-05-05). No curated bug's `affects` intersects this scenario's `sub_features_covered`
  (`chem.analyze.activity-cliffs.*`). The 12 curated bugs target Scaffold Tree / Sketcher /
  Substructure Search / Filter Panel / Descriptors / MMP / R-Groups / Chemical Space / Notations /
  CSV Export / MPO Profile / RDKit cell renderer — none target Activity Cliffs. `related_bugs: []`
  in the migrated frontmatter. Closest neighbor — `GROK-18407` (Chemical Space empty-input silent
  dialog close) — targets `chem.analyze.chemical-space.*` exclusively, NOT
  `chem.analyze.activity-cliffs.*`; it is a sister-pattern (analyze dialog input-validation invariant)
  but a different feature scope. The Activity Cliffs dialog has its own row-count cap (10,000 rows)
  and long-running-confirm threshold (5000 / 2000 rows depending on method) per atlas
  `chem.analyze.activity-cliffs.top-menu`; these are NOT exercised here (dataset volumes below
  threshold) and any future failure-mode invariant for the Activity Cliffs editor would warrant a
  dedicated bug-library entry, not propagation of `GROK-18407`. No related-bug citation in
  frontmatter; awareness only here.
- **Decision log queried:** yes — `decision-log.yaml` grepped for `feature: chem` — **zero matches**.
  This is the FIRST CYCLE for the chem section. No prior `migration_decisions`, `layer_decisions`,
  `manual_only`, or `failed_attempts` entries exist for chem. No "approaches off the table"
  constraints apply.
- **Cross-cutting bug citations (chain YAML `bug_focused_candidates[]`).** Per `migration-prompt.md`
  "Cross-cutting bug citations from chain YAML" §: this scenario is NOT a span in any of the 10
  `bug_focused_candidates[]` (`chem-grok-12758-spec.ts` / `chem-grok-14028-spec.ts` /
  `chem-grok-17595-spec.ts` / `chem-grok-17621-spec.ts` / `chem-grok-17964-spec.ts` /
  `chem-grok-18407-spec.ts` / `chem-grok-19624-spec.ts` / `chem-github-2942-spec.ts` /
  `chem-github-3004-spec.ts` / `chem-grok-16870-spec.ts`). No cross-cutting bug-focused citation
  applies — awareness only that `GROK-18407` is a sister-pattern analyze-dialog input-validation
  invariant owned by `chemical-space.md`. F-BUG-COVERAGE-01 at section-complete is authoritative for
  cross-cutting bug coverage.
- **UI delegation status.** Per chain YAML `ui_coverage_plan.smoke_scenario:
  Advanced/scaffold-tree-functions.md` (NOT this scenario) and chain
  `ui_coverage_plan.delegated_scenarios[] → activity-cliffs.md → delegated_to: null` ("UI coverage
  owned (Activity Cliffs editor dialog + viewer property-panel + Show-only-cliffs toggle +
  cliff-count link + grid row click + scatter zoom + line-click interaction). Disjoint from Scaffold
  Tree smoke. Scenario owns its own coverage."), this scenario's `ui_coverage_responsibility` flows
  are owned in-scenario. No SCOPE_REDUCTION proposal substitutes JS API for any of the 7 UI flows;
  therefore the UI delegation discipline of migration-prompt § "UI delegation in SCOPE_REDUCTION
  proposals" does not apply (no SR proposal to delegate). Reaffirmed as a hard constraint in the
  migrated Notes ("No JS API substitution").
- **Scenario constraint extraction (per migration-prompt § Scenario constraint extraction).**
  - **(a) FORBIDDEN substitutions:** None explicit in the original body (no "must not", "forbidden",
    "not allowed" phrasing). HOWEVER `pyramid_layer: integration` per chain + `ui_coverage_delegated_
    to: null` (scenario owns its own coverage) imply JS API substitution for the 7
    `ui_coverage_responsibility` flows is NOT acceptable — surfaced explicitly in the migrated Notes
    (`No JS API substitution`) as a structured behavioral guard.
  - **(b) REQUIRED actions:** All 7 entries of `ui_coverage_responsibility` are exercised via UI
    driving in the migrated body (add via top menu, editor dialog interaction, Show-only-cliffs
    toggle, cliff-count link click, cliffs-grid row click, scatter-plot zoom + line click,
    property-panel pair-of-molecules display). No flow is deferred.
  - **(c) Missing-selector escalation:** Migrator does NOT have explicit selector definitions for
    the Activity Cliffs scatter-plot viewer container, the Show-only-cliffs property toggle, the
    cliff-count link, the cliffs grid root, or the property-panel pair-of-molecules block. Spec-time
    selector discovery is the Automator's responsibility (per existing sibling
    `activity-cliffs-spec.ts` patterns + `page.evaluate` introspection); a reference-file
    approval-required proposal to add Activity-Cliffs-specific selectors to
    `grok-browser/references/viewers/activity-cliffs.md` is a Phase 2 deliverable — out of this
    migration's scope.
  - **(d) Reference templates:** Code-style anchors for the Automator (cited verbatim, not modified
    here): `public/packages/UsageAnalysis/files/TestTrack/Chem/activity-cliffs-spec.ts` (current
    existing test — sibling shape) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts` (section-mate
    playwright walk; `spec-login` + `page.evaluate` + `page.locator` patterns) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/calculate-spec.ts` (multi-format walk
    precedent for the dataset matrix loop).
  - **(e) Scenario authority clause:** No prompt-vs-scenario conflict surfaced. The original
    scenario is an 11-step (with gap) interaction walk; the migrated body preserves all original
    steps + makes the implicit cross-cell isolation invariant explicit. No `prompt_scenario_conflict:`
    field to flag for orchestrator.
- **Source-text fixes silently applied (per chain rev 2 directive footer note (c) — Olena
  2026-05-11).**
  - **Step-numbering gap closed.** Original numbered steps `1 → 2 → 3 → 4 → 5 → 6 → 7 → 9 → 10 → 11`
    (no step 8). Renumbered cleanly in the migrated body as 12 contiguous steps; no semantic content
    moved or dropped. The gap is documented in the step-mapping table as `(gap — original has no
    Step 8; Step 7 → Step 9 directly) → source-text defect (silently fixed)` per D-SAN-02.
  - **`smiles_2_columns.csv` TODO resolved.** Original step 1 sub-bullet flags
    `smiles_2_columns.csv (TODO: add to linked datasets)`; the JSON footer's path is
    `System:AppData/UsageAnalysis/test_datasets/smiles_2_columns.csv` but the file actually lives at
    `System:AppData/Chem/tests/smiles_2_columns.csv` (verified per repo search at
    `public/packages/Chem/files/tests/smiles_2_columns.csv`, per the `calculate.md` migration finding
    which resolved the same TODO in the same cycle). Migrated body uses the canonical Chem package
    path in Setup step 1 (Variant D); the original JSON footer's incorrect UsageAnalysis path is
    superseded.
  - **SPGI dataset reconciled.** Original step 1 body prose lists only smiles / molV2000 / molV3000
    / 2-column SMILES (4 formats), while the JSON footer's `datasets[]` lists 5 entries including
    `System:DemoFiles/SPGI.csv`. Chain rev 2 dependency_graph notes:
    `notes: "Datasets: SPGI.csv, smiles.csv, mol1K.sdf, ApprovedDrugs2015.sdf, smiles_2_columns.csv
    (the JSON includes SPGI but step 1 lists only the latter four formats — flagged under
    unresolved_ambiguities)"`. Per chain rev 2 directive footer note (c) "dataset list extension"
    explicit for activity-cliffs.md, the SPGI entry is preserved in Setup step 1 as Variant E
    (reference dataset — completeness, not a new format-variant — molV2000 is already covered by
    `mol1K.sdf`).
- **No invented sub_features / helpers.** Every id in `sub_features_covered` is present in
  `feature-atlas/chem.yaml` rev 2 (verified per atlas lines 337-372 for
  `chem.analyze.activity-cliffs.*` ids — 5 ids total: parent + top-menu + transform + init + editor).
  Every reused helper is in `helpers-registry.yaml`; candidate helpers are flagged as candidates, NOT
  referenced by name in the migrated body.

## Opt-outs (SCOPE_REDUCTION proposals)

### SR-01: A-STRUCT-02 carryforward (chain-level edge/perf coverage)

**Cited technical dependency:** A-STRUCT-02 ("at least one scenario is edge case or negative path —
i.e. its frontmatter `coverage_type` is `edge` or `perf`") is a **section-level structural
invariant**. This scenario (`activity-cliffs.md`) is a multi-format viewer + interaction walk —
`coverage_type: regression` is the natural fit per chain `pyramid_layer: integration` + the
multi-format × multi-interaction shape. Forcing `coverage_type: edge` would mis-classify (no specific
failure-mode invariant, no threshold being asserted).

**Resolution path:** A-STRUCT-02 satisfaction is chain-wide via 10 `bug_focused_candidates[]` in
`scenario-chains/chem.yaml` rev 2 (each chem-bug-focused spec naturally carries `coverage_type: edge`
once authored — sister-pattern to predecessor migrations in the same cycle:
`Advanced/scaffold-tree-functions-migration-report.md` SR-01,
`r-group-analysis-migration-report.md` SR-01, `calculate-migration-report.md` SR-01,
`elemental-analysis-migration-report.md` SR-01, `sketcher-migration-report.md` SR-01).

**Scope of waiver:** A-STRUCT-02 verdict on `activity-cliffs.md` is deferred to chain-level
evaluation. Critic A should return SR (carryforward), not FAIL.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through A-STRUCT-06, A-COVERAGE-*, A-MERIT-*) are
expected to PASS on `activity-cliffs.md` without any other SR proposals. No step is opted out for
effort. The single entry above cites a real section-level structural property whose satisfaction path
is owned by the chain, not by this per-scenario migration (D-MERIT-01 compliant)._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist yet. The migrated body's 12 steps per
dataset variant are realizable against the current playwright + helpers-registry + atlas state. The
candidate helpers surfaced in Decisions are convenience abstractions that the Automator may inline
pending registration — they do NOT block spec realization. The selector-reference write proposal
(Decisions § Scenario constraint extraction (c)) is a downstream Phase 2 enhancement, not a blocker
on this scenario's automation (existing `activity-cliffs-spec.ts` already drives the Activity Cliffs
viewer DOM via `page.evaluate` + `page.locator` per the sibling-test pattern; Automator can extend
that pattern without a reference update)._

## Edge cases

The original lists no explicit "edge case" keyword. Implicit edge cases derivable from the scenario
content:

- **Multi-format input variant (smiles / molV2000 / molV3000 / 2-column SMILES / SPGI).** Original
  Step 1 enumerates four format variants in body prose plus a fifth (SPGI) in the JSON footer.
  PRESERVED as scenario steps (Scenarios step 1 + the D1-D5 dataset matrix — every variant is
  exercised as a separate cell in the per-dataset loop).
- **Two-molecule-column variant (`smiles_2_columns.csv`).** Implicit dialog-shape variance: the
  Activity Cliffs editor's column-picker behavior differs when the table has two molecule columns
  vs. one (auto-detect picks the first; user can override). PRESERVED as scenario step (Scenarios
  step 2 — explicit "for Variant D (two molecule columns) select the first molecule column in the
  column picker" guidance).
- **Default-parameters vs custom-parameters pass.** Original Steps 3 (default) and Steps 9-11
  (custom) test the editor dialog under two distinct parameter configurations. PRESERVED as scenario
  steps (Scenarios steps 3 and 9-11 — explicit two-pass invariant). Implicit edge case: re-running
  Activity Cliffs on the same active table produces a fresh scatter plot with potentially different
  cliff distribution; the `activityCliffsParams` table tag is overwritten on re-run. The migrated
  body asserts the "fresh scatter + cliffs overlay reflects edited parameters" verification at step
  11.
- **Cliff count == 0 silently / non-zero.** The cliff-count link's behavior when zero cliffs are
  detected at the given similarity threshold is implicit — the original does not specify. The
  migrated body assumes non-zero cliff count on the default-parameters pass (Steps 3-8 walk the
  cliffs grid + row click — pre-condition is that cliffs exist). Atlas / future regression depth may
  want a "zero-cliff threshold" edge case (e.g. set similarity threshold to 1.0 or 0.0 and verify
  graceful empty-state); OUT OF SCOPE for this regression walk. Atlas curator note (no atlas write
  here): an atlas-level `edge_cases` entry for the zero-cliff threshold behavior could be a
  candidate addition.
- **Row-count cap (10,000 rows) and long-running-confirm threshold (5000 / 2000 rows depending on
  method).** Atlas `chem.analyze.activity-cliffs.top-menu`: the Activity Cliffs run is bounded at
  10,000 rows; a long-running-analysis confirm dialog appears above 5000 / 2000 rows depending on
  method. None of the five linked datasets exceed these thresholds (smiles.csv ≈ 700 rows; mol1K.sdf
  = 1000 rows; ApprovedDrugs2015.sdf ≈ 1900 rows; smiles_2_columns.csv small; SPGI.csv ≈ 1200 rows).
  PRESERVED as Notes-level caveat ("the bounded-row dialog confirm... is NOT exercised because the
  linked datasets sit below the threshold; the cap invariant is observable only on a large-volume
  dataset and is intentionally out of scope for this regression walk"). Flagged for atlas curator
  consideration as a candidate `edge_cases` / `manual_only` entry — out of scope here.
- **Console errors throughout.** Implicit across all 11 original steps: each action ("Open", "Run
  top menu", "Click OK", toggle, "Click on the link", grid row click, line hover, double-click,
  zoom, line click) must complete without console errors. PRESERVED as per-step verification weave
  (Scenarios steps 3, 5, 6, 11 explicitly include "no console errors" assertions). Sister with
  `calculate.md` pattern.

No edge case is moved to atlas, manual_only, deferred, or a separate scenario. All implicit edge
cases are PRESERVED as scenario steps OR Notes-level caveats in the migrated body, with the
row-count cap flagged for atlas curator consideration.

## Unresolved ambiguities

- **"Click on the link with number of cliffs" — exact DOM selector and label format.** Original Step
  5 references a "link with number of cliffs" — implementation surface is in
  `Chem/src/package.ts` (atlas `chem.analyze.activity-cliffs.init`) and may render as a text link in
  the scatter-plot viewer overlay or in the property panel. The exact label format ("N cliffs", "N
  cliff pairs", or similar) and the locator strategy are implementation-defined. Migrated body
  describes the action semantically; Automator at spec time may use a text-match selector
  (`text=/\d+ cliff/i`) or a structural locator based on the scatter-plot viewer's child structure.
  Flag for QA pair review at spec realization.
- **"First row of the grid" — first by DOM-render order or first by cliff-magnitude / pair-index.**
  Original Step 6 says "click on the first row of the grid" without specifying ordering. The cliffs
  grid presents cliff pairs in some default order (likely by cliff magnitude descending, or by pair
  index ascending). Migrated body preserves "first row of the cliffs grid" — DOM-order interpretation.
  Automator at spec time picks a deterministic selector (positional `nth-child(1)` is acceptable; if
  the grid's sort order is stable across runs the DOM-first row IS the semantically-first cliff).
  Flag for QA pair review.
- **"Hovering over line" — line locator in the zoomed scatter plot.** Original Step 6's tooltip
  verification ("when hovering over line the tooltip... should be showed") depends on the cliff
  edges being rendered as discrete SVG / canvas line elements that can receive `mouseover`. The
  Activity Cliffs scatter plot draws cliffs as overlaid lines (`activityCliffsInitFunction`); the
  exact hit-test mechanism (SVG line element vs. canvas with hover hit-detection) is
  implementation-defined. Spec-time hover strategy may use `page.hover()` on a line selector if SVG,
  or `page.mouse.move()` to canvas coordinates if canvas. Flag for QA pair review.
- **"Change the parameters arbitrarily" — which specific parameter to change.** Original Step 10
  asks to "change the parameters arbitrarily" — silent on which parameter(s), or how many. The
  Activity Cliffs editor (`ActivityCliffsEditor`, atlas `chem.analyze.activity-cliffs.editor`)
  exposes similarity threshold + method + metric + preprocessing + options. Migrated body resolves
  by asserting "change at least one parameter from its default value (e.g. similarity threshold,
  method, metric, or preprocessing option)" — any deterministic change is acceptable. Automator at
  spec time may pick a single high-signal parameter (e.g. lower similarity threshold from default
  0.8 to 0.5 → more cliffs detected) for reproducible verification. Flag for QA pair review.
- **"Double click to unzoom" — exact zoom-reset target.** Original Step 7 says "double click to
  unzoom the scatter plot" — implicit that the scatter plot's double-click handler resets the zoom
  to the full domain. The exact double-click target (anywhere on the scatter plot canvas, or
  specifically on the background) is implementation-defined. Migrated body says "Double-click the
  scatter plot to unzoom"; Automator at spec time may double-click on the scatter-plot center or on
  a deliberately empty background area. Flag for QA pair review.
- **Cliffs grid docking position ("on the bottom").** Original Step 5 says "grid with cliffs should
  appear on the bottom" — implicit that the cliffs grid docks below the scatter plot in the same
  view container. The docking surface is implementation-defined (could be a sibling
  TableView, a dock-spawn pane below the scatter plot, or an attached property-panel block).
  Migrated body says "docked below the scatter plot" without specifying the exact dock surface;
  Automator at spec time resolves via the actual implementation. Flag for QA pair review.
