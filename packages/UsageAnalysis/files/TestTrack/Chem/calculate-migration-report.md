# Migration Report — calculate.md

First-cycle migration of `calculate.md` (TestTrack/Chem section) per chain YAML
`scenario-chains/chem.yaml` rev 2 (`output_plan.calculate.md`:
`target_layer: playwright`, strategy `simple`, classification `medium`,
`pyramid_layer: integration`, `depends_on: []`). Multi-format × multi-Calculate-menu-item
happy-path walk; realizes atlas critical path `chem.cp.calculate-descriptors-docker` (p1)
for the Descriptors anchor.

## Step mapping

The original is a 17-line, 4-step (with collision) scenario plus a prose "Do the asme [sic]
for each section in the Calculate menu ()" enumeration directive plus a JSON footer with
`order: 3` and 3 datasets. Body defects per chain rev 2 directive (Olena 2026-05-11,
chain footer note (i) "Source-text defects in scenario bodies — Migrator fixes during
migration. Affected: calculate.md (typos 'asme', 'menu ()')") fixed silently.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Step 1 (first) — "Open: Linked datasets (check for smiles (smiles.csv), molV3000 (ApprovedDrugs2015.sdf), molV2000 (mol1K.sdf) formats)" | Setup step 1 (provisioning — datasets D1, D2, D3 enumerated) + Scenarios > "Calculate menu walk per dataset" step 1 (dataset open as the per-cell loop entry) | preserved (split for clarity — provisioning declaration in Setup; the open action in scenario body) |
| Step 1 (first) sub-bullet — "smiles_2_columns.csv (TODO: add to linked datasets)" | Setup step 1 (provisioning — dataset D4 enumerated at verified path `System:AppData/Chem/tests/smiles_2_columns.csv`) | preserved (split for clarity — TODO resolved: file verified at the `System:AppData/Chem/tests/smiles_2_columns.csv` path per repo search at `public/packages/Chem/files/tests/smiles_2_columns.csv`; provisioned as Variant D) |
| Step 1 (second, collision) — "Run Chem > Calculate > Descriptors" | Scenarios > "Calculate menu walk per dataset" step 2 (instantiated for cell (D, M1) where M1 = Descriptors top-menu) | preserved (renumbered; "Descriptors" anchored as M1 in the menu-item matrix) |
| Step 2 — "Select arbitrary values in the dialog." | Scenarios > "Calculate menu walk per dataset" step 3 (instantiated per cell — selects arbitrary values in each dialog, with column-selection guidance for Variant D and at-least-one-checkbox guidance for checkbox dialogs M1/M2/M3/M10) | preserved (split for clarity — Variant D column-selection branch and checkbox-dialog guidance added per dialog-shape variance) |
| Step 3 — "Click OK." | Scenarios > "Calculate menu walk per dataset" step 4 (instantiated per cell — OK click) | preserved |
| Step 4 — "Make sure that the column with calculated values is added to the table." | Scenarios > "Calculate menu walk per dataset" step 5 (instantiated per cell — verify expected column appended; row count unchanged; no console errors) | preserved as verification |
| Prose "Do the asme for each section in the Calculate menu ()" | Scenarios > "Calculate menu walk per dataset" — enumerated as the M1–M10 Calculate-menu-item matrix; Steps 2–5 instantiate per cell (D, M) | preserved (split for clarity — typo "asme" → "same" silently fixed; empty parens "()" silently dropped; "each section" enumerated explicitly as 10 menu items per atlas `chem.calculate.*` + `chem.mpo.top-menu` + `chem.notation.inchi` + `chem.notation.inchi-keys` ids) |
| JSON footer `{ "order": 3, "datasets": ["System:DemoFiles/chem/smiles.csv", "System:AppData/Chem/mol1K.sdf", "System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf"] }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` + migrated Setup step 1 dataset enumeration — extended to 4 entries including the resolved `smiles_2_columns.csv` Variant D) |

No original step is silently dropped. Step-numbering collision (original had two
distinct `1.` markers and then `2.` / `3.` / `4.`) renumbered cleanly into the
migrated body's `Setup` + `Scenarios > "Calculate menu walk per dataset"` linear
sequence.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/chem.yaml` `output_plan.calculate.md.target_layer = playwright`.
  The scenario requires DOM-level UI driving: top-menu walk (**Chem > Calculate >
  ${item}**) for each of the 10 Calculate menu items, in-dialog input selection
  (column pickers, per-property checkbox toggles), OK click on each dialog, and
  grid-column-appended verification on the active table view — all UI behaviors
  that Playwright drives natively and the apitest layer cannot exercise as a
  matrix walk. Chain YAML reason field: "Multi-format / multi-Calculate-menu-item
  walk. Sibling calculate-spec.ts already exists at playwright per
  existing-test-index. Migrator must enumerate Calculate menu items (currently
  prose 'Do the same for each section' — surfaced under unresolved_ambiguities);
  resolution lands in the migrated .md and downstream spec. Single scenario;
  simple strategy."
- **Why this `coverage_type`:** chose `regression` per chain rev 2 directive
  (`pyramid_layer: integration` — multi-format × multi-menu-item walk; not
  `smoke` since the section's smoke is `Advanced/scaffold-tree-functions.md`
  per chain `ui_coverage_plan.smoke_scenario`; not `edge` / `perf` since the
  scenario walks the happy path only — no specific failure-mode invariant
  being asserted, no timing or volume threshold). The Docker-unavailable /
  timeout failure-mode invariant for `chem.calculate.descriptors`
  (`GROK-17621`) is parallel-coverage owned by the dedicated bug-focused
  candidate spec `chem-grok-17621-spec.ts` per chain
  `bug_focused_candidates[]`; this scenario assumes the `chem-chem` container
  is available. `regression` is the natural fit for a multi-format walk that
  guards against regression-of-the-set across the full Calculate menu.
- **Why this `strategy`:** `simple` per chain YAML
  `output_plan.calculate.md.strategy = simple`. Single scenario, no
  cross-file fixture; chain analyzer classified `medium` (multi-format /
  multi-subsystem) but kept strategy `simple` — single migrated `.md`,
  single downstream spec. Pattern 1 (decomposition) is NOT triggered (chain
  classifier returned `medium`, not `complex-standalone`; no child scenarios
  with `produced_from: decomposed`). Pattern 2 (bug-focused slice) does NOT
  apply to this scenario itself — `calculate.md` IS a span for two
  bug-focused candidates (`chem-grok-17621-spec.ts` and
  `chem-github-2942-spec.ts`) but those are downstream Automator concerns;
  this scenario stays a parent. Pattern 3 (UI-smoke consolidation) is NOT
  triggered (no `ui_consolidation_proposals[]` from Gate F). Standard
  per-scenario migration per "Migrated scenario shape".
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/calculate-spec.ts`
    — existing playwright-layer test ("Chem: Calculate Descriptors",
    category Chem, `features_covered: [chem.calculate]`, helpers `spec-login`,
    patterns `uses-fixture-runners` / `uses-grok.dapi` / `uses-grok.shell` /
    `uses-page.evaluate` / `uses-page.locator` / `uses-playwright-test`).
    Per `existing-test-index.yaml` line 32409-32426. Existing spec covers
    Descriptors only; Automator will extend to the full enumerated menu walk.
    Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts`
    — sibling Chem playwright test ("Chem: Info Panels on smiles.csv").
    Same shape pattern; confirms the section-wide playwright convention.
    Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/scaffold-tree-functions.md`
    + `scaffold-tree-functions-migration-report.md` (2026-05-11 predecessor
    in same migration cycle) — shape anchor for the migrated `.md`
    (frontmatter, Setup / Scenarios / Notes order) and the SR-01 A-STRUCT-02
    carryforward template adopted here. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/scaffold-tree.md`
    + migration report — same-cycle section-mate; confirms
    `coverage_type: regression` precedent for an integration-pyramid scenario
    in this section. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/r-group-analysis.md`
    + migration report — same-cycle section-mate; confirms section-wide
    multi-block scenario house-style + frontmatter shape. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/activity-cliffs.md`
    — adjacent Chem section scenario that references the same
    `smiles_2_columns.csv` TODO note (line 3 of body). Sibling scenario
    keeps the TODO open and uses the
    `System:AppData/UsageAnalysis/test_datasets/smiles_2_columns.csv` path
    in JSON footer (line 21); this migration resolves the TODO via the
    verified file path at `System:AppData/Chem/tests/smiles_2_columns.csv`
    (per repo search confirming the file at
    `public/packages/Chem/files/tests/smiles_2_columns.csv`).
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok`
      (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`)
      — section-standard fixture login; Automator will use in spec
      `beforeEach`.
    - `softStep`
      (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`)
      — per-step soft-fail wrapper; Automator wraps each cell verification
      in `softStep` for granular failure reporting (especially important
      for a 40-cell matrix where individual cells may surface user-visible
      error balloons rather than column-append).
    - `closeAllViews`
      (`public/packages/UITests/playwright/scripts/helpers.ts:129`)
      — cleanup between cells of the matrix walk (per "Close the active
      view before moving to the next (D, M) cell" instruction in the
      migrated body).
  - **Candidate helpers (NOT yet in registry — flagged for addition via
    helpers-registry curator; per migration-prompt Helpers discipline §):**
    - `helpers.playwright.chem.openChemTopMenuItem(page, menuPath)`
      — walks **Chem > Calculate > ${menuItem}** (or any top-menu path).
      Surfaced by Scenarios step 2 per cell. Convenience abstraction; spec
      Author can inline the top-menu walk until the helper lands.
    - `helpers.playwright.chem.runCalculateDialogWithDefaults(page,
      dialogTitle)` — accepts an open Calculate dialog, selects sensible
      defaults (auto-detected molecule column; at-least-one checkbox where
      applicable), clicks OK, and waits for the column-append. Surfaced by
      Scenarios steps 3-5 per cell.
    - `helpers.playwright.chem.assertColumnAppended(page, columnName |
      columnNamePredicate, prevColumnCount)` — asserts the grid column
      count grew (or a specifically-named column is present) without
      changing row count. Surfaced by Scenarios step 5 per cell.
- **Bug library consulted:** yes — `bug-library/chem.yaml` revision 1
  (12 curated_bugs, generated 2026-05-05). One bug intersects this
  scenario's flows but does NOT trigger inclusion in `related_bugs` per
  the chain YAML's frontmatter delta — it is owned by a dedicated
  bug-focused spec in chain `bug_focused_candidates[]`:
  - **GROK-17621** (Chem: Descriptors: Timeout of waiting for container
    status change). `affects` intersects `chem.calculate.descriptors`,
    `chem.calculate.descriptors.top-menu`, `chem.calculate.descriptors.transform`
    — three sub_features in this scenario's `sub_features_covered`.
    However the bug's repro requires the `chem-chem` Docker container
    to be unavailable; `calculate.md` walks the happy path only —
    container assumed available per the migrated Setup step 2. Chain
    emits `chem-grok-17621-spec.ts` as the bug-focused candidate
    (spans `calculate.md:Step 1`); awareness only here. Cited in
    migrated Notes section.
  - **github-2942** (CSV export with "Molecules As Smiles" + active
    filtering fails with rowCount mismatch). `affects` intersects
    `calculate.md` via the `chem.notation.convert-column` /
    `convertMoleculeNotation` dependency invoked by Calculate menu
    walks (chain `bug_focused_candidates[]` spans `calculate.md:Step
    1` + `filter-panel.md:Step 1`). The dedicated
    `chem-github-2942-spec.ts` candidate spec owns the filter+export
    combination; `calculate.md` does not filter or export, so this
    awareness only. Cited in migrated Notes section.
  - Remaining curated bugs (GROK-16329 R-Groups, GROK-18517 MMP,
    GROK-12758 Scaffold Tree integration, GROK-17595 Save+reopen
    deserialization, GROK-18407 Chemical Space silent dialog close,
    GROK-17964 Convert Notations duplicate, GROK-19624 MPO Profile
    refresh, github-3004 Multi-table Scaffold Tree wrong-binding,
    GROK-16870 RDKit cell renderer in non-Chem viewer, GROK-14028
    Filter Panel Reset cleanup gap) — `affects` disjoint from this
    scenario's `sub_features_covered`. No overlap.
- **Decision log queried:** yes — `decision-log.yaml` grepped for
  `feature: chem` — **zero matches** (consistent with the
  scaffold-tree-functions / scaffold-tree / r-group-analysis migration
  reports in the same cycle). This is the FIRST CYCLE for the chem
  section. No prior `migration_decisions`, `layer_decisions`,
  `manual_only`, or `failed_attempts` entries exist for chem. No
  "approaches off the table" constraints apply.
- **Cross-cutting bug citations (chain YAML `bug_focused_candidates[]`).**
  Per `migration-prompt.md` "Cross-cutting bug citations from chain YAML"
  §: this scenario is a span (or affected scope) in two candidates —
  GROK-17621 (`chem-grok-17621-spec.ts`, Descriptors Docker timeout) and
  github-2942 (`chem-github-2942-spec.ts`, CSV export with filtering).
  Both cited in migrated Notes section ("Cross-cutting bug awareness —
  GROK-17621 / github-2942"). Per migration-prompt this is RECOMMENDED,
  not mandatory; F-BUG-COVERAGE-01 at section-complete is authoritative.
- **UI delegation status.** Per chain YAML
  `ui_coverage_plan.smoke_scenario: Advanced/scaffold-tree-functions.md`
  (this scenario is NOT the section smoke) and the chain entry
  `ui_coverage_plan.delegations[].scenario: calculate.md` /
  `delegated_to: null`, this scenario owns its UI coverage directly
  (no upstream delegation). Chain
  `ui_coverage_responsibility: [chem-calculate-descriptors,
  chem-calculate-properties, chem-calculate-toxicity-risks,
  chem-calculate-mpo-score, chem-calculate-bitbirch-clustering,
  chem-calculate-cluster-mcs, chem-calculate-map-identifiers,
  chem-calculate-to-inchi, chem-calculate-to-inchi-keys]`. No
  SCOPE_REDUCTION proposal substitutes JS API for any of these UI
  flows; therefore the UI delegation discipline of migration-prompt §
  "UI delegation in SCOPE_REDUCTION proposals" does not apply (no SR
  proposal to delegate). Reaffirmed as a hard constraint in the
  migrated Notes ("No JS API substitution").
- **Scenario constraint extraction (per migration-prompt § Scenario
  constraint extraction).**
  - **(a) FORBIDDEN substitutions:** None explicit in the original body
    (no "must not", "forbidden", "not allowed" phrasing). HOWEVER
    `pyramid_layer: integration` per chain implicitly forbids JS API
    substitution for the 9 `ui_coverage_responsibility` flows —
    surfaced explicitly in the migrated Notes (`No JS API
    substitution`) as a structured behavioral guard. Additionally the
    Calculate menu walk via JS API
    (`grok.functions.eval('Chem:descriptorsDocker')` etc.) would defeat
    the dialog-driven workflow that the scenario asserts (input
    selection, OK click).
  - **(b) REQUIRED actions:** All 9 entries of
    `ui_coverage_responsibility` are exercised via UI driving in the
    migrated body (top-menu walk per cell). Biochemical Properties
    (M10) is added beyond the chain `ui_coverage_responsibility` list
    per the enumeration directive ("each section in the Calculate
    menu") — atlas `chem.calculate.biochemical-properties` is a
    registered Calculate top-menu entry per source line 2860, so the
    enumeration includes it for fidelity to the original directive.
    No flow is deferred.
  - **(c) Missing-selector escalation:** Migrator does NOT have
    explicit selector definitions for the per-Calculate-menu-item
    dialogs in `grok-browser/references/projects.md`. Spec-time
    selector discovery is the Automator's responsibility (per
    existing sibling `calculate-spec.ts` patterns + `page.evaluate`
    introspection); a reference-file approval-required proposal to
    add Calculate-dialog-specific selectors to
    `grok-browser/references/widgets/calculate-dialogs.md` is a
    Phase 2 deliverable — out of this migration's scope.
  - **(d) Reference templates:** Code-style anchors for the Automator
    (cited verbatim, not modified here):
    `public/packages/UsageAnalysis/files/TestTrack/Chem/calculate-spec.ts`
    (current existing test — section-mate Descriptors-only spec;
    Automator extends per the enumerated menu walk) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts`
    (section-mate playwright walk; `spec-login` + `page.evaluate` +
    `page.locator` patterns).
  - **(e) Scenario authority clause:** No prompt-vs-scenario conflict
    surfaced. The original scenario is a 4-step linear walk anchored
    on Descriptors plus a prose directive to walk "each section in
    the Calculate menu"; the migrated body preserves all 4 anchor
    steps and resolves the prose directive into an explicit 10-item
    menu enumeration. No `prompt_scenario_conflict:` field to flag for
    orchestrator.
- **Source-text fixes (silently applied per chain rev 2 directive,
  Olena 2026-05-11 chain footer note (i)).**
  - **Typo `asme` → `same`.** Original line: "Do the asme for each
    section in the Calculate menu ()". Migrated body uses "same"
    (resolution embedded in the enumerated menu-walk matrix; the
    phrase itself is not preserved verbatim because the directive is
    realized structurally as the M1–M10 table).
  - **Empty parens `()`.** Original line ends "Calculate menu ()" —
    appears to be an editor artifact. Empty parens dropped.
  - **`smiles_2_columns.csv` TODO.** Original sub-bullet: "smiles_2_columns.csv
    (TODO: add to linked datasets)". TODO resolved: the file exists at
    `System:AppData/Chem/tests/smiles_2_columns.csv` (verified by repo
    search at `public/packages/Chem/files/tests/smiles_2_columns.csv`).
    Provisioned in Setup step 1 as dataset D4 (Variant D — 2-column
    SMILES). The sibling activity-cliffs.md scenario flags the same
    TODO and uses a different (UsageAnalysis-side) path; this
    migration prefers the Chem-owned path since the file lives under
    `Chem/files/tests/`.
  - **Step-numbering collision.** Original has two `1.` markers (the
    dataset-open block and the Descriptors-walk anchor) and then
    `2.` / `3.` / `4.`. Renumbered cleanly into the migrated body's
    `Setup` + `Scenarios > "Calculate menu walk per dataset"` structure.
  - **Calculate menu enumeration.** Original prose "Do the same for
    each section in the Calculate menu" enumerates implicitly. Chain
    rev 2 (Olena 2026-05-11 directive) requires Migrator to make this
    explicit. Resolved per atlas `chem.calculate.*` ids + `chem.mpo.top-menu`
    + `chem.notation.inchi` + `chem.notation.inchi-keys` (per source
    `package.ts` `Calculate` top-menu entries at lines 570 / 612 / 695 /
    725 / 1243 / 1265 / 2146 / 2199 / 2732 / 2860): 10 menu items —
    Descriptors / Chemical Properties / Toxicity Risks / MPO Score /
    BitBIRCH Clustering / Cluster MCS / Map Identifiers / To InchI /
    To InchI Keys / Biochemical Properties.
- **Parameterized data_driven strategy.** The migrated body presents
  the matrix as a 4×10 cell walk (4 datasets × 10 menu items) per
  migration-prompt § "Matrix scenarios from the chain". The downstream
  spec author (Automator) is expected to realize this as a Playwright
  `test.describe` block iterating datasets × menu items (or a
  `test.each`-style fixture) with per-cell soft-step verification.
- **No invented sub_features / helpers.** Every id in
  `sub_features_covered` is present in `feature-atlas/chem.yaml` rev 2
  (verified per atlas lines 542-716 for `chem.calculate.*` ids;
  719-732 for `chem.mpo` / `chem.mpo.top-menu`; 875-887 for
  `chem.notation.inchi` / `chem.notation.inchi-keys`). Every reused
  helper is in `helpers-registry.yaml`; candidate helpers are flagged
  as candidates, NOT referenced by name in the migrated body.

## Opt-outs (SCOPE_REDUCTION proposals)

### SR-01: A-STRUCT-02 carryforward (chain-level edge/perf coverage)

**Cited technical dependency:** A-STRUCT-02 ("at least one scenario is
edge case or negative path — i.e. its frontmatter `coverage_type` is
`edge` or `perf`") is a **section-level structural invariant**. This
scenario (calculate.md) is a multi-format × multi-menu-item regression
walk — `coverage_type: regression` is the natural fit per chain rev 2
(`pyramid_layer: integration`, classification `medium`, multi-format /
multi-subsystem walk). Forcing `coverage_type: edge` would mis-classify
(the scenario asserts happy-path column-append across the full menu,
not a specific failure-mode invariant).

**Resolution path:** A-STRUCT-02 satisfaction is chain-wide via 10
`bug_focused_candidates[]` in `scenario-chains/chem.yaml` rev 2 (each
chem-bug-focused spec naturally carries `coverage_type: edge` once
authored — including `chem-grok-17621-spec.ts` for the Descriptors
Docker timeout invariant directly relevant to this scenario's
Descriptors anchor, and `chem-github-2942-spec.ts` for the filter+export
combination intersecting this scenario via convert-column). Additionally
`r-group-analysis.md` (same-cycle, same section, `coverage_type: edge`)
already organically satisfies the section-level invariant. Decision-log
line 8528-8530 records the radar.md 2026-05-07 carryforward precedent.

**Scope of waiver:** A-STRUCT-02 verdict on `calculate.md` is deferred
to chain-level evaluation. Critic A should return SR (carryforward),
not FAIL.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through A-STRUCT-06,
A-COVERAGE-*, A-MERIT-*) are expected to PASS on `calculate.md` without
any other SR proposals. No step is opted out for effort. The single
entry above cites a real section-level structural property whose
satisfaction path is owned by the chain, not by this per-scenario
migration (D-MERIT-01 compliant)._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist yet.
The migrated body's 4-dataset × 10-menu-item walk is realizable against
the current playwright + helpers-registry + atlas state. The candidate
helpers surfaced in Decisions are convenience abstractions that the
Automator may inline pending registration — they do NOT block spec
realization. The selector-reference write proposal (Decisions § Scenario
constraint extraction (c)) is a downstream Phase 2 enhancement, not a
blocker on this scenario's automation (existing `calculate-spec.ts`
already drives the Descriptors dialog via `page.evaluate` +
`page.locator` per the sibling-test pattern; Automator can extend that
pattern across the full enumerated menu walk without a reference
update). The Docker container availability assumption (Setup step 2)
is a runtime precondition, not a deferred step — failure-mode coverage
is parallel via `chem-grok-17621-spec.ts`._

## Edge cases

The original lists no explicit "edge case" keyword. Implicit edge cases
derivable from the scenario content + chain-level matrix:

- **Cross-format calculator invocation.** Each Calculate menu item is
  invoked across 4 distinct notation formats (smiles, molV2000,
  molV3000, 2-column SMILES). Implicit invariant: every calculator
  must support every notation format, OR surface a clear
  user-visible error (the per-cell assertion in the migrated body —
  "either column is appended OR a user-visible error balloon is
  shown — never silent"). PRESERVED as scenario step (Scenarios step
  5 cross-cell invariant clause).
- **Variant D 2-column molecule input.** The `smiles_2_columns.csv`
  dataset has two molecule columns; calculator dialogs must present a
  column-selection input and handle the user's choice correctly.
  PRESERVED as scenario step (Scenarios step 3 — "Variant D — two
  molecule columns — pick the first molecule column").
- **Empty / no-selection dialog OK.** Calculator dialogs with
  per-property / per-descriptor checkboxes (Descriptors / Chemical
  Properties / Toxicity Risks / Biochemical Properties) may allow OK
  to be clicked with no checkboxes selected. PRESERVED as scenario
  step guidance (Scenarios step 3 — "select at least one item to
  ensure a non-empty calculation request; full defaults are acceptable
  when the dialog pre-selects items"). The strict "no-selection OK"
  edge case (does the dialog reject? produce no columns? show an
  error?) is NOT exercised — flagged as Atlas-curator candidate edge
  case ("Calculate dialog with no checkbox selected: expected behavior
  is dialog rejects with user-visible message; not in scenario scope").
- **Docker container unavailability for Descriptors (M1).** The happy
  path assumes the `chem-chem` container is up. The container-down
  invariant is OUT OF SCOPE for this scenario — owned by the dedicated
  `chem-grok-17621-spec.ts` per chain `bug_focused_candidates[]`.
  Cited in migrated Notes (`Cross-cutting bug awareness — GROK-17621`).
- **MPO profile absence (M4).** The MPO Score dialog requires at least
  one MPO profile. On environments without bundled profiles the dialog
  surfaces a user-visible message. PRESERVED in Notes (`MPO Score
  profile dependency`).
- **Chembl connectivity for Map Identifiers (M7).** Map Identifiers
  depends on the Chembl package. On disconnected environments the
  dialog surfaces a user-visible error. PRESERVED in Notes (`Map
  Identifiers dependency`).
- **Biochemical Properties calculator-registry absence (M10).** The
  Biochemical Properties dialog discovers functions tagged
  `function_family: biochem-calculator`. On environments without such
  functions the dialog may surface empty / informational state.
  PRESERVED in Notes (`Biochemical Properties dependency`).
- **CSV export with filter applied (cross-cutting).** Not exercised by
  this scenario (no filter, no export). OUT OF SCOPE — owned by
  `chem-github-2942-spec.ts`. Cited in migrated Notes
  (`Cross-cutting bug awareness — github-2942`).
- **Console errors throughout.** Implicit across all 40 cells of the
  matrix: each cell's dialog open + OK click + column append must
  complete without console errors. PRESERVED as a per-cell verification
  (Scenarios step 5 — "No console errors fire during dialog open, OK
  click, or column-append"). Sister with `info-panels.md` /
  `scaffold-tree-functions.md` patterns.

No edge case is moved to atlas, manual_only, deferred, or a separate
scenario silently. The strict "no-selection OK" calculator dialog edge
case is flagged for atlas curator (no atlas write here).

## Unresolved ambiguities

- **`smiles_2_columns.csv` path canonical location.** The migrated body
  uses `System:AppData/Chem/tests/smiles_2_columns.csv` (verified path
  per repo search at `public/packages/Chem/files/tests/smiles_2_columns.csv`).
  Sibling `activity-cliffs.md` uses
  `System:AppData/UsageAnalysis/test_datasets/smiles_2_columns.csv`
  (line 21 of JSON footer). Both paths may exist on a deployed
  Datagrok instance (or one may be a duplicate / older path). The
  Automator at spec time should verify which path the running test
  server exposes; if both exist, prefer the Chem-owned path since the
  file is bundled with the Chem package. Flag for QA pair review.
- **"Select arbitrary values in the dialog" — non-deterministic UI
  driving across 10 distinct dialogs.** Each Calculate menu item has
  its own dialog shape: Descriptors has a tree of descriptor classes
  with per-class checkboxes (`chem.calculate.descriptors.tree`);
  Chemical Properties has 9 per-property checkboxes; Toxicity Risks
  has 4 per-risk checkboxes; MPO Score has profile + aggregation
  inputs; BitBIRCH has threshold + fingerprint-type inputs; Cluster
  MCS has cluster-column input; Map Identifiers has source / target
  identifier selectors; To InchI / To InchI Keys have minimal inputs
  (column only); Biochemical Properties discovers calculators
  dynamically. The migrated body's step 3 ("accept arbitrary values
  ... select the auto-detected molecule column where prompted ... at
  least one checkbox where applicable") is non-deterministic — the
  Automator at spec time picks per-dialog deterministic input shapes.
  Flag for QA pair review.
- **Expected column name(s) per Calculate menu item.** The migrated
  body's matrix table lists likely-appended columns (e.g. "MolWt",
  "LogP" for M1; "MW, HBA, HBD ..." for M2) but the exact column
  names depend on user input selections in each dialog. The assertion
  shape in step 5 is "the expected column ... is appended", which
  bounds to "the column count grew and a new column with a
  calculator-recognizable name was added" — not "an exactly-named
  column equals X". Automator at spec time may use a column-count-
  delta assertion (`afterCols.length > beforeCols.length`) or a
  predicate assertion on column-name patterns. Flag for QA pair
  review.
- **Cell ordering inside the matrix.** The migrated body presents M1–M10
  in source-`package.ts`-declaration order (Descriptors first, then
  Properties, etc.) but does NOT enforce a strict sequential ordering
  across cells. The Automator at spec time can parallelize cells per
  dataset or run them sequentially; the scenario's invariants are
  per-cell (no inter-cell dependencies). Flag for spec design at
  Automator stage.
- **Variant D molecule-column selection.** The migrated body says
  "pick the first molecule column" for Variant D. Both columns in
  `smiles_2_columns.csv` are molecule columns; the choice is
  arbitrary but should be consistent across cells. Automator at spec
  time picks a deterministic column (e.g. by column name or
  positional index). Flag for QA pair review.
