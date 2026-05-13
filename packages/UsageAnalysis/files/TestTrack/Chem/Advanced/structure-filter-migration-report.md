# Migration Report — Advanced/structure-filter.md

First-cycle migration of `Advanced/structure-filter.md` (TestTrack/Chem
section, scenario #7 in `order_from_files` per chain rev 2) per chain
YAML `scenario-chains/chem.yaml`
(`output_plan.Advanced/structure-filter.md`: `target_layer: playwright`,
strategy `simple`, classification `medium`, `pyramid_layer: integration`,
`depends_on: []`). 44-line scenario with 5 sub-blocks separated by `***`
fences — each block a distinct substructure-filter invariant
(disable/re-enable through panel close+reopen; Use-as-filter after panel
close; hamburger-menu Filter + Draw + Add filter; Remove + Use as filter
ordering; clone-view cross-view sync). Realizes atlas critical path
`chem.cp.substructure-search-with-filter` (p0).

## Step mapping

The original is 44 lines, 5 `***`-fenced sub-blocks, each containing
its own renumbered step list. Source-text defects preserved verbatim
in original: Block 1 numbering glitch `1. 2. 3. 3. 3.`; Block 4 typo
"corectly" → "correctly"; Block 5 numbering glitch `4. 4.`. These
are silently fixed per chain rev 2 directive (Olena 2026-05-11 —
"Source-text defects in scenario bodies — Migrator fixes during
migration"). JSON footer carries `order: 4` + single dataset entry —
metadata, not a step; captured in chain `order_from_files` line 100 and
migrated `## Setup` step 1.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Preamble "Open linked datasets" (Block 1 line 1) | Setup step 1 (dataset enumerated as `System:DemoFiles/chem/SPGI.csv` per JSON footer) | preserved (split for clarity — generic "linked datasets" specialized to the JSON footer's single entry; Block 1 then opens the Filter Panel via `chem-structure-filter-draw` UI flow) |
| Block 1 step "1. Add a structure filter (draw or current value > use as filter)" | Scenarios > "Disable filter + close+reopen panel + re-enable" step 1 | preserved |
| Block 1 step "2. Disable the structure filter" | Scenarios > "Disable filter + close+reopen panel + re-enable" step 2 | preserved |
| Block 1 step "3. Close the filter panel" | Scenarios > "Disable filter + close+reopen panel + re-enable" step 3 | preserved (numbering glitch `1. 2. 3. 3. 3.` silently fixed — original repeated `3.` thrice; renumbered cleanly 1-5) |
| Block 1 step "3. Open the filter panel" (with stray double space in original "Open  the") | Scenarios > "Disable filter + close+reopen panel + re-enable" step 4 | preserved as verification (added implicit verification: filter state preserved across panel close+reopen; double-space typo silently fixed) |
| Block 1 step "3. Turn on the structure filter" | Scenarios > "Disable filter + close+reopen panel + re-enable" step 5 | preserved as verification (added implicit verification: table re-filters identically to pre-disable state, no console errors) |
| Block 2 step "1. Set a filter" | Scenarios > "Set filter then Use as filter after panel close" step 1 | preserved |
| Block 2 step "2. Close the filter panel" | Scenarios > "Set filter then Use as filter after panel close" step 2 | preserved |
| Block 2 step "3. Use Current value > Use as filter. Check that filter UI is correct" | Scenarios > "Set filter then Use as filter after panel close" step 3 (right-click → cell context-menu → verify sketcher UI shows molecule sketched correctly + filter re-applies + no console errors) | preserved (split for clarity — action + verification rolled into one numbered step per the original's compact form) |
| Block 3 step "1. Set a filter" | Scenarios > "Hamburger menu Filter + Draw + Add filter" step 1 | preserved |
| Block 3 step "2. Close the filter panel" | Scenarios > "Hamburger menu Filter + Draw + Add filter" step 2 | preserved |
| Block 3 step "3. Go to the mol column's hamburger menu > filter" | Scenarios > "Hamburger menu Filter + Draw + Add filter" step 3 | preserved (atlas alignment — "filter" lower-case in original; canonical menu phrasing **Filter** in migrated body) |
| Block 3 step "4. Draw another structure" | Scenarios > "Hamburger menu Filter + Draw + Add filter" step 4 | preserved |
| Block 3 step "5. Click Add filter" | Scenarios > "Hamburger menu Filter + Draw + Add filter" step 5 (action + verification — combined filter applied, second Structure filter appears on panel, table re-filters, no errors) | preserved as verification |
| Block 4 step "1. Remove Structure filter" | Scenarios > "Remove + Use as filter ordering" step 1 | preserved |
| Block 4 step "2. Current value > Use as filter - check:" (with sub-bullets) | Scenarios > "Remove + Use as filter ordering" step 2 (right-click → cell context-menu → two-bullet verification — molecule sketched correctly AND Structure filter first on Filter Panel) | preserved as verification (sub-bullets preserved as bullet list under step 2; "corectly" → "correctly" silent typo fix per chain directive; numbering glitch `4. 4.` in original — there is no Block 4 "step 4" duplicate — analyzed as a labeling mix-up between blocks 4 and 5, both numbered "4." at the block-header level; renumbered cleanly within Block 4 as steps 1-2 + an added implicit "no console errors" verification step 3) |
| Block 4 sub-bullet "the molecule is sketched corectly" | Scenarios > "Remove + Use as filter ordering" step 2 bullet 1 (rendered as nested bullet under the right-click verification step) | preserved (typo silently fixed) |
| Block 4 sub-bullet "the Structure filter is the first on the Filter Panel" | Scenarios > "Remove + Use as filter ordering" step 2 bullet 2 | preserved |
| (added implicit verification — "no console errors fire during the remove → re-add transition") | Scenarios > "Remove + Use as filter ordering" step 3 | added as implicit verification per scenario-wide no-console-errors invariant (Block 4 numbered cleanly 1-3 after silent fix) |
| Block 5 step "1. Open linked datasets" | Scenarios > "Clone view + cross-view filter sync" step 1 (open SPGI.csv + open Filter Panel + confirm structure column added) | preserved (split for clarity — the Setup-level "Open SPGI.csv" provisioning is shared, but Block 5 explicitly re-opens the panel to confirm the structure column is added before the clone-view action; treated as a Block-5-local re-assertion) |
| Block 5 step "2. Open filters and make sure some structure column is added to the filter panel" | Scenarios > "Clone view + cross-view filter sync" step 1 (rolled into the open + verify per the original's wording) | preserved as verification |
| Block 5 step "3. Add another view with the same table (e.g. clone current view)" | Scenarios > "Clone view + cross-view filter sync" step 2 | preserved |
| Block 5 step "4. Open filters and make sure same structure column is added to the filter panel in the 2nd view" | Scenarios > "Clone view + cross-view filter sync" step 3 | preserved as verification |
| Block 5 step "4. Right-click some structure > Current value > Use as filter" | Scenarios > "Clone view + cross-view filter sync" step 4 (action + verification — filter applies and synchronizes across both views per the trailing Expected result line) | preserved (numbering glitch `4. 4.` silently fixed — original numbered two consecutive steps as "4."; renumbered cleanly 1-4) |
| Block 5 trailing line "Expected result: filter is applied, filters are synchronized in two open views" | Scenarios > "Clone view + cross-view filter sync" step 4 (folded into the action step's verification clause) | preserved as verification (D-STEP-02 — Expected result becomes verification step) |
| JSON footer `{ "order": 4, "datasets": ["System:DemoFiles/chem/SPGI.csv"] }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` line 100 + migrated Setup step 1 dataset enumeration) |

No original step or implicit verification is silently dropped.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/chem.yaml`
  `output_plan.Advanced/structure-filter.md.target_layer = playwright`.
  The scenario requires DOM-level UI driving across five distinct
  surfaces: Filter Panel toggle / sketcher draw / cell context-menu
  (Current value > Use as filter) / column header hamburger menu /
  view-clone affordance / cross-view filter inspection — all UI
  behaviors that Playwright drives natively and the apitest layer
  cannot exercise. Chain YAML reason field: "Substructure Filter cycle
  states + clone-view sync. Sibling Advanced/structure-filter-spec.ts
  already exists at playwright per existing-test-index. Clone-view
  + cross-view filter sync (Block 5) requires real view-state DOM
  that playwright provides. Single scenario across 5 blocks; simple
  strategy."
- **Why this `coverage_type`:** chose `regression` per task directive
  ("multi-invariant substructure-filter cycle walk; not bug-focused").
  The scenario walks five distinct substructure-filter invariants on
  a single viewer/panel surface with `pyramid_layer: integration`. It
  is not a single-bug repro (`related_bugs: []` — GROK-14028 is
  routed to the dedicated bug-focused spec; see Decisions bullet
  "Bug library consulted" and "Cross-cutting bug citations" below),
  not a one-flow smoke (five sub-scenarios), and not a perf test (no
  thresholds). `regression` matches a multi-invariant walk protecting
  against feature surface erosion as Filter Panel state machine,
  hamburger-menu wiring, cell-action ordering, or cross-view sync
  semantics evolve. Precedent: `similarity-search.md` (regression —
  multi-property-knob walk) + `info-panels.md` (regression —
  multi-format Context Panel walk) + `calculate.md` (regression —
  multi-format × multi-menu-item walk), all same-cycle same-section.
- **Why this `strategy`:** `simple` per chain YAML
  `output_plan.Advanced/structure-filter.md.strategy = simple`. Single
  scenario across 5 sub-blocks → single migrated `.md` with 5
  `### <name>` sub-scenarios, single downstream spec. Pattern 1
  (decomposition) is NOT triggered (chain classifier returned
  `medium`, not `complex-standalone`; no child scenarios with
  `produced_from: decomposed`; preserved as 5 sub-scenarios in one
  file per task directive). Pattern 2 (bug-focused slice) applies
  for GROK-14028 but the slice is owned by the dedicated
  `chem-grok-14028-spec.ts` candidate, NOT authored here (Migrator's
  pattern-2 role is surfacing only — see "Cross-cutting bug citations"
  below). Pattern 3 (UI-smoke consolidation) is NOT triggered (no
  `ui_consolidation_proposals[]` from Gate F; the section's UI smoke
  is `Advanced/scaffold-tree-functions.md`). Standard per-scenario
  migration per "Migrated scenario shape".
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/structure-filter-spec.ts`
    — existing playwright-layer test ("Chem: Structure Filter
    substructure search + panel toggle + view sync", category Chem,
    helpers_called `spec-login`, patterns `uses-fixture-runners` /
    `uses-grok.dapi` / `uses-grok.shell` / `uses-page.evaluate` /
    `uses-page.locator` / `uses-playwright-test`). Per
    `existing-test-index.yaml` line 32392-32407. Confirms the
    section's playwright + `spec-login` + `page.evaluate` / `page.locator`
    convention; the spec is the Automator's code-style anchor for the
    Substructure Filter UI driving. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/similarity-search.md`
    + `similarity-search-migration-report.md` — adjacent same-section
    same-cycle migration with `coverage_type: regression`,
    `pyramid_layer: integration`, multi-axis walk on a single viewer.
    Shape anchor for the migrated `.md` (frontmatter, Setup /
    Scenarios / Notes order), the SR-01 A-STRUCT-02 carryforward
    template, the Notes section "Cross-cutting bug awareness"
    phrasing convention, and the source-text-defect silent-fix
    convention. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/scaffold-tree-functions.md`
    + `scaffold-tree-functions-migration-report.md` — adjacent same-
    section migrated scenario. Anchor for the SR-01 A-STRUCT-02
    carryforward template, the cross-cutting bug-awareness Notes
    convention, and the `pyramid_layer: ui-smoke`-vs-`integration`
    distinction surfacing in the report's Decisions section.
    Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/filter-panel.md`
    — sibling Chem scenario (root level; not in Advanced/). Per chain
    `bug_focused_candidates[]` GROK-14028 spans both `filter-panel.md:Step 2`
    and `Advanced/structure-filter.md:Step 1` — i.e. they share the
    Filter Panel Reset cleanup invariant. The two scenarios are
    distinct: `filter-panel.md` walks the Filter Panel widget set
    breadth (`chem-substructure-filter-draw` +
    `chem-substructure-filter-settings-{contains,included-in,exact,similar}`
    + drag-drop-to-panel + cross-filter-other-types + sketcher-modify-
    from-hamburger), while `Advanced/structure-filter.md` walks the
    Substructure Filter cycle states (disable/re-enable, panel close+
    reopen, hamburger Filter + Draw + Add, Remove + Use as filter
    ordering, clone-view sync). Chain rev 2 notes the two scenarios
    "complement rather than duplicate". Read-only.
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok`
      (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`)
      — section-standard fixture login; Automator will use in spec
      `beforeEach`.
    - `softStep`
      (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`)
      — per-step soft-fail wrapper; Automator wraps each sub-scenario
      step in `softStep` for granular failure reporting (relevant
      for the 5-sub-scenario walk where soft-fails per block aid
      diagnosis).
    - `closeAllViews`
      (`public/packages/UITests/playwright/scripts/helpers.ts:129`) —
      cleanup at end of scenario / between sub-scenarios (especially
      Block 5 after the clone-view).
  - **Candidate helpers (NOT yet in registry — flagged for addition
    via helpers-registry curator; per migration-prompt Helpers
    discipline §):**
    - `helpers.playwright.chem.openFilterPanel(page)` /
      `helpers.playwright.chem.closeFilterPanel(page)` — drives the
      `.grok-icon-filter` toggle on the Table View toolbar
      (`table_view.dart:647` per MEMORY.md Datagrok UI Selector
      Notes); container is `.d4-filter-panel`. Surfaced across all 5
      blocks. Spec author can inline the toggle click pattern until
      helper lands.
    - `helpers.playwright.chem.drawStructureFilter(page, smiles)` —
      opens / focuses the Substructure Filter sketcher widget on the
      Filter Panel and draws (or pastes SMILES into) the given
      structure. Surfaced by Blocks 1, 2, 3.
    - `helpers.playwright.chem.useAsFilterFromCell(page, columnName,
      rowIndex)` — drives the cell context-menu **Current value >
      Use as filter** action on the named molecule cell (atlas
      `chem.search.use-as-filter`, `Chem/src/package.ts#L1957`).
      Surfaced by Blocks 1, 2, 4, 5.
    - `helpers.playwright.chem.hamburgerMenuFilter(page, columnName)` —
      drives the column header hamburger menu's **Filter** action
      to open the column-scoped filter sketcher dialog. Surfaced by
      Block 3 step 3.
    - `helpers.playwright.views.cloneCurrentView(page)` — drives
      the *Clone view* affordance (view-tab context menu) to add a
      second view onto the active table. Surfaced by Block 5 step 2.
- **Bug library consulted:** yes — `bug-library/chem.yaml` rev 1
  (12 curated_bugs, generated_date 2026-05-05). **GROK-14028**
  (Filter Panel Reset cleanup gap — Reset doesn't clear sketcher
  input) `affects` intersects this scenario's coverage
  (`chem.search.substructure.filter`, `chem.sketcher`). However,
  per the task directive and `pattern 2` (bug-focused slice — no
  separate `.md` authored), GROK-14028 is **owned by the dedicated
  bug-focused candidate spec** `chem-grok-14028-spec.ts` (chain
  `bug_focused_candidates[]` entry; `spans:
  [filter-panel.md:Step 2, Advanced/structure-filter.md:Step 1,
  sketcher.md:Step 6]`). Per the task directive
  (`related_bugs=[] (GROK-14028 owned by bug-focused candidate spec,
  noted in Decisions)`) and consistent with `pattern 2` Migrator role
  ("under Pattern 2, no separate `<bug-name>.md` is authored.
  Migrator's role is to … confirm `related_bugs` in the parent scenario
  frontmatter references the relevant GROK ticket(s) — if the bug
  coverage spans multiple steps in the parent, they are listed once
  at parent level; spec authoring (Automator stage) creates a
  satellite `*-spec.ts` with self-documenting header per the
  bug-focused slice pattern"), `related_bugs: []` is set
  intentionally because (a) the bug invariant is owned by the
  dedicated bug-focused spec, NOT this parent scenario's spec, and
  (b) the migrated body's Block 1 IS a defense-in-depth sister
  check ("does panel close+reopen preserve filter state?") rather
  than the GROK-14028 invariant itself ("does Reset clear the
  sketcher UI?") — the two flows test adjacent state-cleanup
  surfaces but distinct invariants. Documented in migrated Notes
  ("Cross-cutting bug awareness — GROK-14028") + this report
  Decisions section per migration-prompt § "Cross-cutting bug
  citations from chain YAML".
  - Remaining curated bugs (GROK-12758 Scaffold Tree × Sketcher ×
    Substructure-Search integration, GROK-17595 Project save+reopen
    Scaffold Tree deserialization, GROK-17621 Descriptors Docker
    timeout, GROK-18517 MMP failure, GROK-16329 R-Groups null,
    GROK-18407 Chemical Space silent dialog close, GROK-17964
    Convert Notations duplicate, github-3004 multi-table Scaffold
    Tree wrong-binding, github-2942 CSV export with filtering,
    GROK-19624 MPO Profile refresh, GROK-16870 RDKit cell renderer
    in non-Chem viewer tooltip) — `affects` disjoint from this
    scenario's `sub_features_covered`. GROK-12758 mentions
    `chem.search.substructure` and `chem.search.substructure.api`
    in `affects` (intersect with this scenario's coverage at the
    API level) but the bug invariant is Scaffold-Tree-specific
    (opening a scaffold node in the sketcher corrupts subsequent
    substructure search); this scenario does NOT exercise
    Scaffold Tree at all. No overlap with this scenario's
    flow-level coverage.
- **Decision log queried:** yes — `decision-log.yaml` grepped for
  `feature: chem` AND `Advanced/structure-filter.md` — no
  `failed_attempts` entries matching this scenario. No "approaches
  off the table" constraints apply.
- **Cross-cutting bug citations (chain YAML
  `bug_focused_candidates[]`).** Per `migration-prompt.md`
  "Cross-cutting bug citations from chain YAML" § + Pattern 2: this
  scenario is a `spans` entry in **GROK-14028**
  (`chem-grok-14028-spec.ts`); the cross-cutting candidate is cited
  in the migrated Notes section ("Cross-cutting bug awareness —
  GROK-14028"). Per migration-prompt this is RECOMMENDED, not
  mandatory; F-BUG-COVERAGE-01 at section-complete is authoritative.
- **UI delegation status.** Per chain YAML
  `ui_coverage_plan.delegated_scenarios.Advanced/structure-filter.md`
  (lines 749-757): `delegated_to: null` (UI coverage owned, not
  delegated to `Advanced/scaffold-tree-functions.md` ui-smoke).
  Rationale verbatim from chain: "UI coverage owned (Substructure
  Filter cycle states: disable / re-enable / close-and-reopen panel
  / hamburger-menu Filter / Remove + Use as filter ordering /
  clone-view sync). Filter-subsystem UI not witnessed by Scaffold
  Tree smoke. Scenario owns its own coverage." All 8
  `ui_coverage_responsibility` flows (chain lines 403-411:
  `chem-structure-filter-draw`,
  `chem-substructure-use-as-filter`,
  `chem-substructure-disable-filter`,
  `chem-substructure-remove-filter`,
  `chem-substructure-hamburger-menu-add-filter`,
  `filter-panel-close-open`, `clone-view`,
  `filter-cross-view-sync`) are owned by this scenario. No UI
  coverage gap introduced; no SCOPE_REDUCTION proposal moves UI
  driving to a sibling — therefore the UI delegation discipline of
  migration-prompt § "UI delegation in SCOPE_REDUCTION proposals"
  does not apply here (no SR proposal to delegate).
- **Source-text defect silent fixes** (per chain rev 2 directive,
  Olena 2026-05-11 — "Source-text defects in scenario bodies —
  Migrator fixes during migration"):
  - **Block 1 numbering glitch** `1. 2. 3. 3. 3.` — original
    numbered five consecutive steps as 1/2/3/3/3; renumbered
    cleanly 1-5 in the migrated `### Disable filter + close+reopen
    panel + re-enable` sub-scenario.
  - **Block 1 double-space typo** "Open  the filter panel" → "Open
    the Filter Panel" silently fixed.
  - **Block 4 typo** "the molecule is sketched **corectly**" →
    "the molecule is sketched **correctly**" (per task directive).
  - **Block 4 numbering glitch** "4. 4." — the block was labeled
    as block-header "4" then steps numbered 1-2 with one stray "4"
    label drift between Blocks 4 and 5; renumbered cleanly within
    Block 4 as steps 1-3 (the third step being the added implicit
    "no console errors" verification per Block-4 wrap-up).
  - **Block 5 numbering glitch** "4. 4." — original numbered two
    consecutive steps as "4."; renumbered cleanly 1-4 in the
    migrated `### Clone view + cross-view filter sync` sub-scenario.
  - **Block 3 menu phrasing** "mol column's hamburger menu > filter"
    → canonical phrasing "molecule column header → hamburger menu →
    **Filter**" (capitalization + clarity) per atlas
    `chem.search.substructure.filter` / hamburger-menu convention.
  - Precedent: `similarity-search.md` migration silently fixed two
    "Similatity" → "Similarity" typos and `calculate.md` migration
    fixed "asme" → "same" + empty `()` typo — see respective
    migration reports' Decisions § "Source-text defect silent fixes".
- **Scenario constraint extraction (per migration-prompt §
  Scenario constraint extraction).**
  - **(a) FORBIDDEN substitutions:** None explicit in the original
    body (no "must not", "forbidden", "not allowed" phrasing).
    HOWEVER `pyramid_layer: integration` + 8
    `ui_coverage_responsibility` flows per chain rev 2 mandate UI
    driving for the substructure-filter cycle, cell context-menu
    *Use as filter*, hamburger-menu *Filter*, view-clone, and
    cross-view inspection — surfaced explicitly in the migrated
    Notes (`No JS API substitution for UI flows`) as a structured
    behavioral guard.
  - **(b) REQUIRED actions:** All 8 entries of
    `ui_coverage_responsibility` are exercised via UI driving in
    the migrated body — Block 1 covers
    `chem-substructure-disable-filter` + `filter-panel-close-open`
    + `chem-structure-filter-draw` (the draw entry-point alternative)
    + `chem-substructure-use-as-filter` (the alternative entry-point);
    Block 2 covers `chem-substructure-use-as-filter`; Block 3 covers
    `chem-substructure-hamburger-menu-add-filter` +
    `chem-structure-filter-draw`; Block 4 covers
    `chem-substructure-remove-filter` +
    `chem-substructure-use-as-filter`; Block 5 covers `clone-view`
    + `filter-cross-view-sync` + `chem-substructure-use-as-filter`.
    No flow is deferred.
  - **(c) Missing-selector escalation:** Migrator does NOT have
    explicit selector definitions for the Substructure Filter
    widget sketcher, Filter Panel toggle DOM tree (`.grok-icon-filter`
    + `.d4-filter-panel`), cell context-menu *Current value > Use
    as filter*, column header hamburger menu, or the view-tab
    context menu *Clone view* affordance in
    `grok-browser/references/projects.md`. Spec-time selector
    discovery is the Automator's responsibility (per existing
    sibling `structure-filter-spec.ts` patterns +
    `page.evaluate` introspection); a reference-file approval-
    required proposal to add Substructure-Filter-specific selectors
    to `grok-browser/references/widgets/substructure-filter.md` is
    a Phase 2 deliverable — out of this migration's scope.
  - **(d) Reference templates:** Code-style anchors for the
    Automator (cited verbatim, not modified here):
    `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/structure-filter-spec.ts`
    (current existing test — sibling shape with `spec-login` +
    `page.evaluate` + `page.locator` patterns) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/similarity-search-spec.ts`
    (adjacent-folder same-cycle viewer-property-panel walk pattern,
    if helpful) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/scaffold-tree-functions-spec.ts`
    (adjacent-folder same-section viewer + property-panel pattern).
  - **(e) Scenario authority clause:** No prompt-vs-scenario
    conflict surfaced. The original scenario is a 5-block multi-
    invariant walk; the migrated body preserves all 5 blocks as
    `### <name>` sub-scenarios with renumbered steps + explicit
    verifications. No `prompt_scenario_conflict:` field to flag
    for orchestrator.
- **No invented sub_features / helpers.** Every id in
  `sub_features_covered` (`chem.search.substructure`,
  `chem.search.substructure.filter`,
  `chem.search.substructure.top-menu`,
  `chem.search.use-as-filter`, `chem.sketcher`) is present in
  `feature-atlas/chem.yaml` (verified at atlas lines 141-186 for
  `chem.search.substructure.*` ids and `chem.search.use-as-filter`
  at line 180; `chem.sketcher` present in atlas). Every reused
  helper is in `helpers-registry.yaml`; candidate helpers are
  flagged as candidates, NOT referenced by name in the migrated
  body.

## Opt-outs (SCOPE_REDUCTION proposals)

### SR-01: A-STRUCT-02 carryforward (chain-level edge/perf coverage)

**Cited technical dependency:** A-STRUCT-02 ("at least one scenario is
edge case or negative path — i.e. its frontmatter `coverage_type` is
`edge` or `perf`") is a **section-level structural invariant**. This
scenario (`Advanced/structure-filter.md`) is a multi-invariant
substructure-filter cycle walk across five sub-scenarios on a single
viewer/panel surface — `coverage_type: regression` is the natural fit
per chain rev 2 (`pyramid_layer: integration`, classification
`medium`, multi-invariant walk). Forcing `coverage_type: edge` would
mis-classify (the scenario asserts happy-path filter-cycle
transitions across five independent invariants, not a specific
failure-mode invariant).

**Resolution path** (mirrors `info-panels.md`, `calculate.md`, and
`similarity-search.md` 2026-05-11 carryforward precedents in this same
migration cycle, and the `radar.md` 2026-05-07 origin precedent at
decision-log line 8528-8530): A-STRUCT-02 will be satisfied
**chain-wide** by the 10 `bug_focused_candidates[]` already emitted
under `scenario-chains/chem.yaml` rev 2 (each chem-bug-focused spec
naturally carries `coverage_type: edge` because it tests a specific
bug repro path — sister-pattern to charts `radar-save-reopen-bug.md`,
`sunburst-date-column-bug.md`, etc., all of which carry
`coverage_type: edge`). Additionally `r-group-analysis.md` (same-
cycle, same section, `coverage_type: edge`) already organically
satisfies the section-level invariant. **GROK-14028**
(`chem-grok-14028-spec.ts`) is one of the 10 candidates and will
carry `coverage_type: edge` once authored, organically extending
coverage for the substructure-filter UI/data cleanup invariant
adjacent to this scenario's Block 1.

**Why this is a real dependency, not effort:** the per-scenario
edge|perf marker cannot be invented mid-migration without violating
D-SAN-02 (no silent additions). The chain-wide satisfaction path is
already wired (10 candidate specs to be authored downstream by Test
Designer / Automator; plus `r-group-analysis.md` organic edge
coverage). `Advanced/structure-filter.md` must not block the migrate
cycle on a structural property that is the **chain's**
responsibility, not the per-scenario responsibility.

**Scope of waiver:** A-STRUCT-02 verdict on
`Advanced/structure-filter.md` is **deferred to chain-level
evaluation**. Critic A should return SCOPE_REDUCTION (carryforward),
not FAIL, citing this proposal.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through
A-STRUCT-06, A-COVERAGE-*, A-MERIT-*) are expected to PASS on
`Advanced/structure-filter.md` without any other SR proposals. No
step is opted out for effort. The single entry above cites a real
section-level structural property whose satisfaction path is owned
by the chain, not by this per-scenario migration (D-MERIT-01
compliant)._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist
yet. The migrated body's 5 sub-scenarios are realizable against the
current playwright + helpers-registry + atlas state. The candidate
helpers surfaced in Decisions are convenience abstractions that the
Automator may inline pending registration — they do NOT block spec
realization. The existing sibling
`Advanced/structure-filter-spec.ts` already drives the Substructure
Filter via `page.evaluate` + `page.locator` per the sibling-test
pattern; Automator can extend that pattern without a reference
update._

## Edge cases

The original lists no explicit "edge case" keyword. Implicit edge
cases derivable from the scenario content + chain-level
considerations:

- **Filter state survival across panel close + reopen.** Block 1's
  core invariant — disabling the filter + closing the panel +
  re-opening the panel preserves the filter shape (sketcher UI,
  disabled state) so re-enabling restores the original filtered
  row set identically. PRESERVED as scenario steps (Block 1 steps
  2-5). This is a defense-in-depth sister to GROK-14028's Reset-
  cleanup invariant — adjacent surfaces (panel re-render vs.
  Reset action), distinct failure modes.
- **Use as filter cell action with panel closed.** Block 2 contracts
  on the cell context-menu action working uniformly whether the
  Filter Panel is open or closed (panel auto-opens or the filter
  applies behind a closed panel — either is an acceptable
  rendering). PRESERVED as scenario step (Block 2 step 3) with
  the verification bounded to "filter UI correct" (the original's
  phrasing). Atlas does not enumerate the auto-open vs. behind-
  closed-panel rule; the original assertion governs.
- **Hamburger-menu Filter + Add filter combinatorics.** Block 3
  contracts on the hamburger-menu **Filter** action opening a
  sketcher (dialog or in-place widget) for the molecule column,
  and on **Add filter** combining the new structure with the
  existing filter (additive, not replacing). PRESERVED as
  scenario steps (Block 3 steps 3-5). The exact combining
  semantics (AND vs. OR of substructure matches) is not specified
  by the original; the assertion is bounded to "combined filter
  applied + table re-filters". Flag for QA pair review at spec
  realization (Automator picks a verification surface — bitset
  count delta, hit-tile re-render — at spec time).
- **First-on-panel ordering after Remove + Use as filter.** Block 4
  contracts on the newly-added Structure filter appearing first
  on the Filter Panel after the *Remove → Use as filter* sequence.
  PRESERVED as scenario step (Block 4 step 2 bullet 2). The exact
  ordering rule (LIFO? Structure-filter-always-on-top? insertion-
  order?) is not enumerated by the original or atlas; the
  assertion is on the post-condition for this specific sequence.
  Flag for QA pair review.
- **Clone view + cross-view filter sync.** Block 5's core
  invariant — a filter applied on one view of a cloned table-pair
  synchronizes the filtered row set across both views. PRESERVED
  as scenario step (Block 5 step 4). Datagrok's filter state is
  a property of the underlying DataFrame rather than per-view, so
  the sync should be implicit; the scenario explicitly contracts
  on observable cross-view consistency to guard against
  regression if the filter-state ownership semantic ever drifts
  to per-view.
- **"Some structure column is added to the filter panel" (Block 5
  step 1 + step 3).** Both views' Filter Panels must auto-include
  at least one structure column (the molecule column from
  SPGI.csv) — implicit filter-panel-auto-population invariant
  that the scenario contracts on. PRESERVED as verification
  (Block 5 step 1 + step 3 verification clauses). The exact
  auto-population rule (semType=Molecule → add filter? always-on-
  for-Molecule-columns?) is implementation-defined; the assertion
  is bounded to "the structure column is added".
- **Console errors throughout.** Implicit across all 5 blocks:
  each transition (toggle, panel-cycle, cell-action, hamburger-
  menu, view-clone, cross-view inspection) must complete without
  console errors. PRESERVED as a per-step verification weave
  (Blocks 1 step 4-5, Block 2 step 3, Block 3 step 5, Block 4
  step 3, Block 5 step 4 explicitly include "no console errors"
  assertions). Sister with `info-panels.md` /
  `scaffold-tree-functions.md` / `similarity-search.md` patterns.

No edge case is moved to atlas, manual_only, deferred, or a
separate scenario. All implicit edge cases are PRESERVED as
scenario steps OR Notes-level caveats in the migrated body.

## Unresolved ambiguities

- **"Set a filter" (Blocks 2, 3, 4 preamble) — entry-point freedom.**
  Blocks 2 and 3 (and Block 4 implicitly) start with "Set a filter"
  / "Remove Structure filter" without specifying the entry-point
  (draw in sketcher? *Current value > Use as filter*? hamburger
  menu? *Add filter* button?). Block 1 step 1 documents the
  alternatives ("draw OR current value > use as filter"). Migrated
  body preserves the entry-point freedom — each block's Step 1
  reads "set a substructure filter (draw a structure or via
  **Current value > Use as filter**)". Automator at spec time picks
  one entry-point per block — likely the cell-action *Use as
  filter* (smaller selector surface) — or alternates entry-points
  across blocks to broaden coverage. Flag for QA pair review at
  spec realization.
- **"Filter UI is correct" (Block 2 step 3) — verification depth.**
  Block 2's terminating assertion is qualitative ("check that
  filter UI is correct") — silent on which specific UI properties
  to check (sketcher renders the molecule? Filter Panel auto-opens?
  filter is added to the panel? table re-filters?). Migrated body
  resolves by enumerating four concrete sub-checks (panel re-opens
  OR filter applies behind closed panel; sketcher shows molecule
  sketched correctly; table re-filters; no console errors). Flag
  for QA pair review — Automator may add or omit sub-checks at
  spec time.
- **"Draw another structure" (Block 3 step 4) — sketcher widget
  identity.** Block 3 step 3 says "Go to the mol column's hamburger
  menu > filter" — the resulting sketcher could be (a) a modal
  sketcher dialog, (b) an in-place sketcher widget inside the
  Filter Panel, or (c) the column-header-attached sketcher. The
  original is silent on the rendering. Migrated body phrases as
  "a sketcher dialog (or in-place filter sketcher widget) opens
  for the molecule column" — preserves the ambiguity. Automator
  pins to actual rendering at spec time.
- **"Some structure" cell selection (Block 5 step 4) — which row?**
  Block 5 step 4 says "Right-click some structure > Current value >
  Use as filter" — silent on which row's molecule to click.
  Migrated body preserves "right-click a molecule cell" — generic.
  Automator picks a deterministic row (e.g. row 0 / first molecule
  cell in the rendered grid) at spec time. Flag for QA pair review.
- **Block 4 "4. 4." numbering interpretation.** The original Block
  4 has two "4." labels at the block-header level (Block 4 is
  numbered "4" and Block 5 is numbered "5"; internally Block 4
  contains "1. Remove Structure filter" + "2. Current value > Use
  as filter - check:" with sub-bullets — i.e. 2 numbered steps + 2
  sub-bullets, NOT 4 numbered steps). The chain rev 2 directive
  flags "4. 4." duplicate; this is interpreted as Block 4's two
  numbered steps + Block 5's two `4.` labels (Block 5's "4. Open
  filters..." + "4. Right-click..."). The migrated body renumbers
  Block 4 cleanly 1-3 (steps 1-2 + the added implicit "no console
  errors" verification step 3) and Block 5 cleanly 1-4.
  Documented in Step mapping. Flag for QA pair review if the
  interpretation is contested.
