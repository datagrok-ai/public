# Migration Report — Advanced/scaffold-tree-functions.md

First-cycle migration of `Advanced/scaffold-tree-functions.md`
(TestTrack/Chem section) per chain YAML `scenario-chains/chem.yaml`
rev 2 (`output_plan.Advanced/scaffold-tree-functions.md`:
`target_layer: playwright`, strategy `simple`, classification
`simple`, `pyramid_layer: ui-smoke`, `depends_on: []`,
`ui_coverage_plan.smoke_scenario: Advanced/scaffold-tree-functions.md`).
This scenario is the section's UI smoke and the canonical realization
of atlas critical path `chem.cp.scaffold-tree-add-filter` (p0).

## Step mapping

The original is a 5-step linear scenario walking the Scaffold Tree
viewer create + configure smoke surface on `smiles-50.csv`. Step 2
contains an implicit "empty-tree state opens" verification (bullet
sub-step), and a separate magic-wand sub-step that triggers tree
generation. Step 5's "check them" is an implicit verification (render
+ interact). JSON footer carries `order: 1` + `datasets[1]` —
metadata, not a step.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Step 1 — "Open smiles-50.csv" | Setup step 1 (dataset enumerated) + Scenarios > "Scaffold Tree add + generate + filter + inspect" step 1 (open + verify Molecule auto-detection) | preserved (split for clarity — provisioning declaration in Setup; the open action + auto-detection verification in scenario body) |
| Step 2 main — "Top panel. Run Chem > Analyze > Scaffold tree" | Scenarios > "Scaffold Tree add + generate + filter + inspect" step 2 (top-menu run) | preserved |
| Step 2 sub-bullet — "A view with 'Scaffold Tree is empty' opens" | Scenarios > "Scaffold Tree add + generate + filter + inspect" step 3 (verify viewer added in empty state with placeholder text) | preserved as verification |
| Step 2 sub-bullet — "Press magic wand icon (it's tooltip: Generate from molecular column), the scaffold tree creation should start" | Scenarios > "Scaffold Tree add + generate + filter + inspect" step 4 (press magic-wand) + step 5 (verify generation completes + no console errors) | preserved (split for clarity — action + verification separated) |
| Step 3 — "After the scaffold tree is generated, click the first scaffold - the table should be filtered" | Scenarios > "Scaffold Tree add + generate + filter + inspect" step 6 (click first scaffold node) + step 7 (verify table filter) | preserved (split for clarity — action + verification separated) |
| Step 4 — "On the top of the viewer, check the toolbox of the scaffold tree viewer" | Scenarios > "Scaffold Tree add + generate + filter + inspect" step 8 (inspect toolbox + verify) | preserved |
| Step 5 — "Go to the scaffold tree viewer properties (on the Context Panel) and check them" | Scenarios > "Scaffold Tree add + generate + filter + inspect" step 9 (open property panel + verify properties visible / interactable) | preserved |
| JSON footer `{ "order": 1, "datasets": ["System:AppData/Chem/tests/smiles-50.csv"] }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` + migrated Setup step 1 dataset enumeration) |

No original step or implicit verification bullet is silently dropped.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/chem.yaml` `output_plan.Advanced/scaffold-tree-
  functions.md.target_layer = playwright`. The scenario requires
  DOM-level UI driving: top-menu walk (**Chem > Analyze > Scaffold
  Tree**), magic-wand click inside the empty viewer, scaffold-node
  click on the rendered tree, in-viewer toolbox inspection, and
  Context Panel property-panel inspection — all UI behaviors that
  Playwright drives natively and the apitest layer cannot exercise.
  Chain YAML reason field: "Section's ui-smoke (chosen per Rule 1).
  Sibling `Advanced/scaffold-tree-functions-spec.ts` already exists
  at playwright layer per existing-test-index. Add Viewer + Magic
  Wand + click-to-filter + property-panel inspection — DOM-driven
  flow routed to playwright. Single scenario; simple strategy."
- **Why this `coverage_type`:** chose `smoke`. Per chain
  `pyramid_layer: ui-smoke` — the scenario IS the section's UI
  smoke (`ui_coverage_plan.smoke_scenario` points at this file).
  Canonical create + configure walk for a single viewer, no
  parameter-axis matrix, no bug-repro path, no perf threshold,
  no negative path. `regression` would not fit (single dataset,
  single notation format — not a multi-format walk). `edge` /
  `perf` would not fit (no specific failure-mode or threshold). The
  A-STRUCT-02 edge|perf coverage requirement is a section-level
  invariant; section-level satisfaction is via the chain's 10
  `bug_focused_candidates[]` (see SR-01 in Opt-outs).
- **Why this `strategy`:** `simple` per chain YAML
  `output_plan.Advanced/scaffold-tree-functions.md.strategy =
  simple`. Single scenario, no cross-file fixture; chain analyzer
  classified `simple`, not `complex-standalone`. Pattern 1
  (decomposition) is NOT triggered (5-step linear walk, no
  multi-subsystem split). Pattern 2 (bug-focused slice) is NOT
  triggered for this scenario itself (it is a parent for two
  bug-focused candidates listed in chain `bug_focused_candidates[]`
  — `chem-grok-12758-spec.ts` and `chem-github-3004-spec.ts` — but
  those are downstream Automator concerns; this scenario stays a
  parent smoke). Pattern 3 (UI-smoke consolidation) is NOT
  triggered (no `ui_consolidation_proposals[]` from Gate F).
  Standard per-scenario migration per "Migrated scenario shape".
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/scaffold-tree-functions-spec.ts`
    — existing playwright-layer test ("Chem: Scaffold Tree basic
    functions", category Chem, currently has `test.skip` pattern
    plus `uses-fixture-runners` / `uses-grok.dapi` / `uses-grok.shell`
    / `uses-page.evaluate` / `uses-page.locator` / `uses-playwright-test`;
    helpers_called `spec-login`). Per `existing-test-index.yaml`
    line 32339-32354. Confirms section convention: playwright +
    `spec-login` helper + `page.evaluate` for JS API setup +
    `page.locator` for DOM assertions. The `test.skip` status is
    a downstream Automator concern (re-enabling after migration
    is the Automator's responsibility — out of Migrator scope).
    Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts`
    — sibling Chem playwright test ("Chem: Info Panels on
    smiles.csv"). Same shape pattern; confirms the section-wide
    playwright convention. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels.md`
    + `info-panels-migration-report.md` (2026-05-11 predecessor
    pilot, same section, same migration cycle) — shape anchor for
    the migrated `.md` (frontmatter, Setup / Scenarios / Notes
    order) and the SR-01 A-STRUCT-02 carryforward template adopted
    here. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Charts/sunburst.md`
    — adjacent-section playwright-layer migrated scenario for
    Sunburst viewer create + configure + property-panel walk;
    cross-section parity confirms the viewer-smoke shape (Setup
    + scenario blocks + Notes). Read-only.
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok`
      (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`)
      — section-standard fixture login; Automator will use in spec
      `beforeEach`.
    - `softStep`
      (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`)
      — per-step soft-fail wrapper; Automator wraps each scenario
      verification in `softStep` for granular failure reporting.
    - `closeAllViews`
      (`public/packages/UITests/playwright/scripts/helpers.ts:129`)
      — cleanup at the end of the scenario / between Playwright
      tests.
    - `expandAccordionPane`
      (`public/packages/UITests/playwright/scripts/helpers.ts:141`)
      — applicable in step 9 (expand property-panel accordion
      groups in the Context Panel when the property panel exposes
      grouped knobs).
  - **Candidate helpers (NOT yet in registry — flagged for
    addition via helpers-registry curator; per migration-prompt
    Helpers discipline §):**
    - `helpers.playwright.chem.addScaffoldTreeViewer(page, dataset)`
      — opens dataset + walks **Chem > Analyze > Scaffold Tree**
      top-menu; returns a locator for the empty Scaffold Tree
      viewer. Surfaced by Setup step 1 + Scenarios step 2.
      Convenience abstraction; spec Author can inline the
      top-menu walk until the helper lands.
    - `helpers.playwright.chem.generateScaffoldTreeViaMagicWand(viewer)`
      — clicks the magic-wand icon inside the empty Scaffold Tree
      viewer (tooltip *Generate from molecular column*) and waits
      for the empty-state placeholder to be replaced by the
      rendered tree. Surfaced by Scenarios steps 4-5. Helper
      encapsulates the empty-state → generated-tree transition.
    - `helpers.playwright.chem.clickScaffoldTreeNode(viewer,
      nodeIndex)` — clicks the Nth scaffold node in the rendered
      tree. Surfaced by Scenarios step 6. Companion helper for
      Scenarios step 7 filter verification.
- **Bug library consulted:** yes — `bug-library/chem.yaml`
  revision 1 (12 curated_bugs, generated 2026-05-05). Three bugs
  intersect this scenario's flows but do NOT trigger inclusion in
  `related_bugs` per the chain YAML's frontmatter delta — all
  three are parallel-coverage candidates owned by dedicated
  bug-focused specs in chain `bug_focused_candidates[]`:
  - **GROK-12758** (Scaffold Tree × Sketcher × Substructure-Search
    integration: opening a scaffold node in the Sketcher corrupts
    subsequent `searchSubstructure` state). `affects` intersects
    `chem.analyze.scaffold-tree`, `chem.analyze.scaffold-tree.viewer`,
    `chem.analyze.scaffold-tree.filter` — three sub_features in
    this scenario's `sub_features_covered`. However the bug's repro
    requires (a) opening a scaffold node in the Sketcher and
    (b) a parallel Substructure-Search invocation —
    `scaffold-tree-functions.md` does NEITHER. Chain emits
    `chem-grok-12758-spec.ts` as the bug-focused candidate (spans
    `Advanced/scaffold-tree.md:Step 1` + `filter-panel.md:Step 1`
    + `sketcher.md:Step 2`); awareness only here. Cited in
    migrated Notes section.
  - **GROK-17595** (Project save+reopen deserialization error for
    Scaffold Tree blob with uploaded custom tree + Substructure
    Search + datasync). `affects` intersects
    `chem.analyze.scaffold-tree.add`,
    `chem.analyze.scaffold-tree.generate`, etc.
    `scaffold-tree-functions.md` does NOT save / close / reopen a
    project — the tree blob is in-session only per chain
    `produces`. Chain emits `chem-grok-17595-spec.ts` as the
    bug-focused candidate. Awareness only here. Cited in migrated
    Notes section.
  - **github-3004** (Multi-table Scaffold Tree wrong-binding —
    opening Scaffold Tree from the second table binds to the
    first-opened table). `affects` intersects
    `chem.analyze.scaffold-tree`, `chem.analyze.scaffold-tree.viewer`,
    `chem.analyze.scaffold-tree.add` — three sub_features in this
    scenario's coverage. Chain's `bug_focused_candidates[]`
    explicitly lists `Advanced/scaffold-tree-functions.md:Step 1`
    as a span. However the bug invariant requires multi-table
    state (open table1 → open table2 without closing table1 →
    Scaffold Tree from active table) — `scaffold-tree-functions.md`
    opens only `smiles-50.csv` (one table). Chain emits
    `chem-github-3004-spec.ts` as the bug-focused candidate.
    Awareness only here. Cited in migrated Notes section.
  - **GROK-14028** (Filter Panel Reset cleanup gap) was reviewed
    and rejected for inclusion — `affects` intersects
    `chem.search.substructure.filter`, `chem.sketcher`,
    `chem.actions.copy-smiles`, NONE of which appear in this
    scenario's `sub_features_covered`. No overlap.
  - Remaining curated bugs (GROK-17621 Descriptors Docker timeout,
    GROK-18517 MMP failure, GROK-16329 R-Groups null, GROK-18407
    Chemical Space silent dialog close, GROK-17964 Convert
    Notations duplicate, github-2942 CSV export with filtering,
    GROK-19624 MPO Profile refresh, GROK-16870 RDKit cell renderer
    in non-Chem viewer tooltip) — `affects` disjoint from this
    scenario's `sub_features_covered`. No overlap.
- **Decision log queried:** yes — `decision-log.yaml` grepped
  for `feature: chem` — **zero matches**. This is the FIRST
  CYCLE for the chem section. No prior `migration_decisions`,
  `layer_decisions`, `manual_only`, or `failed_attempts` entries
  exist for chem. No "approaches off the table" constraints apply.
- **Cross-cutting bug citations (chain YAML
  `bug_focused_candidates[]`).** Per `migration-prompt.md`
  "Cross-cutting bug citations from chain YAML" §: this scenario
  is a span (or affected scope) in three candidates — GROK-12758
  (`chem-grok-12758-spec.ts`), GROK-17595
  (`chem-grok-17595-spec.ts`), and github-3004
  (`chem-github-3004-spec.ts`). All three cited in migrated Notes
  section ("Cross-cutting bug awareness — GROK-12758 /
  GROK-17595 / github-3004"). Per migration-prompt this is
  RECOMMENDED, not mandatory; F-BUG-COVERAGE-01 at section-
  complete is authoritative.
- **UI delegation status.** Per chain YAML
  `ui_coverage_plan.smoke_scenario:
  Advanced/scaffold-tree-functions.md` and
  `ui_coverage_responsibility: [chem-add-scaffold-tree,
  chem-scaffold-tree-magic-wand-generate,
  chem-scaffold-tree-node-click-filter, chem-scaffold-tree-toolbox,
  chem-scaffold-tree-property-panel]`,
  `ui_coverage_delegated_to: null` — this scenario IS the section's
  UI smoke (no upstream delegation; it is the delegated-to target
  for sibling scenarios' UI flows where applicable, but its own
  responsibilities are owned). No SCOPE_REDUCTION proposal
  substitutes JS API for any of the 5 UI flows; therefore the UI
  delegation discipline of migration-prompt § "UI delegation in
  SCOPE_REDUCTION proposals" does not apply (no SR proposal to
  delegate). Reaffirmed as a hard constraint in the migrated
  Notes ("No JS API substitution").
- **Scenario constraint extraction (per migration-prompt §
  Scenario constraint extraction).**
  - **(a) FORBIDDEN substitutions:** None explicit in the
    original body (no "must not", "forbidden", "not allowed"
    phrasing). HOWEVER `pyramid_layer: ui-smoke` per chain
    implicitly forbids JS API substitution for the 5
    `ui_coverage_responsibility` flows — surfaced explicitly in
    the migrated Notes (`No JS API substitution`) as a structured
    behavioral guard.
  - **(b) REQUIRED actions:** All 5 entries of
    `ui_coverage_responsibility` are exercised via UI driving in
    the migrated body (add via top menu, magic-wand click,
    scaffold-node click, toolbox inspection, property-panel
    inspection). No flow is deferred.
  - **(c) Missing-selector escalation:** Migrator does NOT have
    explicit selector definitions for `.d4-scaffold-tree-viewer`,
    `.d4-scaffold-tree-viewer .grok-icon-magic-wand`, or the
    in-viewer toolbox / property-panel locators in
    `grok-browser/references/projects.md`. Spec-time selector
    discovery is the Automator's responsibility (per existing
    sibling `scaffold-tree-functions-spec.ts` patterns +
    `page.evaluate` introspection); a reference-file
    approval-required proposal to add Scaffold-Tree-specific
    selectors to `grok-browser/references/widgets/scaffold-tree.md`
    is a Phase 2 deliverable — out of this migration's scope.
  - **(d) Reference templates:** Code-style anchors for the
    Automator (cited verbatim, not modified here):
    `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/scaffold-tree-functions-spec.ts`
    (current existing test, `test.skip` — sibling shape) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts`
    (section-mate playwright walk; `spec-login` + `page.evaluate`
    + `page.locator` patterns) + `public/packages/UsageAnalysis/files/TestTrack/Charts/sunburst-spec.ts`
    (adjacent-section viewer-add + property-panel pattern, if
    helpful).
  - **(e) Scenario authority clause:** No prompt-vs-scenario
    conflict surfaced. The original scenario is a 5-step linear
    smoke walk; the migrated body preserves all 5 with
    sub-bullet splits surfaced as explicit verifications. No
    `prompt_scenario_conflict:` field to flag for orchestrator.
- **No invented sub_features / helpers.** Every id in
  `sub_features_covered` is present in
  `feature-atlas/chem.yaml` rev 2 (verified per atlas lines
  429-463 for scaffold-tree.* ids). Every reused helper is in
  `helpers-registry.yaml`; candidate helpers are flagged as
  candidates, NOT referenced by name in the migrated body.

## Opt-outs (SCOPE_REDUCTION proposals)

### SR-01: A-STRUCT-02 carryforward (chain-level edge/perf coverage)

**Cited technical dependency:** A-STRUCT-02 ("at least one scenario is
edge case or negative path — i.e. its frontmatter `coverage_type` is
`edge` or `perf`") is a **section-level structural invariant**. This
scenario (Advanced/scaffold-tree-functions.md) is a viewer-add+configure
smoke walk — coverage_type: smoke is the natural fit. Forcing
coverage_type: edge would mis-classify.

**Resolution path:** A-STRUCT-02 satisfaction is chain-wide via 10
bug_focused_candidates[] in scenario-chains/chem.yaml rev 2 (each
chem-bug-focused spec naturally carries coverage_type: edge once
authored — sister-pattern to charts radar-save-reopen-bug.md /
sunburst-date-column-bug.md / etc.). Decision-log line 8528-8530
records the radar.md 2026-05-07 carryforward precedent.

**Scope of waiver:** A-STRUCT-02 verdict on
Advanced/scaffold-tree-functions.md is deferred to chain-level
evaluation. Critic A should return SR (carryforward), not FAIL.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through
A-STRUCT-06, A-COVERAGE-*, A-MERIT-*) are expected to PASS on
`Advanced/scaffold-tree-functions.md` without any other SR proposals.
No step is opted out for effort. The single entry above cites a real
section-level structural property whose satisfaction path is owned by
the chain, not by this per-scenario migration (D-MERIT-01 compliant)._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist yet.
The migrated body's 9 steps are realizable against the current
playwright + helpers-registry + atlas state. The candidate helpers
surfaced in Decisions are convenience abstractions that the Automator
may inline pending registration — they do NOT block spec realization.
The selector-reference write proposal (Decisions § Scenario
constraint extraction (c)) is a downstream Phase 2 enhancement, not
a blocker on this scenario's automation (existing
`scaffold-tree-functions-spec.ts` already drives the Scaffold Tree
viewer DOM via `page.evaluate` + `page.locator` per the sibling-test
pattern; Automator can extend that pattern without a reference
update)._

## Edge cases

The original lists no explicit "edge case" keyword. Implicit edge
cases derivable from the scenario content:

- **Empty-state to generated-state transition.** Original Step 2's
  sub-bullet "A view with 'Scaffold Tree is empty' opens" implies
  the viewer starts in an empty state when no tree blob is available
  (i.e. the user has not generated nor uploaded a tree). The
  magic-wand click then transitions the viewer from empty-state to
  generated-state. PRESERVED as scenario step (Scenarios step 3 =
  empty-state verification; Scenarios step 5 = post-generation
  verification — explicit empty-state → generated-state transition
  is exercised in the body).
- **Per-scaffold filter selection (first node click).** Original
  Step 3 says "click the first scaffold - the table should be
  filtered" — implicit "any scaffold node's checkbox / click
  filters by that scaffold" invariant. PRESERVED as scenario step
  (Scenarios step 6 + 7 — first-node click + filter verification).
  The broader "ANY scaffold node click filters by THAT node"
  invariant is implicitly tested but bounded to "first node" by the
  original wording — the scenario tests one node click, not all
  nodes. Atlas / future regression depth may want a "click multiple
  nodes" extension; OUT OF SCOPE for this smoke. Atlas curator note
  (no atlas write here): an atlas-level `edge_cases` entry for
  multi-node-click cumulative-filter behavior could be a candidate
  addition (e.g. "checkbox selection on scaffold tree nodes
  accumulates filters AND vs OR"); flag for atlas curator.
- **Toolbox surface vs. application-level Toolbox.** Original Step
  4 says "check the toolbox of the scaffold tree viewer" — ambiguous
  with the application-level **Toolbox** sidebar (`.d4-toolbox`).
  RESOLVED in migrated body as the in-viewer toolbar surface (top of
  the viewer), not the sidebar. Documented explicitly in migrated
  Notes ("Toolbox surface"). Implementation surface is in
  `Chem/src/widgets/scaffold-tree.ts`. PRESERVED as scenario step
  (Scenarios step 8); resolution documented; flag in Unresolved
  ambiguities for QA pair review.
- **Property panel content (Step 5 implicit verification).**
  Original Step 5 says "check them" (the properties) — implicit
  verification that the properties are present + renderable +
  interactable. PRESERVED as scenario step (Scenarios step 9 —
  property-panel renders + property groups expand + each control is
  visible / interactable). Specific property names / value ranges
  are NOT enumerated by the original — the verification is bounded
  to "renders + interactable", not "values equal X". This bounds
  the spec's assertion surface and keeps the smoke fast.
- **Console errors throughout.** Implicit across all 5 original
  steps: each action ("Open", "Run top menu", "Press magic wand",
  "click", "go to properties") must complete without console errors.
  PRESERVED as a per-step verification weave (Scenarios steps 5, 7,
  8, 9 explicitly include "no console errors" assertions). Sister
  with `info-panels.md` pattern.

No edge case is moved to atlas, manual_only, deferred, or a separate
scenario. All implicit edge cases are PRESERVED as scenario steps OR
Notes-level caveats in the migrated body.

## Unresolved ambiguities

- **"toolbox of the scaffold tree viewer" — exact icon / action set
  to assert.** Original Step 4 is silent on what the toolbox
  contains. The viewer toolbox is implementation-defined in
  `Chem/src/widgets/scaffold-tree.ts` (atlas
  `chem.analyze.scaffold-tree.viewer`) and may evolve. Migrated body
  resolves by asserting "toolbox renders without errors and its
  action set is visible" — non-enumerative. Automator at spec time
  may choose to enumerate specific actions (e.g. "Generate", "Add
  Scaffold", "Clone", "Settings") for higher-fidelity coverage, OR
  keep the assertion bounded to "renders + visible" — spec-time
  decision. Flag for QA pair review at spec realization. (See
  Edge cases section "Toolbox surface vs. application-level
  Toolbox" for the disambiguation in the migrated body.)
- **"first scaffold" — first by DOM-render order or first by
  generation-order / scaffold-tree-root.** Original Step 3 says
  "click the first scaffold" without specifying ordering. The
  Scaffold Tree viewer presents scaffolds in a hierarchy; "first"
  could mean (a) the topmost-rendered scaffold in DOM order
  (positional), (b) the scaffold tree root / Generation 0
  scaffold (semantic), or (c) the first scaffold in the
  generated-blob array order (algorithmic). Migrated body
  preserves "first scaffold node in the rendered tree" — DOM-order
  interpretation. Automator at spec time picks a deterministic
  selector (positional or by-text); if multiple scaffolds are
  equivalent the first-in-DOM-order is acceptable. Flag for QA
  pair review.
- **"check them" (Step 5 property-panel verification depth).**
  Original asks to "go to the scaffold tree viewer properties (on
  the Context Panel) and check them" — silent on which properties,
  expected values, or specific interactions. Migrated body resolves
  by asserting "property panel renders, property groups expand,
  each property control is visible / interactable" — non-specific.
  Automator at spec time may pick representative properties for
  per-knob verification (e.g. assert `Title` editable, assert
  `Ring cutoff` slider visible) but the original is bounded to
  "check them" / visibility. Flag for QA pair review at spec
  realization.
- **Magic-wand icon DOM locator.** Migrated body says "press the
  magic-wand icon (tooltip *Generate from molecular column*)" —
  the canonical discovery is by tooltip text. Datagrok uses
  `Icons.fa('magic')` / `Icons.fa('wand-magic-sparkles')` (per
  MEMORY.md ~Datagrok UI Selector Notes: icons render as
  `.grok-icon i.fal.fa-{name}`). Spec-time selector may use
  `button[aria-label*="Generate from molecular column"]` or
  `.grok-icon-magic` or `[role="button"][title*="Generate"]`
  depending on actual implementation. The reference doc
  `grok-browser/references/widgets/scaffold-tree.md` does NOT
  currently exist (per Decisions § Scenario constraint extraction
  (c) — missing-selector escalation flagged for Phase 2 enhancement).
  Migrator preserves the tooltip-text discovery as the contract;
  spec-time selector adapts.
- **Scaffold Tree viewer empty-state placeholder text.** Migrated
  body says "Scaffold Tree is empty" (exact phrasing from original
  Step 2 sub-bullet) — Automator should verify against the actual
  placeholder text in
  `Chem/src/widgets/scaffold-tree.ts`. Localization variants are
  unlikely (Datagrok is English-only for the smoke surface) but
  the assertion should tolerate minor punctuation / casing
  drift. Flag for QA pair review.
