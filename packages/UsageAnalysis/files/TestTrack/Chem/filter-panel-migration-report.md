# Migration Report — filter-panel.md

First-cycle migration of `filter-panel.md` (TestTrack/Chem section) per chain YAML
`scenario-chains/chem.yaml` rev 2 (`output_plan.filter-panel.md`: `target_layer: playwright`,
strategy `simple`, classification `medium`, `pyramid_layer: integration`, `depends_on: []`). Filter
Panel deep-dive walk on `System:DemoFiles/SPGI.csv` across 11 `ui_coverage_responsibility` entries
spanning Filter Panel basics + four settings tabs (Contains / Included in / Exact / Similar) +
**Use as filter** cell action + hamburger-menu Filter route + drag-and-drop column header +
cross-filter with other column types + sketcher↔filter synchronization via hamburger-menu re-edit.

## Step mapping

The original is a 22-line scenario with severe step-numbering chaos in the body — numerals
`1.1.1.2.1.2.1.1.2.2.2.3.4.` interleave across the 11 numbered bullets without semantic grouping
(markdown auto-renumbers display but the source is unintelligible). Renumbered with explicit block
labels per chain rev 2 footer note (c) (Olena 2026-05-11, chain footer note (c) "source-text defects
in scenario bodies — Migrator fixes during migration"); the body parses semantically into 5 blocks
(A — Filter Panel basics; B — Use as filter; C — Hamburger menu Filter; D — Drag-and-drop column
header; E — Cross-filter + sketcher↔filter sync). Typo `sketcherr` → `sketcher`. `Drag'n'drop` →
`Drag-and-drop` (canonical phrasing). Defects fixed silently per chain rev 2 directive.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Line 1 — "Open linked datasets" (chaotic `1.`) | Setup step 1 (provisioning — `System:DemoFiles/SPGI.csv` per JSON footer order 11) + Scenarios step 1 (open dataset as the scenario entry) | preserved (split for clarity — provisioning in Setup; open action in scenario body) |
| Line 2 — "Open the Filter Panel." (chaotic `1.`) | Scenarios > Block A step 2 (Filter Panel toggle + verify panel renders the substructure filter widget) | preserved (expected-result verification added per D-STEP-02 — panel appears + filter widget rendered) |
| Line 3 — "Sketch a substructure in the structure filter." (chaotic `1.`) | Scenarios > Block A step 3 (sketch substructure inside the in-Filter-Panel OCL sketcher + verify table filters + substructure-highlight overlay rendered) | preserved (expected-result verification added per D-STEP-02 — table filters + atoms visually highlighted) |
| Line 4 — "Check all filter settings (Contains, Included in, Exact, Similar,...)." (chaotic `2.`) | Scenarios > Block A step 4 (cycle four settings tabs Contains / Included in / Exact / Similar; per-tab verify re-evaluation + row-set update + no console errors) | preserved (the four tabs enumerated explicitly + each tab's semantic meaning called out + per-tab expected-result verification added per D-STEP-02) |
| Line 5 — "Close the Filter Panel." (chaotic `1.`) | Scenarios > Block B step 5 (close Filter Panel between Block A and Block B) | preserved |
| Line 6 — "On the grid, right-click a molecule and select **Current Value > Use as filter**." (chaotic `2.`) | Scenarios > Block B step 6 (right-click molecule cell → Current Value > Use as filter; verify Filter Panel reopens populated + orientation preserved on SMILES + table filters) | preserved (expected-result verification added per D-STEP-02 — panel reopens + filter populated + orientation preserved per atlas + table filters) |
| Line 7 — "Close the Filter Panel." (chaotic `1.`) | Scenarios > Block C step 7 (close Filter Panel between Block B and Block C) | preserved |
| Line 8 — "Open the column's hamburger menu > Filter." (chaotic `1.`) | Scenarios > Block C step 8 (open column header hamburger menu → Filter; verify sketcher dialog opens) | preserved (expected-result verification added per D-STEP-02 — sketcher dialog opens per atlas `chem.search.substructure.editor`) |
| Line 9 — "Draw another structure and click **Add filter**." (chaotic `2.`) | Scenarios > Block C step 9 (draw structure in dialog + Add filter; verify substructure filter added to Filter Panel + table filters) | preserved (expected-result verification added per D-STEP-02 — filter added + table filters) |
| Line 10 — "Remove the structure filter from the Filter Panel." (chaotic `2.`) | Scenarios > Block C step 10 (remove substructure filter from Filter Panel; verify filter widget removed + table reverts to no-substructure-filter state) | preserved (expected-result verification added per D-STEP-02 — widget removed + table reverts) |
| Line 11 — "Drag'n'drop the Structure column's header to the Filter Panel - filter should be added to the top of the Filter Panel." (chaotic `2.`) | Scenarios > Block D step 11 (drag-and-drop column header onto Filter Panel; verify filter added at the TOP of the panel + empty until query drawn) | preserved (expected-result verification preserved verbatim — "filter should be added to the top of the Filter Panel"; typo-canonicalized `Drag'n'drop` → `Drag-and-drop`) |
| Line 13 — "Check the filtering of chemical columns with any other filter." (chaotic `3.`) | Scenarios > Block E step 12 (apply additional non-molecule filter; verify substructure filter and additional filter compose via logical AND; cross-filter invariant) | preserved (expanded into concrete actionable form — "any other filter" resolved to a numeric range filter or categorical filter on a non-molecule SPGI column, since SPGI.csv has both numeric and string columns alongside the molecule column; expected-result verification added per D-STEP-02 — visible-row set is the intersection of both filters) |
| Line 15 — "Changes in the sketcher through the hamburger menu modify the filter, and it applies the filtering. **In the filter and the sketcherr [sic] we always should observe the same structure.**" (chaotic `4.`) | Scenarios > Block E step 13 (re-enter sketcher dialog via hamburger menu with current filter structure pre-loaded; modify structure; apply; verify filter widget updates + table re-filters + same structure observable in filter widget and sketcher dialog) | preserved (verbatim "we must always observe the same structure" assertion mapped to the filter↔sketcher synchronization invariant per D-STEP-02; typo `sketcherr` → `sketcher` silently fixed) |
| (chaotic line-numbering gap — original has interleaved blank lines between bullets at lines 12 and 14) | n/a | source-text defect (silently fixed — migrated body numbers steps cleanly in five labeled blocks; chain rev 2 footer note (c) "source-text defects" — Migrator fixes for step renumbering) |
| Typo `sketcherr` (line 15) | n/a | source-text typo (silently fixed — `sketcherr` → `sketcher` per chain rev 2 directive) |
| Phrasing `Drag'n'drop` (line 11) | n/a | source-text non-canonical phrasing (silently fixed — `Drag'n'drop` → `Drag-and-drop` per chain rev 2 directive) |
| JSON footer `{ "order": 11, "datasets": ["System:DemoFiles/SPGI.csv"] }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` + migrated Setup step 1 dataset enumeration) |

No original step is silently dropped. The 11 numbered bullets in the original body all land in either
Setup step 1 (Line 1 provisioning half) or Scenarios > Block A-E (Lines 2-15). The chaotic numbering
is documented above per D-SAN-02.

## Decisions

- **Why this `target_layer`:** chose `playwright` per `scenario-chains/chem.yaml`
  `output_plan.filter-panel.md.target_layer = playwright`. The scenario requires DOM-level UI
  driving: Filter Panel toggle on the table view toolbar, in-Filter-Panel sketcher draw via the OCL
  inline canvas, settings-tabs cycle through the four match-mode tabs (gear / context-menu
  interaction), molecule-cell right-click context menu **Current Value > Use as filter**, column
  header hamburger-menu **Filter** entry, sketcher dialog draw + **Add filter** button, substructure
  filter widget remove control, drag-and-drop of the column header onto the Filter Panel,
  cross-filter composition with another column-type filter, and hamburger-menu re-edit propagation
  to the filter widget — all UI behaviors that Playwright drives natively and the apitest layer
  cannot exercise as an interaction walk. The filter↔sketcher synchronization invariant ("we must
  always observe the same structure") specifically requires real DOM observation across two
  widgets (Filter Panel widget + sketcher dialog) that no JS API substitution can mirror. Chain YAML
  reason field: "Filter Panel deep-dive — 4 settings tabs + multiple entry points + sketcher↔
  sketcher sync invariant. Sibling filter-panel-spec.ts already exists at playwright per existing-
  test-index. Filter↔sketcher coordination requires real DOM. Single scenario; simple strategy."
- **Why this `coverage_type`:** chose `regression` per chain rev 2 (`pyramid_layer: integration` —
  multi-entry-point × settings-tabs × cross-filter × filter↔sketcher-sync walk; not `smoke` since
  the section's smoke is `Advanced/scaffold-tree-functions.md` per chain
  `ui_coverage_plan.smoke_scenario`; not `edge` / `perf` since the scenario walks the happy path
  only — no Reset-filter-doesn't-clear-sketcher-input failure-mode invariant exercised here, no
  `MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT` threshold tested). The Reset-cleanup invariant per
  `GROK-14028` is owned by the dedicated bug-focused spec `chem-grok-14028-spec.ts` per chain
  `bug_focused_candidates[]`, not this happy-path walk. `regression` is the natural fit for a multi-
  entry-point × multi-interaction walk that guards against regression-of-the-set across the
  substructure filter surface.
- **Why this `strategy`:** `simple` per chain YAML `output_plan.filter-panel.md.strategy = simple`.
  Single scenario, no cross-file fixture; chain analyzer classified `medium` (multi-entry-point /
  multi-interaction) but kept strategy `simple` — single migrated `.md`, single downstream spec.
  Pattern 1 (decomposition) is NOT triggered (chain classifier returned `medium`, not
  `complex-standalone`; no child scenarios with `produced_from: decomposed`). Pattern 2 (bug-focused
  slice) applies for awareness only — `filter-panel.md` IS a span in three of the chain's 10
  `bug_focused_candidates[]` (`chem-grok-12758-spec.ts` Step 1; `chem-grok-14028-spec.ts` Step 2;
  `chem-grok-17595-spec.ts` Step 1; `chem-github-2942-spec.ts` Step 1), but the spec authoring is
  the Automator's downstream concern, not the Migrator's; `related_bugs` cites `GROK-14028` per
  bug-library directly-affecting overlap (see Bug library bullet below). Pattern 3 (UI-smoke
  consolidation) is NOT triggered (no `ui_consolidation_proposals[]` from Gate F). Standard per-
  scenario migration per "Migrated scenario shape".
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/filter-panel-spec.ts` — existing
    playwright-layer test ("Chem: Filter Panel", category Chem, `features_covered:
    [chem.filter-panel]`, helpers `spec-login`, patterns `test.skip` + `uses-fixture-runners` +
    `uses-grok.dapi` + `uses-grok.shell` + `uses-page.evaluate` + `uses-page.locator` +
    `uses-playwright-test`). Per `existing-test-index.yaml` line 32480. Note the `test.skip` marker
    in the existing spec — Automator at extension time will need to assess whether to remove the
    skip while extending coverage to the full block-A-through-E walk per this migrated scenario.
    Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/chemical-space-migration-report.md` and
    `chemical-space.md` (2026-05-11 predecessor in same migration cycle) — shape anchor for the
    migrated `.md` (multi-block scenario body + `coverage_type: regression` precedent + SR-01
    carryforward block + bug-awareness `related_bugs` precedent). Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/activity-cliffs-migration-report.md` and
    `activity-cliffs.md` (2026-05-11 predecessor in same migration cycle) — shape anchor for the
    step-numbering-gap silent-fix discipline + JSON footer dataset reconciliation pattern + chain-
    rev-2 footer note (c) precedent. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/sketcher.md` and `sketcher-migration-
    report.md` (2026-05-11 predecessor in same migration cycle) — shape anchor for the sketcher-
    based scenario format and `chem.sketcher.*` sub_feature coverage. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts` — sibling Chem
    playwright test. Same shape pattern; confirms the section-wide playwright convention. Read-only.
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`) —
      section-standard fixture login; Automator will use in spec `beforeEach`.
    - `softStep` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`) — per-step
      soft-fail wrapper; Automator wraps each block / step verification in `softStep` for granular
      failure reporting (especially important for a 5-block × multi-step walk where individual
      blocks may surface entry-point-specific dialog branches).
    - `closeAllViews` (`public/packages/UITests/playwright/scripts/helpers.ts:129`) — cleanup at the
      end of the scenario.
  - **Candidate helpers (NOT yet in registry — flagged for addition via helpers-registry curator;
    per migration-prompt Helpers discipline §):**
    - `helpers.playwright.chem.toggleFilterPanel(page)` — toggles the Filter Panel side pane on the
      active table view. Surfaced by Scenarios > Block A step 2, Block B step 5, Block C step 7.
      Convenience abstraction; spec Author can inline the toolbar-icon click until the helper
      lands.
    - `helpers.playwright.chem.drawSubstructureInFilterPanelSketcher(page, smiles)` — sketches a
      substructure inside the in-Filter-Panel OCL sketcher tile (using SMILES → MOLFILE conversion +
      injection through the sketcher's molecular-input field). Surfaced by Scenarios > Block A
      step 3 and Block D step 11 (post-drag-and-drop fill).
    - `helpers.playwright.chem.cycleSubstructureFilterSettingsTabs(filterWidget)` — cycles through
      the four match-mode tabs (Contains / Included in / Exact / Similar) on a substructure filter
      widget. Surfaced by Scenarios > Block A step 4.
    - `helpers.playwright.chem.useAsFilter(grid, rowIndex, columnName)` — right-clicks a molecule
      cell and selects **Current Value > Use as filter**. Surfaced by Scenarios > Block B step 6.
    - `helpers.playwright.chem.openColumnHamburgerFilter(grid, columnName)` — opens the column
      header hamburger menu and selects **Filter** to open the sketcher dialog. Surfaced by
      Scenarios > Block C step 8 and Block E step 13.
    - `helpers.playwright.chem.dragColumnHeaderToFilterPanel(grid, columnName)` — drag-and-drops
      the column header onto the Filter Panel. Surfaced by Scenarios > Block D step 11.
- **Bug library consulted:** yes — `bug-library/chem.yaml` revision 1 (12 curated_bugs, generated
  2026-05-05). Curated bugs whose `affects` intersects this scenario's `sub_features_covered`
  (`chem.search.substructure`, `chem.search.substructure.filter`, `chem.search.substructure.editor`,
  `chem.search.use-as-filter`, `chem.sketcher`, `chem.sketcher.ocl`):
  - **`GROK-14028`** — "Chem | Filter Panel: Reset filter doesn't clear the input line in sketcher".
    `affects: [chem.search.substructure.filter, chem.search.substructure, chem.sketcher,
    chem.actions.copy-smiles]` — intersects 3 of this scenario's 6 sub_features_covered. CITED in
    migrated frontmatter `related_bugs: [GROK-14028]`. Defense-in-depth coverage of the three-layer
    cleanup invariant (UI clear + BitSet null + visual cleanup) is owned by the dedicated bug-
    focused spec `chem-grok-14028-spec.ts` per chain `bug_focused_candidates[]` Step 2 span on
    `filter-panel.md`; this scenario's Notes block calls out the awareness explicitly.
  - **`GROK-12758`** — Scaffold Tree × Sketcher × Substructure Search cross-cutting integration bug.
    `affects` intersects `chem.search.substructure` + `chem.search.substructure.api` + downstream
    Scaffold Tree sub_features. Bug `affects` intersects this scenario's `chem.search.substructure`
    sub_feature, but the cross-cutting span (Scaffold Tree → sketcher → substructure search) is
    NOT exercised by this scenario's body — the cross-cutting invariant is owned by the dedicated
    bug-focused spec `chem-grok-12758-spec.ts` per chain `bug_focused_candidates[]` Step 1 span on
    `filter-panel.md`. NOT cited in `related_bugs` (the bug span attached to this scenario in the
    chain is for the proposed bug-focused spec, not for the happy-path migration; the Scaffold Tree
    cross-coupling is out of scope here).
  - **`GROK-17595`** — Project save+reopen deserialization error for Scaffold Tree blob + Substructure
    Search + datasync. `affects` intersects `chem.search.substructure.top-menu` +
    `chem.search.substructure` + Scaffold Tree sub_features. Bug `affects` intersects this
    scenario's `chem.search.substructure` sub_feature, but the save/reopen + datasync triple
    invariant is NOT exercised by this scenario (no project save/reopen in this scenario's body).
    Owned by `chem-grok-17595-spec.ts` per chain `bug_focused_candidates[]`. NOT cited in
    `related_bugs`.
  - **`GitHub-2942`** (`github 2942`) — CSV export with "Molecules As Smiles" + active filtering
    fails with rowCount mismatch. `affects` intersects `chem.notation.convert-column` +
    `chem.notation.convert-mol` (filter-then-export class of flow); this scenario does NOT exercise
    CSV export. Owned by `chem-github-2942-spec.ts` per chain `bug_focused_candidates[]`. NOT cited
    in `related_bugs`.
  - No other curated bug's `affects` intersects this scenario's sub_features.
- **Decision log queried:** yes — `decision-log.yaml` grepped for `feature: chem` — **zero matches**.
  This is the FIRST CYCLE for the chem section. No prior `migration_decisions`, `layer_decisions`,
  `manual_only`, or `failed_attempts` entries exist for chem. No "approaches off the table"
  constraints apply.
- **Cross-cutting bug citations (chain YAML `bug_focused_candidates[]`).** Per `migration-prompt.md`
  "Cross-cutting bug citations from chain YAML" §: this scenario IS a span in four of the 10
  `bug_focused_candidates[]`:
  - `GROK-12758` → proposed_spec `chem-grok-12758-spec.ts` (Step 1 span).
  - `GROK-14028` → proposed_spec `chem-grok-14028-spec.ts` (Step 2 span).
  - `GROK-17595` → proposed_spec `chem-grok-17595-spec.ts` (Step 1 span).
  - `GitHub-2942` → proposed_spec `chem-github-2942-spec.ts` (Step 1 span).
  Of these, only `GROK-14028` has a directly-overlapping affected sub_feature set with this
  scenario's happy-path coverage (Filter Panel Reset-cleanup gap is the closest behavioral neighbor
  to the substructure-filter happy-path walk); cited in `related_bugs`. The other three are cross-
  cutting / out-of-scope-here invariants whose specs will be authored independently. F-BUG-
  COVERAGE-01 at section-complete is authoritative for cross-cutting bug coverage.
- **UI delegation status.** Per chain YAML `ui_coverage_plan.smoke_scenario:
  Advanced/scaffold-tree-functions.md` (NOT this scenario) and chain
  `ui_coverage_plan.delegated_scenarios[] → filter-panel.md → delegated_to: null` ("UI coverage
  owned (Filter Panel deep-dive: 4 settings tabs + Use as filter + hamburger-menu Filter + drag-
  drop to panel + cross-filter-other-types + sketcher-modify-from-hamburger + filter↔sketcher sync
  invariant). Disjoint from Scaffold Tree smoke. Scenario owns its own coverage."), this scenario's
  11 `ui_coverage_responsibility` flows are owned in-scenario. No SCOPE_REDUCTION proposal
  substitutes JS API for any of the 11 UI flows; therefore the UI delegation discipline of
  migration-prompt § "UI delegation in SCOPE_REDUCTION proposals" does not apply (no SR proposal to
  delegate). Reaffirmed as a hard constraint in the migrated Notes ("No JS API substitution").
- **Scenario constraint extraction (per migration-prompt § Scenario constraint extraction).**
  - **(a) FORBIDDEN substitutions:** None explicit in the original body (no "must not", "forbidden",
    "not allowed" phrasing). HOWEVER `pyramid_layer: integration` per chain + `ui_coverage_delegated_
    to: null` (scenario owns its own coverage) + the filter↔sketcher synchronization invariant ("we
    must always observe the same structure") imply JS API substitution for the 11
    `ui_coverage_responsibility` flows is NOT acceptable — surfaced explicitly in the migrated
    Notes (`No JS API substitution`) as a structured behavioral guard. The synchronization
    invariant cannot be exercised without real DOM observation across two widgets.
  - **(b) REQUIRED actions:** All 11 entries of `ui_coverage_responsibility` are exercised via UI
    driving in the migrated body — Filter Panel toggle (Block A step 2), in-panel sketcher draw
    (Block A step 3), 4 settings-tabs cycle (Block A step 4), Use as filter cell action (Block B
    step 6), hamburger-menu Filter route + Add filter (Block C steps 8-9), filter widget remove
    (Block C step 10), drag-and-drop column header (Block D step 11), cross-filter with other
    column type (Block E step 12), sketcher-modify-from-hamburger + filter↔sketcher sync invariant
    (Block E step 13). No flow is deferred.
  - **(c) Missing-selector escalation:** Migrator does NOT have explicit selector definitions for
    the Filter Panel side-pane root, the in-panel substructure filter widget container, the
    settings tabs (Contains / Included in / Exact / Similar), the **Current Value > Use as filter**
    context-menu entry, the column header hamburger menu **Filter** entry, the substructure filter
    widget remove control, the column header drag handle (for drag-and-drop), or the cross-filter
    composition assertion locator. Spec-time selector discovery is the Automator's responsibility
    (per existing sibling `filter-panel-spec.ts` patterns + `page.evaluate` introspection); a
    reference-file approval-required proposal to add Filter-Panel-specific selectors to
    `grok-browser/references/viewers/filter-panel.md` is a Phase 2 deliverable — out of this
    migration's scope.
  - **(d) Reference templates:** Code-style anchors for the Automator (cited verbatim, not modified
    here): `public/packages/UsageAnalysis/files/TestTrack/Chem/filter-panel-spec.ts` (current
    existing test — sibling shape; `test.skip` marker to reassess) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts` (section-mate
    playwright walk; `spec-login` + `page.evaluate` + `page.locator` patterns) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/sketcher-spec.ts` (in-Filter-Panel /
    in-dialog sketcher manipulation precedent — sketcher draw + commit pattern).
  - **(e) Scenario authority clause:** No prompt-vs-scenario conflict surfaced. The original
    scenario is an 11-bullet (chaotic-numbered) interaction walk; the migrated body preserves all
    original bullets + organizes into 5 labeled blocks. No `prompt_scenario_conflict:` field to
    flag for orchestrator.
- **Source-text fixes silently applied (per chain rev 2 directive footer note (c) — Olena
  2026-05-11).**
  - **Step-numbering chaos resolved with block labels.** Original body's 11 bullets carry numerals
    `1.1.1.2.1.2.1.1.2.2.2.3.4.` interleaved without semantic grouping. The markdown renderer
    auto-renumbers display (so a reader sees 1-11 in the rendered output) but the source is
    unintelligible. Renumbered with explicit block labels (Block A — Filter Panel basics, Block B
    — Use as filter, Block C — Hamburger menu Filter, Block D — Drag-and-drop column header, Block
    E — Cross-filter + sketcher↔filter sync) per chain rev 2 footer note (c) "filter-panel.md (step
    renumber with block labels)". No semantic content moved or dropped; the 5 blocks reflect the
    natural semantic grouping of the original.
  - **Typo `sketcherr` → `sketcher`.** Original line 15 ("In the filter and the sketcherr we always
    should observe the same structure.") has a typo; silently fixed in Block E step 13's verbatim
    assertion preservation. Per chain rev 2 footer note (c).
  - **`Drag'n'drop` → `Drag-and-drop` (canonical phrasing).** Original line 11 uses informal
    `Drag'n'drop`; the canonical phrasing in TestTrack scenario house style (per sibling Chem and
    Projects migrations) is `Drag-and-drop`. Silently canonicalized in Block D step 11. Per chain
    rev 2 footer note (c).
- **No invented sub_features / helpers.** Every id in `sub_features_covered` is present in
  `feature-atlas/chem.yaml` rev 2 (verified per atlas lines 119, 126, 142, 165, 173, 180 for
  `chem.search.substructure` + `chem.search.substructure.filter` + `chem.search.substructure.editor`
  + `chem.search.use-as-filter` + `chem.sketcher` + `chem.sketcher.ocl` — 6 ids total). Every reused
  helper is in `helpers-registry.yaml`; candidate helpers are flagged as candidates, NOT referenced
  by name in the migrated body.

## Opt-outs (SCOPE_REDUCTION proposals)

### SR-01: A-STRUCT-02 carryforward (chain-level edge/perf coverage)

**Cited technical dependency:** A-STRUCT-02 ("at least one scenario is edge case or negative path —
i.e. its frontmatter `coverage_type` is `edge` or `perf`") is a **section-level structural
invariant**. This scenario (`filter-panel.md`) is a multi-entry-point × settings-tabs × cross-filter
× filter↔sketcher-sync walk — `coverage_type: regression` is the natural fit per chain
`pyramid_layer: integration` + the multi-entry-point × multi-interaction shape. Forcing
`coverage_type: edge` would mis-classify (no specific failure-mode invariant being asserted here —
the Reset-filter-doesn't-clear-sketcher-input invariant per `GROK-14028` is owned by the dedicated
bug-focused spec `chem-grok-14028-spec.ts`, not by this happy-path migration).

**Resolution path:** A-STRUCT-02 satisfaction is chain-wide via 10 `bug_focused_candidates[]` in
`scenario-chains/chem.yaml` rev 2 (each chem-bug-focused spec naturally carries `coverage_type: edge`
once authored — `chem-grok-14028-spec.ts` directly targets the filter↔sketcher state-cleanup
invariant and will organically satisfy the section-level invariant in the chem-section bug-focused
authoring phase). Additionally `r-group-analysis.md` (same-cycle, same section, `coverage_type:
edge`) already organically satisfies the section-level invariant. Same SR-01 carried by sibling-cycle
reports `activity-cliffs-migration-report.md`, `calculate-migration-report.md`,
`scaffold-tree-functions-migration-report.md`, `scaffold-tree-migration-report.md`,
`info-panels-migration-report.md`, `Advanced/structure-filter-migration-report.md`,
`Advanced/similarity-search-migration-report.md`, `elemental-analysis-migration-report.md`,
`chemical-space-migration-report.md`, `sketcher-migration-report.md`.

**Scope of waiver:** A-STRUCT-02 verdict on `filter-panel.md` is deferred to chain-level evaluation.
Critic A should return SR (carryforward), not FAIL.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through A-STRUCT-06, A-COVERAGE-*, A-MERIT-*) are
expected to PASS on `filter-panel.md` without any other SR proposals. No step is opted out for
effort. The single entry above cites a real section-level structural property whose satisfaction path
is owned by the chain (and reinforced by the directly-targeting bug-focused spec
`chem-grok-14028-spec.ts`), not by this per-scenario migration (D-MERIT-01 compliant)._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist yet. The migrated body's 13 steps
(Setup step 1 reused as Scenarios step 1 entry + 12 block-organized scenario steps) are realizable
against the current playwright + helpers-registry + atlas state. The candidate helpers surfaced in
Decisions are convenience abstractions that the Automator may inline pending registration — they do
NOT block spec realization. The selector-reference write proposal (Decisions § Scenario constraint
extraction (c)) is a downstream Phase 2 enhancement, not a blocker on this scenario's automation
(existing `filter-panel-spec.ts` already drives the Filter Panel DOM via `page.evaluate` +
`page.locator` per the sibling-test pattern; Automator can extend that pattern without a reference
update). The `test.skip` marker in the existing `filter-panel-spec.ts` is an Automator-side
extension-time assessment, NOT a Migrator-side deferral — this scenario's body is realizable in full;
whether the existing spec's skip should be removed while extending is the Automator's call.

## Edge cases

The original lists no explicit "edge case" keyword. Implicit edge cases derivable from the scenario
content:

- **Four settings tabs (Contains / Included in / Exact / Similar) on the substructure filter
  widget.** Original Line 4 enumerates all four match modes. PRESERVED as scenario step (Block A
  step 4 — explicit per-tab verification of re-evaluation + row-set update + each tab's semantic
  meaning called out). Each tab is its own behavioral branch — Contains is the default forward
  substructure match, Included in is the inverse direction, Exact is structural identity, Similar
  is Tanimoto similarity above threshold.
- **Use as filter from molecule cell context menu — orientation-preservation invariant.** Atlas
  `chem.search.use-as-filter` notes "preserves orientation when source is SMILES". Implicit edge
  case: the filter widget displays the clicked molecule with the same orientation as the source cell
  (not normalized). PRESERVED as scenario step verification (Block B step 6 — explicit "orientation
  preserved when the source is SMILES, per atlas").
- **Drag-and-drop insertion position — top of Filter Panel.** Original Line 11 specifies "filter
  should be added to the top of the Filter Panel" — implicit invariant: drag-and-drop insertion
  position is the TOP of the panel, not the bottom. PRESERVED as scenario step verification (Block
  D step 11 — explicit "filter added at the TOP of the panel"). Distinguishes drag-and-drop from
  other entry routes (Use as filter / hamburger menu Filter — those routes don't specify panel
  insertion position; drag-and-drop is the only route with an explicit position invariant).
- **Cross-filter composition (logical AND with other column-type filter).** Original Line 13
  ("Check the filtering of chemical columns with any other filter") is implicit on which filter type
  to add. PRESERVED as scenario step (Block E step 12 — resolved to "any non-molecule column filter
  e.g. numeric range filter or categorical filter" since SPGI.csv has both column types alongside
  the molecule column; logical-AND composition explicit; toggle invariant explicit). Implicit edge
  case: removing one filter must not affect the other; toggling each filter independently must
  update the row set accordingly.
- **Filter↔sketcher synchronization invariant via hamburger-menu re-edit.** Original Line 15 ("In
  the filter and the sketcher we always should observe the same structure", typo `sketcherr`
  silently fixed) is the critical invariant — the highest-signal step in the scenario. PRESERVED as
  scenario step (Block E step 13 — explicit "we must always observe the same structure" verbatim
  assertion; modify-then-apply-then-verify-three-things sequence). Implicit edge case: the hamburger-
  menu Filter route pre-loads the current filter structure into the sketcher dialog (round-trip
  read), and applying the modified structure propagates back to the filter widget (round-trip
  write). This is the round-trip invariant.
- **Reset-filter cleanup gap (GROK-14028) — three-layer cleanup invariant.** Implicit edge case
  surfaced by `bug-library/chem.yaml` GROK-14028. The scenario does NOT exercise Reset directly —
  flagged in migrated Notes as bug awareness. PRESERVED-AS-NOTE (no scenario step here; defense-in-
  depth coverage owned by `chem-grok-14028-spec.ts` per chain `bug_focused_candidates[]`). Atlas
  critical_path `chem.cp.substructure-search-with-filter` already cites the three-layer cleanup
  invariant ("Reset filter clears all three layers — sketcher UI, BitSet, visual highlighting
  (defense-in-depth per GROK-14028)") — no atlas write required.
- **`MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT` bound.** Atlas `chem.search.substructure`: substructure
  search is bounded by `MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT`. SPGI.csv (≈ 1200 rows) sits below the
  bound; the bound invariant is observable only on a large-volume dataset and is intentionally out
  of scope for this regression walk. PRESERVED as Notes-level caveat. Flagged for atlas curator
  consideration as a candidate `edge_cases` / `manual_only` entry — out of scope here.
- **Console errors throughout.** Implicit across all 11 original bullets: each action (open dataset,
  toggle Filter Panel, sketch, settings tabs, right-click, hamburger menu, drag-and-drop,
  cross-filter) must complete without console errors. PRESERVED as per-step verification weave
  (Block A step 4 + Block E step 12 + Block E step 13 explicitly include "without console errors"
  assertions; the other steps implicitly assume same).

No edge case is moved to atlas, manual_only, deferred, or a separate scenario. All implicit edge
cases are PRESERVED as scenario steps OR Notes-level caveats in the migrated body, with the
row-count bound flagged for atlas curator consideration.

## Unresolved ambiguities

- **"Sketch a substructure" — which specific substructure.** Original Line 3 ("Sketch a substructure
  in the structure filter") is silent on which substructure to draw. Migrated body resolves by
  proposing "any deterministic small fragment — e.g. a benzene ring or pyridine — drawn via the OCL
  sketcher inline canvas". Automator at spec time picks a deterministic SMILES for reproducible
  verification (e.g. `c1ccccc1` benzene → known matching SPGI row count). Flag for QA pair review.
- **"Check all filter settings" — exact verification granularity per tab.** Original Line 4 ("Check
  all filter settings (Contains, Included in, Exact, Similar,...)") with trailing ellipsis is silent
  on whether the four enumerated tabs are exhaustive or whether additional tabs exist. Atlas
  `chem.search.substructure.filter` does NOT enumerate the tabs; the substructure filter widget
  source (`public/packages/Chem/src/package.ts:L337`) is the implementation surface. Migrated body
  enumerates the four named tabs (Contains / Included in / Exact / Similar) and resolves "..." as
  "any additional tabs present in the current implementation are exercised by extension of the
  cycle". Automator at spec time may introspect the tab control's DOM children to discover
  additional tabs (defensive — current implementation may have evolved). Flag for QA pair review.
- **"Use as filter" molecule cell selection.** Original Line 6 ("right-click a molecule and select
  **Current Value > Use as filter**") is silent on which row to right-click. Migrated body resolves
  by "any molecule cell" (semantic — first row is acceptable for determinism). Automator at spec
  time picks `row 0` (or first visible row) for reproducibility. Flag for QA pair review.
- **"Hamburger menu > Filter" exact menu label.** Original Line 8 ("Open the column's hamburger
  menu > Filter") is silent on whether the menu entry label is exactly "Filter" or a variant
  (e.g. "Filter...", "Add filter..."). Atlas `chem.search.substructure.editor` references the
  editor function (`SearchSubstructureEditor`) that "prompts for molecule column when ambiguous";
  the column hamburger menu Filter entry presumably invokes this editor for the column's
  substructure search. Migrated body uses "Filter" as the menu label. Automator at spec time may
  use a text-match selector that accommodates the trailing ellipsis. Flag for QA pair review.
- **"Add filter" button label in sketcher dialog.** Original Line 9 ("click **Add filter**") is
  silent on whether the dialog commits via an "Add filter" button or an "OK" / "Apply" button.
  Migrated body uses "Add filter" as the explicit label, with "OK / equivalent commit control" as
  fallback in Block E step 13's modification-apply step. Automator at spec time picks the actual
  button label by introspection. Flag for QA pair review.
- **"Check the filtering of chemical columns with any other filter" — which "other filter".**
  Original Line 13 is silent on which other filter type to add. SPGI.csv has both numeric columns
  (range filter applicable) and string columns (categorical / search filter applicable). Migrated
  body resolves to "any non-molecule column filter (e.g. add a numeric range filter on a numeric
  SPGI column, or a categorical filter on a string column)". Automator at spec time picks one
  deterministic combination for reproducibility. Flag for QA pair review.
- **"Changes in the sketcher through the hamburger menu modify the filter, and it applies the
  filtering" — modification mechanics.** Original Line 15 is silent on whether the sketcher dialog
  re-opened via the hamburger menu starts BLANK (user re-draws from scratch) or PRE-LOADED with the
  current filter structure (user edits the existing structure). Migrated body resolves to "re-enters
  the sketcher dialog with the current filter structure pre-loaded" based on the round-trip
  expectation implied by the "we must always observe the same structure" assertion (if the dialog
  started blank, the synchronization invariant would be trivially satisfied only after re-drawing).
  This is the more semantically meaningful interpretation. Automator at spec time verifies by
  observation. Flag for QA pair review.
