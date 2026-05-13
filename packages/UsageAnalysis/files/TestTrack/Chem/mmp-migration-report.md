# Migration Report — mmp.md

First-cycle migration of `mmp.md` (TestTrack/Chem section, `order: 14` — LAST scenario of
the Chem section per chain YAML `order_from_files`) per chain YAML
`scenario-chains/chem.yaml` rev 2 (`output_plan.mmp.md`: `target_layer: playwright`,
strategy `simple`, classification `simple`, `pyramid_layer: bug-focused`,
`depends_on: []`). End-to-end Matched Molecular Pairs analysis walk on the bundled
`System:DemoFiles/chem/mmp_demo.csv` demo dataset covering six chain
`ui_coverage_responsibility` entries spanning Chem > Analyze > MMP top-menu entry,
MMPEditor dialog activity-column selection, and all four MMP viewer tabs (Transformation,
Fragments, Cliffs, Generation). **Scenario body IS the GROK-18517 reproduction path
verbatim** — bug-focused regression guard against the minified JS runtime error
("`J.aS(...).b7 is not a function`") that broke MMP generation on the bundled demo
dataset.

## Step mapping

The original is a 17-line scenario (6 numbered steps + sub-bullet on step 3 +
metadata-not-step JSON footer) with one source-text defect: original steps 4 and 5
reference "all three activities" but the bundled `mmp_demo.csv` (read from
`public/packages/Chem/files/demo_files/mmp_demo.csv`) has only TWO activity columns
(`CYP3A4`, `hERG_pIC50`). Step 3 correctly references "two activities". This defect is
silently fixed per chain rev 2 directive footer note (c) (Olena 2026-05-11). The body
otherwise parses semantically into a single linear 7-step walk (open → top menu → editor
+ OK → Transformation tab → Fragments tab → Cliffs tab → Generation tab).

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Line 1 — "Open linked dataset (DemoFiles > chem > mmp_demo)" | Setup step 1 (provisioning — `System:DemoFiles/chem/mmp_demo.csv` per JSON footer `order: 14`) + Scenarios step 1 (open dataset as the scenario entry; verify table view renders + 4 columns + numeric activity columns) | preserved (split for clarity — provisioning in Setup; open action in scenario body; expected-result verification added per D-STEP-02 — table renders + 4 columns enumerated + RDKit cell renderer + 2 numeric activity columns confirmed) |
| Line 2 — "Run top menu Chem -> Analyze -> Matched Molecular Pairs" | Scenarios step 2 (top-menu walk → MMPEditor dialog opens; verify dialog renders without console errors) | preserved (expected-result verification added per D-STEP-02 — MMPEditor surfaces per atlas `chem.analyze.mmp.editor` line 421) |
| Line 3 — "Select all two activities, press 'OK'" | Scenarios step 3 (MMPEditor dialog: select both `CYP3A4` and `hERG_pIC50` activity columns + click OK; **GROK-18517 discriminator inline**: no minified-runtime-error fires; MMP viewer renders) | preserved + AUGMENTED (the GROK-18517 discriminator paragraph is the primary bug-focused invariant of this scenario; expected-result verification added per D-STEP-02 — viewer renders per atlas `chem.analyze.mmp.viewer` line 414; explicit "two activities" framing surfaces the exact dataset shape after silently fixing the inconsistency with steps 4-5) |
| Line 6 (sub-bullet on Line 3 in original — "click first 7 entries and observe changes in Transformation tab") | Scenarios step 4 (click first 7 entries in Transformation tab; verify per-row update of transformation pair + activity-difference values across both activity columns; no console errors per click) | preserved (split as standalone step per migrated 7-step walk; "first 7 entries" preserved verbatim as specific walk depth; expected-result verification added per D-STEP-02 — per-row transformation pair update + both activity-difference values surfaced) |
| Line 7 — "Go to Fragments tab and make sure all three activities are shown" | Scenarios step 5 (Fragments tab; verify both activity columns shown — `CYP3A4` and `hERG_pIC50`; no console errors) | **silently fixed** ("three" → "both"; only 2 activities exist in `mmp_demo.csv`) + preserved (expected-result verification added per D-STEP-02 — fragments grid renders + both columns surface in summary) |
| Line 8 — "Go to cliffs tab and make sure filters for all three activities work" | Scenarios step 6 (Cliffs tab; verify filter controls for both activity columns present + functional + cliffs scatter renders; no console errors) | **silently fixed** ("three" → "both"; only 2 activities exist in `mmp_demo.csv`) + preserved (expected-result verification added per D-STEP-02 — filter controls per column + cliffs scatter visualization renders) |
| Line 9 — "Go to Generation tab and check if results are there" | Scenarios step 7 (Generation tab; verify generation grid populates with candidate molecule variants; no console errors) | preserved (expected-result verification added per D-STEP-02 — generation grid populates + candidate variants surface) |
| JSON footer `{ "order": 14, "datasets": ["System:DemoFiles/chem/mmp_demo.csv"] }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` + migrated Setup step 1 dataset enumeration) |

No original step is silently dropped. All 6 numbered bullets + sub-bullet land in either
Setup step 1 (Line 1 provisioning half) or Scenarios steps 1-7. One source-text defect is
silently fixed with full disclosure (see § Source-text fixes below) per D-SAN-02.

## Decisions

- **Why this `target_layer`:** chose `playwright` per `scenario-chains/chem.yaml`
  `output_plan.mmp.md.target_layer = playwright`. The scenario requires DOM-level UI
  driving: top-menu **Chem > Analyze > Matched Molecular Pairs...** entry click,
  MMPEditor dialog interaction (activity-column multi-selector + OK button), and four MMP
  viewer tab clicks with per-tab content verification (Transformation entry clicks,
  Fragments / Cliffs / Generation tab introspection) — all UI behaviors that Playwright
  drives natively. The GROK-18517 discriminator (no minified-runtime error during MMP
  generation) is observable only through the editor-driven flow; direct `mmpAnalysis`
  JS API invocation bypasses the editor surface where the bug reproduced. Chain YAML
  reason field: "Bug-focused scenario (GROK-18517 mmp_demo reproduction) layered on
  canonical MMP viewer walk (4 tabs). Sibling mmp-spec.ts already exists at playwright
  per existing-test-index. MMP viewer tabs + Transformation entry-click + Fragments/
  Cliffs/Generation DOM interaction routed to playwright. Single scenario; simple
  strategy."
- **Why this `coverage_type`:** chose `edge` per chain rev 2
  (`pyramid_layer: bug-focused` — scenario body IS GROK-18517 reproduction;
  `classification: simple` — single linear walk, no matrix). Bug-focused regression
  guard against a specific failure mode (minified runtime error tree-shook out of the
  MMP code path on the bundled demo dataset). `coverage_type: edge` naturally fits
  bug-focused migrations per migration-prompt § A-STRUCT-02 ("edge case or negative
  path") — GROK-18517's failure mode IS the negative path being guarded against. Not
  `regression` (the natural label for happy-path cross-subsystem walks like
  `chemprop.md` / `chemical-space.md`); the bug-focused intent dominates the canonical-
  walkthrough framing per chain YAML `mmp.md.notes`: "Rule 3 wins over Rule 1 (the bug-
  repro intent dominates the canonical-walkthrough framing) and over Rule 4 (single-
  subsystem — MMP viewer only)." Not `smoke` (the section's smoke is
  `Advanced/scaffold-tree-functions.md`); not `perf` (no perf invariant asserted).
- **Why this `strategy`:** `simple` per chain YAML `output_plan.mmp.md.strategy =
  simple`. Single scenario, single dataset, no cross-file fixture; chain analyzer
  classified `simple` (single linear walk on `mmp_demo.csv`, no matrix axes, no
  decomposition proposal). Pattern 1 (decomposition) is NOT triggered (chain classifier
  returned `simple`, not `complex-standalone`; no child scenarios with `produced_from:
  decomposed`). **Pattern 2 (bug-focused slice) IS the dominant pattern here** —
  `mmp.md` is implicitly the bug-focused slice for GROK-18517 per chain YAML's
  bug_coverage_status comment: "GROK-18517 → mmp.md (entire scenario IS the mmp_demo
  repro path: Open DemoFiles>chem>mmp_demo → Chem | Analyze | MMP → select activities
  → OK → exercise Transformation/Fragments/Cliffs/Generation tabs)". This means
  GROK-18517 is NOT emitted under chain `bug_focused_candidates[]` (only 10 of 12
  curated_bugs require dedicated specs — GROK-18517 and GROK-16329 are covered by
  existing scenario reproduction paths). Pattern 3 (UI-smoke consolidation) is NOT
  triggered. Standard per-scenario migration with bug-focused framing per chain rev 2.
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/Chem/src/tests/mmpa-tests.ts` — Chem package test family for MMPA
    (per existing-test-index.yaml lines 1044+, 3286+, 5528+, multiple references — the
    canonical MMPA test surface; `mmpa` + `mmpaOpens` test entries). Read-only.
  - `public/packages/Chem/src/analysis/molecular-matched-pairs/mmp-function-editor.ts` —
    `MMPEditor` source; activity columns are a `column_list` input per atlas
    `chem.analyze.mmp.editor` line 421. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/chemprop-migration-report.md`
    and `chemprop.md` (2026-05-11 predecessor in same migration cycle) — shape anchor
    for the migrated `.md` (single-scenario body + structured Setup / Scenarios /
    Notes / cross-step invariants). Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/r-group-analysis-migration-report.md`
    (2026-05-11 cycle) — shape anchor for `coverage_type: edge` bug-focused migration
    (R-Group's Block 1 covers GROK-16329 by exercising the no-rgroups balloon as an
    in-line discriminator; mirrors this scenario's GROK-18517 discriminator pattern).
    Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/activity-cliffs-migration-report.md`
    (2026-05-11 cycle) — shape anchor for multi-tab Chem analyze viewer walk
    (Activity Cliffs uses a similar scatter + tabs pattern). Read-only.
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`)
      — section-standard fixture login; Automator will use in spec `beforeEach`.
    - `softStep` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`) —
      per-step soft-fail wrapper; Automator wraps each of the 7 walk steps in
      `softStep` for granular failure reporting (step 3 — MMP generation — is the
      primary discriminator step where soft-fail visibility matters most).
    - `closeAllViews` (`public/packages/UITests/playwright/scripts/helpers.ts:129`) —
      cleanup at the end of the scenario.
  - **Candidate helpers (NOT yet in registry — flagged for addition):**
    - `helpers.playwright.chem.openMmpEditor(page)` — top-menu walk **Chem > Analyze >
      Matched Molecular Pairs...** and waits for the MMPEditor dialog to render.
      Surfaced by Scenarios step 2.
    - `helpers.playwright.chem.selectMmpActivityColumns(page, columnNames)` — multi-
      select activity columns in MMPEditor and click OK. Surfaced by Scenarios step 3.
    - `helpers.playwright.chem.switchMmpViewerTab(page, tabName)` — switch MMP viewer
      tab to one of Transformation / Fragments / Cliffs / Generation. Surfaced by
      Scenarios steps 4-7.
- **Bug library consulted:** yes — `bug-library/chem.yaml` revision 1 (12
  curated_bugs, generated 2026-05-05). `GROK-18517` (priority `p2`, status `fixed`)
  directly intersects this scenario — `affects` list
  `[chem.analyze.mmp, chem.analyze.mmp.top-menu, chem.analyze.mmp.editor,
  chem.demos.mmpa]` is a complete overlap with this scenario's
  `sub_features_covered`. `related_bugs: [GROK-18517]` in frontmatter. The
  reproduction text from the bug entry is preserved verbatim as the scenario body
  (the bug entry steps 1-3 ARE the scenario's steps 1-3) plus the canonical viewer-
  tab walk (steps 4-7). Per chain YAML bug_coverage_status comment: GROK-18517 is
  NOT emitted under `bug_focused_candidates[]` because this scenario already covers
  it organically.
- **Decision log queried:** yes — `decision-log.yaml` grepped for `feature: chem` —
  zero per-scenario `mmp`-specific entries. The chem section is in its first
  migration cycle. No prior `migration_decisions`, `layer_decisions`, `manual_only`,
  or `failed_attempts` entries exist for mmp. No "approaches off the table"
  constraints apply.
- **Cross-cutting bug citations (chain YAML `bug_focused_candidates[]`).** Per
  `migration-prompt.md` "Cross-cutting bug citations from chain YAML" §: `mmp.md` is
  NOT a span in any of the chain's 10 `bug_focused_candidates[]` (verified via grep
  on the chain YAML — mmp does not appear in any `spans:` list). `mmp.md` itself is
  the bug-focused slice for GROK-18517 (per chain bug_coverage_status comment) — but
  GROK-18517 is omitted from `bug_focused_candidates[]` precisely BECAUSE `mmp.md`
  organically covers it. No cross-cutting bug citation applicable; this scenario IS
  the bug coverage.
- **Atlas critical-path adjacency.** This scenario directly covers atlas critical
  path `chem.cp.mmp-analysis` (p1) — sub-features used:
  `[chem.analyze.mmp, chem.analyze.mmp.top-menu, chem.analyze.mmp.viewer,
  chem.analyze.mmp.editor, chem.demos.mmpa]` — 5 of 5 atlas critical-path
  sub-features. Atlas description (line 1493-1496): "OK runs analysis → MMP viewer
  renders fragments / generations / cliffs tabs without minified-runtime errors
  (GROK-18517 regression surface). Demo dataset `mmp_demo` is the golden run." This
  scenario IS the critical-path walk verbatim.
- **UI delegation status.** Per chain YAML
  `ui_coverage_plan.delegated_scenarios[] → mmp.md → delegated_to: null` ("UI
  coverage owned (MMP editor select-activities + viewer Transformation / Fragments
  / Cliffs / Generation tabs). Bug-focused (GROK-18517 mmp_demo reproduction).
  Disjoint from Scaffold Tree smoke. Scenario owns its own coverage."), this
  scenario's 6 `ui_coverage_responsibility` flows are owned in-scenario. No
  SCOPE_REDUCTION proposal substitutes JS API for any of the 6 UI flows; therefore
  the UI delegation discipline does not apply (no SR proposal to delegate).
  Reaffirmed as a hard constraint in the migrated Notes ("No JS API substitution").
- **Scenario constraint extraction (per migration-prompt § Scenario constraint
  extraction).**
  - **(a) FORBIDDEN substitutions:** None explicit in the original body (no "must
    not" phrasing). However `pyramid_layer: bug-focused` + the GROK-18517
    discriminator invariant (the bug reproduced through editor-driven flow with
    minified runtime error) imply JS API substitution for the MMPEditor /
    `mmpAnalysis` flow is NOT acceptable — the bug's reproduction surface is the
    editor-driven flow; direct `mmpAnalysis` invocation bypasses the editor and
    misses the reproduction surface. Surfaced explicitly in the migrated Notes
    (`No JS API substitution`).
  - **(b) REQUIRED actions:** All 6 entries of `ui_coverage_responsibility` are
    exercised via UI driving in the migrated body — Chem > Analyze > MMP top-menu
    walk (Scenarios step 2), MMPEditor activity-column selection + OK (Scenarios
    step 3), Transformation tab walk (Scenarios step 4), Fragments tab walk
    (Scenarios step 5), Cliffs tab walk (Scenarios step 6), Generation tab walk
    (Scenarios step 7). No flow is deferred.
  - **(c) Missing-selector escalation:** Migrator does NOT have explicit selector
    definitions for the MMPEditor activity-column multi-selector, OK button, MMP
    viewer tab strip, or per-tab content surfaces. Sibling `mmpa-tests.ts` drives
    the `mmpa` / `mmpaOpens` flows via package test infrastructure (NOT spec-time
    selectors). The Automator will at spec time discover selectors via `page.
    evaluate` introspection on the open MMPEditor dialog and MMP viewer tab strip.
    Reference-file approval-required proposal to add MMP-specific selectors to
    `grok-browser/references/mmp-viewer.md` (or equivalent) is a Phase 2
    deliverable — out of this migration's scope.
  - **(d) Reference templates:** Code-style anchors for the Automator (cited
    verbatim, not modified here): `public/packages/Chem/src/tests/mmpa-tests.ts`
    (Chem package MMPA test surface — `mmpa` / `mmpaOpens` entries) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/mmp-spec.ts` (if exists
    per existing-test-index) + sibling Chem-section playwright walks (e.g.
    `chemical-space-spec.ts`, `activity-cliffs-spec.ts` from same migration cycle).
  - **(e) Scenario authority clause:** No prompt-vs-scenario conflict surfaced.
    The original scenario has the "two vs three activities" internal inconsistency
    which is silently fixed per chain rev 2 (the dataset has 2 activities, period).
    No `prompt_scenario_conflict:` field to flag for orchestrator.
- **Source-text fixes silently applied (per chain rev 2 directive footer note (c) —
  Olena 2026-05-11).**
  - **"three activities" → "both activities" / "two activities"** in original
    steps 4 and 5. The bundled `mmp_demo.csv` has exactly TWO activity columns
    (`CYP3A4`, `hERG_pIC50`) — verified by reading the CSV header at
    `public/packages/Chem/files/demo_files/mmp_demo.csv` line 1
    (`smiles,CMPD_CHEMBLID,CYP3A4,hERG_pIC50`). Original step 3 correctly references
    "two activities". The "three" in steps 4 and 5 is an authoring defect (likely
    counting the molecules column or `CMPD_CHEMBLID` as an "activity" — but neither
    is selectable in MMPEditor as an activity column). Silent fix preserves the
    scenario's clear intent (verify both activity columns show in Fragments and
    Cliffs tabs) while correcting the dataset-shape mismatch. Full disclosure here
    per D-SAN-02; this resolution is the inconsistency-resolution item flagged in
    the prompt.
- **No invented sub_features / helpers.** Every id in `sub_features_covered` is
  present in `feature-atlas/chem.yaml` rev 2 (verified per atlas lines 399, 406,
  414, 421, 1245 for `chem.analyze.mmp`, `chem.analyze.mmp.top-menu`,
  `chem.analyze.mmp.viewer`, `chem.analyze.mmp.editor`, `chem.demos.mmpa` — 5 ids
  total). Every reused helper is in `helpers-registry.yaml`; candidate helpers
  are flagged as candidates, NOT referenced by name in the migrated body.

## Opt-outs (SCOPE_REDUCTION proposals)

### (none)

This scenario carries `coverage_type: edge` (bug-focused per chain
`pyramid_layer: bug-focused`). It **organically satisfies A-STRUCT-02** ("at least one
scenario is edge case or negative path — i.e. its frontmatter `coverage_type` is `edge`
or `perf`") for the chem section — no SR-01 carryforward is required from this
migration. The GROK-18517 discriminator (no minified-runtime error during MMP
generation) IS the negative-path invariant being guarded against, satisfying the
section-level structural property naturally.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through A-STRUCT-06, A-COVERAGE-*,
A-MERIT-*) are expected to PASS on `mmp.md` without any SR proposals. No step is opted
out for effort. D-MERIT-01 compliant._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist yet. The migrated
body's 7 steps (Setup step 1 reused as Scenarios step 1 entry + 6 walk steps from
original 6 numbered bullets + sub-bullet split) are realizable against the current
playwright + helpers-registry + atlas state. The candidate helpers surfaced in
Decisions are convenience abstractions that the Automator may inline pending
registration — they do NOT block spec realization. The selector-reference write
proposal (Decisions § Scenario constraint extraction (c)) is a downstream Phase 2
enhancement, not a blocker on this scenario's automation._

## Edge cases

The original lists no explicit "edge case" keyword, but the entire scenario IS
itself an edge case (bug-focused regression guard against GROK-18517). Implicit
edge cases derivable from the scenario content:

- **Minified runtime error reproduction.** GROK-18517 surfaces as
  "`J.aS(...).b7 is not a function`" — a tree-shook-function signature from
  minified production build. PRESERVED as scenario discriminator at Scenarios step
  3 (the GROK-18517 discriminator paragraph). The Automator at spec time
  introspects `page.on('pageerror')` / `page.on('console')` for minified-function-
  signature errors during MMP generation. This is the primary verifiable invariant.
- **Demo dataset criticality.** Bundled demo datasets are the primary onboarding
  path for Chem features. PRESERVED-AS-NOTE in migrated Notes ("Bundled demo
  dataset criticality") — the regression-here-is-high-visibility framing
  reinforces the GROK-18517 discriminator's importance.
- **Both activities propagation across tabs.** Both `CYP3A4` and `hERG_pIC50`
  selected at editor time must appear in the Fragments tab summary and the
  Cliffs tab filter controls. One activity column dropping out of either tab is
  a coverage gap regression. PRESERVED as cross-step invariant (migrated body
  Implicit cross-step invariants).
- **Console errors throughout.** Implicit across all 7 walk steps: each action
  (open dataset, top-menu walk, dialog OK, four tab switches + entry clicks) must
  complete without console errors. The minified runtime error of GROK-18517
  surfaces in the JS console — console-error introspection is the most sensitive
  detector. PRESERVED as per-step verification weave (Scenarios steps 3-7
  explicitly include "no console errors" assertions; steps 1-2 implicitly assume
  same).
- **Single-dataset coverage.** Original is single-dataset (`mmp_demo.csv` only).
  PRESERVED as single-dataset migration body (no multi-format extension applied
  — distinct from `activity-cliffs.md` / `chemical-space.md` / `calculate.md` /
  `elemental-analysis.md` multi-format extensions which were driven by source-
  body prose). mmp's original body has no multi-format prose; the scenario's
  source defines the dataset surface as `mmp_demo.csv` only (and this is
  intentional — `mmp_demo.csv` IS the GROK-18517 reproduction dataset).

No edge case is moved to atlas, manual_only, deferred, or a separate scenario. All
implicit edge cases are PRESERVED as scenario steps, cross-step invariants, OR
Notes-level caveats in the migrated body.

## Unresolved ambiguities

- **"two vs three activities" inconsistency in source body.** Original step 3
  says "two activities" — correct. Original steps 4 and 5 say "all three
  activities" — incorrect (only 2 activities in `mmp_demo.csv`). Resolved by
  reading the dataset CSV header and silently fixing "three" → "both" /
  "two" per the actual dataset shape (see § Source-text fixes). No QA pair
  review required — dataset is concrete.
- **"first 7 entries" in Transformation tab — what is an "entry"?** Original
  step 3 sub-bullet says "click first 7 entries and observe changes in
  Transformation tab". An "entry" most naturally maps to a row in the
  transformation grid (one MMP pair per row, surfaced via the MMP viewer's
  transformation tab). Migrated body resolves by treating each "entry" as a
  row click in the transformation grid surface. Automator at spec time
  picks the actual grid-row selector by introspection. Flag for QA pair
  review.
- **MMPEditor activity-column selector UI shape — single-select multi-list vs.
  multi-checkbox vs. column-list dropdown.** Original step 3 says "Select all
  two activities" — the exact UI invocation pattern depends on the
  `column_list` input renderer of `MMPEditor` per atlas
  `chem.analyze.mmp.editor` line 421. Common patterns: a multi-select listbox,
  a multi-checkbox list, or a column-picker dropdown with multi-select
  semantics. Migrated body resolves by deferring to spec-time introspection
  (Automator inspects `MMPEditor` dialog DOM and drives the actual control).
  Flag for QA pair review.
- **Diff types / scaling methods / fragment cutoff defaults.** Atlas
  `chem.analyze.mmp.editor` line 421 describes additional `MMPEditor` fields
  ("diff types, scaling methods, fragment cutoff") beyond molecules / activity
  columns. The original scenario does not specify any non-default values —
  the migrated body assumes defaults for these fields (the scenario's
  pristine-defaults walk). Automator at spec time picks defaults. Implicit
  assumption: the defaults are valid and lead to a successful MMP generation.
  Flag for QA pair review.

