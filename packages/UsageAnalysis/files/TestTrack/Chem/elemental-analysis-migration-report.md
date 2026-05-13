# Migration Report — elemental-analysis.md

First-cycle migration of `elemental-analysis.md` (TestTrack/Chem section) per chain YAML
`scenario-chains/chem.yaml` rev 2 (`output_plan.elemental-analysis.md`:
`target_layer: playwright`, strategy `simple`, classification `simple`,
`pyramid_layer: integration`, `depends_on: []`). Multi-format happy-path walk over
`chem.analyze.elemental`; realizes atlas top-menu surface
(`chem.analyze.elemental.top-menu`) + transform (`chem.analyze.elemental.run`) across the
three notation variants the original body promises.

## Step mapping

The original is a 9-line, 4-numbered-step scenario with a step-numbering glitch (`1.` / `2.` /
`2.` / `3.`) plus a JSON footer with `order: 4` and a single dataset entry (`smiles.csv`).
Body defects per chain rev 2 directive (Olena 2026-05-11, chain footer note (i)
"Source-text defects in scenario bodies — Migrator fixes during migration. Affected: …
elemental-analysis.md (step 1/2/2/3 renumber + dataset list extension from smiles.csv to
smi+molV2000+molV3000)…") fixed silently.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Step 1 — "Open linked datasets (check on smiles, molV2000 and molV3000 formats)" | Setup step 1 (provisioning — datasets D1, D2, D3 enumerated) + Scenarios > "Elemental Analysis per format" step 1 (dataset open as the per-variant loop entry) | preserved (split for clarity — provisioning declaration in Setup; the open action in scenario body) |
| Step 2 (first) — "On the menu ribbon, open **Chem > Analyze > Elemental Analysis**. A dialog opens." | Scenarios > "Elemental Analysis per format" step 2 (instantiated per variant — top-menu walk + dialog open) | preserved (renumbered cleanly; "menu ribbon" generalized to "top menu" matching the section-wide playwright convention used in sibling `calculate.md` migration) |
| Step 2 (second, collision) — "Turn all checkboxes on." | Scenarios > "Elemental Analysis per format" step 3 (instantiated per variant — toggle every per-element checkbox on) | preserved (renumbered from collided `2.` → `3.`; phrasing made explicit "every per-element checkbox" for the calculator's checkbox shape per atlas `chem.analyze.elemental` description "per-element atom-count columns") |
| Step 3 — "Click OK to run a function with default parameters." | Scenarios > "Elemental Analysis per format" step 4 (instantiated per variant — OK click) | preserved (renumbered `3.` → `4.`; "default parameters" preserved as semantic — the all-checkboxes-on selection overrides defaults for elements; remaining inputs default per the original) |
| _(implicit)_ "Make sure the result is correct" — original body asserts the calculator runs and adds columns by the structural promise of Step 4 + atlas description "per-element atom-count columns are added to the table" | Scenarios > "Elemental Analysis per format" step 5 (instantiated per variant — verify per-element atom-count columns appended; row count unchanged; no console errors) | preserved as verification (silent fix — added explicit verification step concretizing the original's structural assertion per the section's playwright-test convention; sibling `calculate.md` step 5 + `info-panels.md` Phase B step 6 use the same per-cell verification pattern) |
| JSON footer `{ "order": 4, "datasets": ["System:DemoFiles/chem/smiles.csv"] }` | (dropped from body; extended in chain + migrated Setup) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` + migrated Setup step 1 dataset enumeration — extended from 1 entry to 3 entries per chain rev 2 directive) |

No original step is silently dropped. Step-numbering collision (original had `1.` / `2.` / `2.`
/ `3.`) renumbered cleanly into the migrated body's `Setup` + `Scenarios > "Elemental Analysis
per format"` linear sequence.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/chem.yaml` `output_plan.elemental-analysis.md.target_layer = playwright`.
  The scenario requires DOM-level UI driving: top-menu walk
  (**Chem > Analyze > Elemental Analysis...**), per-element checkbox toggling in the dialog,
  OK click on the dialog, and grid-column-appended verification on the active table view — all
  UI behaviors that Playwright drives natively and the apitest layer cannot exercise. Chain YAML
  reason field: "Standalone Chem | Analyze | Elemental Analysis run with all element checkboxes
  toggled. Sibling elemental-analysis-spec.ts already exists at playwright per
  existing-test-index. Single scenario; simple strategy."
- **Why this `coverage_type`:** chose `regression` per chain rev 2 directive
  (`pyramid_layer: integration` — multi-format walk; not `smoke` since the section's smoke is
  `Advanced/scaffold-tree-functions.md` per chain `ui_coverage_plan.smoke_scenario` — and chain
  rev 1 notes a strict-Rule-1 deviation: strict step-count tie-break would pick this scenario
  (3 main steps = shortest) but the analyzer prefers the Scaffold Tree viewer smoke because
  Elemental Analysis is a Chem|Analyze calculator-style flow, not a viewer create+configure
  smoke; not `edge` / `perf` since the scenario walks the happy path only — no specific
  failure-mode invariant being asserted, no timing or volume threshold). The optional Radar
  viewer / per-row Radar grid cell warning behavior (atlas `chem.analyze.elemental.top-menu`:
  "Warns if Charts / PowerGrid are missing") is out of scope for this scenario.
  `regression` is the natural fit for a multi-format walk that guards against regression-of-
  the-set across the three notation formats.
- **Why this `strategy`:** `simple` per chain YAML
  `output_plan.elemental-analysis.md.strategy = simple`. Single scenario, no cross-file
  fixture; chain analyzer classified `simple` (single calculator across 3 format variants).
  Pattern 1 (decomposition) is NOT triggered (chain classifier returned `simple`, not
  `complex-standalone`; no child scenarios with `produced_from: decomposed`). Pattern 2
  (bug-focused slice) does NOT apply (this scenario is NOT listed in chain
  `bug_focused_candidates[].spans`; no GROK ticket invariant overlay). Pattern 3 (UI-smoke
  consolidation) is NOT triggered (no `ui_consolidation_proposals[]` from Gate F — chain
  `consolidation_history: []` and Edit 8 emitted empty per references/chem.md absence).
  Standard per-scenario migration per "Migrated scenario shape".
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/elemental-analysis-spec.ts` —
    existing playwright-layer test ("Chem: Elemental Analysis", category Chem,
    `features_covered: [chem.elemental-analysis]`, helpers `spec-login`, patterns
    `uses-fixture-runners` / `uses-grok.dapi` / `uses-grok.shell` / `uses-page.evaluate` /
    `uses-page.locator` / `uses-playwright-test`). Per `existing-test-index.yaml` line 32462.
    Existing spec likely covers a single format; Automator will extend to the 3-format
    enumerated walk. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/calculate.md`
    + `calculate-migration-report.md` (same-cycle predecessor) — shape anchor for migrated
    `.md` + report structure (frontmatter, Setup / Scenarios / Notes order, per-variant data-
    driven matrix walk; SR-01 A-STRUCT-02 carryforward template). Same Calculate-menu-style
    calculator flow (run dialog → assert columns added). Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels.md` +
    `info-panels-migration-report.md` (same-cycle predecessor) — section-wide multi-format
    walk house-style anchor (Phase A / Phase B data-driven format matrix). Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/activity-cliffs.md` — adjacent
    same-cycle multi-format scenario; the dataset enumeration pattern in its JSON footer
    (line 16-22: 5 datasets including `smiles.csv`, `mol1K.sdf`, `ApprovedDrugs2015.sdf`) is
    the precedent for the multi-format dataset extension applied here. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/scaffold-tree-functions.md`
    + migration report — chain ui-smoke scenario; cited to triangulate the
    smoke-vs-regression coverage_type decision (chain rev 1 notes the Rule 1 strict
    step-count deviation favoring scaffold-tree-functions.md as smoke).
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok`
      (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`)
      — section-standard fixture login; Automator will use in spec `beforeEach`.
    - `softStep`
      (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`)
      — per-step soft-fail wrapper; Automator wraps each per-variant verification in
      `softStep` for granular failure reporting (3 variants × ~5 sub-actions per variant).
    - `closeAllViews`
      (`public/packages/UITests/playwright/scripts/helpers.ts:129`)
      — cleanup between variants of the per-format walk (per "Close the active view before
      moving to the next dataset" instruction in the migrated body).
  - **Candidate helpers (NOT yet in registry — flagged for addition via helpers-registry
    curator; per migration-prompt Helpers discipline §):**
    - `helpers.playwright.chem.openChemTopMenuItem(page, menuPath)` — walks
      **Chem > Analyze > ${menuItem}** (or any top-menu path). Surfaced by Scenarios step 2
      per variant. Convenience abstraction shared with `calculate.md`'s candidate helper
      surface; spec author can inline the top-menu walk until the helper lands.
    - `helpers.playwright.chem.toggleAllCheckboxes(page, dialogTitle)` — accepts an open
      Chem analysis dialog and toggles every checkbox on. Surfaced by Scenarios step 3.
    - `helpers.playwright.chem.assertColumnsAppended(page, prevColumnCount)` — asserts the
      grid column count grew without changing row count (cardinality-only assertion suitable
      for an N-checkbox-driven calculator where the exact column names depend on the dataset
      content's element distribution). Surfaced by Scenarios step 5.
- **Bug library consulted:** yes — `bug-library/chem.yaml` revision 1
  (12 curated_bugs, generated 2026-05-05). Grepped `curated_bugs[]` for "elemental" — zero
  matches. None of the 12 curated bugs have `affects` intersecting `chem.analyze.elemental` /
  `chem.analyze.elemental.top-menu` / `chem.analyze.elemental.run`. `related_bugs: []` in
  frontmatter.
- **Decision log queried:** yes — `decision-log.yaml` grepped for `feature: chem` —
  **zero matches** (consistent with the calculate.md / info-panels.md / scaffold-tree /
  r-group-analysis migration reports in the same cycle). This is the FIRST CYCLE for the chem
  section. No prior `migration_decisions`, `layer_decisions`, `manual_only`, or
  `failed_attempts` entries exist for chem. No "approaches off the table" constraints apply.
- **Cross-cutting bug citations (chain YAML `bug_focused_candidates[]`).**
  Per `migration-prompt.md` "Cross-cutting bug citations from chain YAML" §: this scenario is
  NOT a span in any chain `bug_focused_candidates[]` entry (verified per chain rev 2:
  `chem-grok-17621-spec.ts` spans `calculate.md:Step 1`; `chem-github-2942-spec.ts` spans
  `calculate.md:Step 1` + `filter-panel.md:Step 1`; remaining bug-focused candidates target
  scaffold-tree / sketcher / mmp / chemprop — none reference `elemental-analysis.md`).
  No cross-cutting citation required. Per migration-prompt this is RECOMMENDED, not
  mandatory; F-BUG-COVERAGE-01 at section-complete is authoritative.
- **UI delegation status.** Per chain YAML `ui_coverage_plan.smoke_scenario:
  Advanced/scaffold-tree-functions.md` (this scenario is NOT the section smoke) and the chain
  entry `ui_coverage_plan.delegations[].scenario: elemental-analysis.md` /
  `delegated_to: null`, this scenario owns its UI coverage directly (no upstream delegation).
  Chain `ui_coverage_responsibility: [chem-elemental-analysis-top-menu,
  chem-elemental-analysis-checkboxes]`. No SCOPE_REDUCTION proposal substitutes JS API for
  either of these UI flows; therefore the UI delegation discipline of migration-prompt §
  "UI delegation in SCOPE_REDUCTION proposals" does not apply (no SR proposal to delegate).
  Reaffirmed as a hard constraint in the migrated Notes ("No JS API substitution").
- **Scenario constraint extraction (per migration-prompt § Scenario constraint extraction).**
  - **(a) FORBIDDEN substitutions:** None explicit in the original body (no "must not",
    "forbidden", "not allowed" phrasing). HOWEVER `pyramid_layer: integration` per chain
    implicitly forbids JS API substitution for the 2 `ui_coverage_responsibility` flows —
    surfaced explicitly in the migrated Notes (`No JS API substitution`) as a structured
    behavioral guard. Additionally invoking `chem.analyze.elemental.run` (atlas
    `runElementalAnalysis(table, molecules)`, role `transform`) directly via
    `grok.functions.eval` or `DG.Func.call` would defeat the dialog-driven workflow
    (per-element checkbox selection) that the scenario asserts.
  - **(b) REQUIRED actions:** Both entries of `ui_coverage_responsibility`
    (`chem-elemental-analysis-top-menu`, `chem-elemental-analysis-checkboxes`) are exercised
    via UI driving in the migrated body (top-menu walk per variant + all-elements-on
    checkbox toggle). No flow is deferred.
  - **(c) Missing-selector escalation:** Migrator does NOT have explicit selector
    definitions for the Elemental Analysis dialog's per-element checkbox set in
    `grok-browser/references/projects.md`. Spec-time selector discovery is the Automator's
    responsibility (per existing sibling `elemental-analysis-spec.ts` patterns +
    `page.evaluate` introspection); a reference-file approval-required proposal to add
    Elemental-Analysis-dialog-specific selectors to
    `grok-browser/references/widgets/elemental-analysis-dialog.md` is a Phase 2
    deliverable — out of this migration's scope.
  - **(d) Reference templates:** Code-style anchors for the Automator (cited verbatim, not
    modified here):
    `public/packages/UsageAnalysis/files/TestTrack/Chem/elemental-analysis-spec.ts`
    (current existing test — section-mate spec; Automator extends per the enumerated
    3-format walk) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/calculate-spec.ts` (section-mate
    Calculate-menu-style calculator walk; `spec-login` + `page.evaluate` + `page.locator`
    patterns — directly parallel given the calculator-style dialog shape) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts` (section-mate
    multi-format playwright walk).
  - **(e) Scenario authority clause:** No prompt-vs-scenario conflict surfaced. The
    original scenario is a 4-step linear walk anchored on Elemental Analysis with all
    element checkboxes on; the migrated body preserves all 4 anchor steps (renumbered) and
    adds a concrete verification step per the prompt directive ("Add explicit verification
    step 'verify per-element atom-count columns added to table'"). No
    `prompt_scenario_conflict:` field to flag for orchestrator.
- **Source-text fixes (silently applied per chain rev 2 directive, Olena 2026-05-11 chain
  footer note (i)).**
  - **Step-numbering glitch `1.` / `2.` / `2.` / `3.`.** Original has a collided `2.`
    (between "open Chem > Analyze > Elemental Analysis" and "Turn all checkboxes on") and
    then `3.` (the OK click). Renumbered cleanly to `1.` / `2.` / `3.` / `4.` in the
    migrated body's `Setup` + `Scenarios > "Elemental Analysis per format"` structure.
  - **Dataset list extension.** Original JSON footer declared a single dataset
    (`System:DemoFiles/chem/smiles.csv`) but the body promises "check on smiles, molV2000
    and molV3000 formats". Per chain rev 2 directive (Olena 2026-05-11 chain footer note
    (i)) the JSON dataset list is extended during migration to include
    `System:AppData/Chem/mol1K.sdf` (molV2000) and
    `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` (molV3000) — same paths used by
    sibling `activity-cliffs.md` (lines 19-20 of JSON footer) and `calculate.md` (Variants
    B / C). Provisioned in Setup step 1 as D1 / D2 / D3.
  - **Explicit verification step added.** Original Step 4 ("Click OK to run a function
    with default parameters") implicitly relied on the structural promise of the calculator
    + atlas `chem.analyze.elemental` description ("per-element atom-count columns are added
    to the table") for verification. Per the prompt directive ("Add explicit verification
    step 'verify per-element atom-count columns added to table' (preserves original step 4
    intent more concretely)") and the section-wide playwright-test convention used in
    sibling `calculate.md` step 5 + `info-panels.md` Phase B step 6, the migrated body's
    Scenarios step 5 surfaces the explicit verification ("per-element atom-count columns
    are appended to the active table view's grid (one column per element that was
    checked, with integer atom counts per row). The grid row count is unchanged. No
    console errors fire during dialog open, OK click, or column-append."). Same fix shape
    that the activity-cliffs.md and calculate.md migrations applied to their respective
    column-append verifications.
- **Parameterized data_driven strategy.** The migrated body presents the multi-format walk
  as a 3-variant cell loop (D1 / D2 / D3) per migration-prompt § "Matrix scenarios from the
  chain". The downstream spec author (Automator) is expected to realize this as a Playwright
  `test.describe` block iterating datasets (or a `test.each`-style fixture) with per-variant
  soft-step verification — directly mirroring the calculate.md / info-panels.md
  data-driven shape.
- **No invented sub_features / helpers.** Every id in `sub_features_covered`
  (`chem.analyze.elemental`, `chem.analyze.elemental.top-menu`, `chem.analyze.elemental.run`)
  is present in `feature-atlas/chem.yaml` rev 2 (verified per atlas lines 465-486). Every
  reused helper is in `helpers-registry.yaml` (verified per lines 3543 / 3597 / 3603);
  candidate helpers are flagged as candidates, NOT referenced by name in the migrated body.

## Opt-outs (SCOPE_REDUCTION proposals)

### SR-01: A-STRUCT-02 carryforward (chain-level edge/perf coverage)

**Cited technical dependency:** A-STRUCT-02 ("at least one scenario is edge case or negative
path — i.e. its frontmatter `coverage_type` is `edge` or `perf`") is a **section-level
structural invariant**. This scenario (elemental-analysis.md) is a multi-format happy-path
calculator walk — `coverage_type: regression` is the natural fit per chain rev 2
(`pyramid_layer: integration`, classification `simple`, multi-format walk). Forcing
`coverage_type: edge` would mis-classify (the scenario asserts happy-path column-append
across 3 notation formats, not a specific failure-mode invariant such as missing-Charts /
missing-PowerGrid warning, malformed-molecule parse error, or empty-table behavior).

**Resolution path:** A-STRUCT-02 satisfaction is chain-wide via 10 `bug_focused_candidates[]`
in `scenario-chains/chem.yaml` rev 2 (each chem-bug-focused spec naturally carries
`coverage_type: edge` once authored — none of the 10 directly target Elemental Analysis but
the section-level invariant is satisfied chain-wide). Additionally `r-group-analysis.md`
(same-cycle, same section, `coverage_type: edge`) already organically satisfies the
section-level invariant. Decision-log line 8528-8530 records the radar.md 2026-05-07
carryforward precedent. Same SR-01 carried by sibling-cycle reports `calculate-migration-
report.md`, `scaffold-tree-functions-migration-report.md`, `scaffold-tree-migration-
report.md`, `info-panels-migration-report.md`, `Advanced/structure-filter-migration-
report.md`, `Advanced/similarity-search-migration-report.md`.

**Scope of waiver:** A-STRUCT-02 verdict on `elemental-analysis.md` is deferred to
chain-level evaluation. Critic A should return SR (carryforward), not FAIL.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through A-STRUCT-06, A-COVERAGE-*,
A-MERIT-*) are expected to PASS on `elemental-analysis.md` without any other SR proposals.
No step is opted out for effort. The single entry above cites a real section-level
structural property whose satisfaction path is owned by the chain, not by this per-scenario
migration (D-MERIT-01 compliant)._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist yet. The migrated body's
3-dataset Elemental Analysis walk is realizable against the current playwright +
helpers-registry + atlas state. The candidate helpers surfaced in Decisions are convenience
abstractions that the Automator may inline pending registration — they do NOT block spec
realization. The selector-reference write proposal (Decisions § Scenario constraint
extraction (c)) is a downstream Phase 2 enhancement, not a blocker on this scenario's
automation (existing `elemental-analysis-spec.ts` already drives the Elemental Analysis
dialog per the sibling-test pattern; Automator can extend that pattern across the full
3-format walk without a reference update). The optional Radar viewer / per-row Radar grid
cell add-on surfaces (atlas `chem.analyze.elemental` Charts + PowerGrid follow-ons) are
out of scope of the original body, not deferred — they belong to a separate scenario
(e.g. an atlas-curator-flagged future `elemental-analysis-radar.md` if those surfaces are
deemed test-worthy)._

## Edge cases

The original lists no explicit "edge case" keyword. Implicit edge cases derivable from the
scenario content + atlas:

- **Cross-format calculator invocation.** Elemental Analysis is invoked across 3 distinct
  notation formats (smiles, molV2000, molV3000). Implicit invariant: the calculator must
  support every notation format, OR surface a clear user-visible error (the per-variant
  assertion in the migrated body — "either per-element columns are appended OR a
  user-visible error balloon is shown — never silent"). PRESERVED as scenario step
  (Scenarios cross-variant invariant clause).
- **All-elements-on selection.** The original explicitly directs "Turn all checkboxes on"
  — exercises the maximum cardinality of per-element columns the calculator can append.
  PRESERVED as scenario step (Scenarios step 3 — "Toggle every per-element checkbox in
  the dialog **on** (select all elements).").
- **Element-set variance across formats.** The exact set of per-element atom-count columns
  appended depends on which elements appear in the dataset's molecules. For
  `smiles.csv` / `mol1K.sdf` / `ApprovedDrugs2015.sdf` the element set is dataset-specific
  (e.g. C/H/N/O/S/P/halogens for organic / drug-like sets). The migrated body's step 5
  asserts "per-element atom-count columns appended" without enumerating exact column
  names — the assertion is on column-count delta, not on exact names. PRESERVED as
  scenario step (Scenarios step 5).
- **Missing Charts / PowerGrid packages (Radar follow-on warning).** Atlas
  `chem.analyze.elemental.top-menu` documents "Warns if Charts / PowerGrid are missing".
  OUT OF SCOPE for this scenario (no Radar viewer assertion, no per-row Radar grid cell
  assertion). Flagged as atlas-curator candidate edge case ("Elemental Analysis Radar
  add-on warning when Charts/PowerGrid absent: scenario coverage opportunity"); not in
  this migration's scope.
- **Empty-table / no-molecule-column dataset.** Not exercised (all 3 variants have
  populated molecule columns). Flagged as atlas-curator candidate edge case ("Elemental
  Analysis dialog behavior on dataset with no molecule column: expected behavior is the
  top-menu item is disabled OR the dialog surfaces a user-visible message"); not in this
  migration's scope.
- **Malformed-molecule rows in dataset.** Not exercised (all 3 variants have well-formed
  molecules per the platform's bundled demo data). Flagged as atlas-curator candidate
  edge case; not in this migration's scope.
- **Console errors throughout.** Implicit across all 3 variants: each variant's dataset
  open + dialog open + OK click + column-append must complete without console errors.
  PRESERVED as a per-variant verification (Scenarios step 5 — "No console errors fire
  during dialog open, OK click, or column-append"). Sister with `info-panels.md` /
  `calculate.md` / `scaffold-tree-functions.md` patterns.

No edge case is moved to atlas, manual_only, deferred, or a separate scenario silently.
The optional Radar viewer / empty-table / malformed-molecule edge cases are flagged for
atlas curator (no atlas write here).

## Unresolved ambiguities

- **"Turn all checkboxes on" — exact checkbox set on the Elemental Analysis dialog.** The
  Elemental Analysis dialog's per-element checkbox set is dictated by the calculator's
  implementation (`elementalAnalysis` at `public/packages/Chem/src/package.ts#L961`). The
  migrated body's step 3 ("Toggle every per-element checkbox in the dialog **on** (select
  all elements).") is non-deterministic on the exact list of elements presented — the
  Automator at spec time iterates whatever checkbox set the dialog renders. If the dialog
  exposes only "Select all" / "Deselect all" master controls, Automator clicks "Select
  all". Flag for QA pair review.
- **Expected column name pattern after Elemental Analysis.** The columns appended after
  Elemental Analysis are typically named after the element symbol (e.g. "C", "H", "N",
  "O") or atom-count-prefixed (e.g. "atoms_C"). The migrated body's step 5 asserts column
  count delta + integer atom-count cell values rather than exact column names. Automator
  at spec time may use a column-count-delta assertion (`afterCols.length >
  beforeCols.length`) or a predicate assertion on column-name patterns. Flag for QA pair
  review (cf. same flag in `calculate-migration-report.md` Unresolved ambiguities §
  "Expected column name(s) per Calculate menu item").
- **Per-variant column count delta.** The atom-count column delta varies by dataset (each
  dataset has a different element-set distribution). The migrated body's step 5 asserts
  the delta is positive (`> 0`) but does not specify an exact count per variant.
  Automator at spec time computes the delta dynamically; the assertion is monotonic, not
  exact. Flag for spec design at Automator stage.
- **Cell ordering inside the matrix.** The migrated body presents D1–D3 in dataset
  enumeration order (smiles / molV2000 / molV3000) but does NOT enforce a strict
  sequential ordering across variants. The Automator at spec time can parallelize
  variants or run them sequentially; the scenario's invariants are per-variant (no
  inter-variant dependencies). Flag for spec design at Automator stage.
- **Original "Make sure the result is correct" implicit verification.** The original Step
  4 ("Click OK to run a function with default parameters") provides no explicit "expected
  result" clause; the assertion is implicit ("the result is correct"). The migrated body's
  step 5 concretizes this as "per-element atom-count columns are appended … with integer
  atom counts per row" per the atlas description. If the original intent was a stricter
  semantic check (e.g. "atom counts match the molecule's actual atomic composition"), the
  migrated assertion is weaker (cardinality-only, not value-correctness). Flag for QA
  pair review — value-correctness verification (e.g. for a known SMILES like benzene
  `c1ccccc1`, the "C" column should equal 6 and "H" column should equal 6) is a stronger
  future-work edge case for an atlas-curator-flagged dedicated scenario.
