# Migration Report — radar.md

Migration date: 2026-05-07
Source: `public/packages/UsageAnalysis/files/TestTrack/Charts/radar.md`
Migrated by: grok-migrate-scenario Step 2 (section invocation, `--proceed-to-migration` for charts.yaml chain rev 1)

## Step mapping

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Step 1: "Open earthquakes.csv: On the Menu Ribbon, click the Add viewer icon, and select > Radar viewer." | Setup (clean session) + Scenarios > "Add Radar viewer to two datasets..." steps 1–2 | preserved (split for clarity) |
| Step 1 expected: "Radar viewer opens for earthquakes.csv." | Scenarios > step 3 | preserved as verification |
| Step 2: "Open demog.csv: On the Menu Ribbon, click the Add viewer icon, and select > Radar viewer." | Scenarios > steps 4–5 | preserved (split for clarity) |
| Step 2 expected: "Radar viewer opens for demog.csv." | Scenarios > step 6 | preserved as verification |
| Step 3 (action): "On the Radar viewer, click the Gear icon." | Scenarios > step 7 | preserved |
| Step 3 expected: "The Context Panel opens with viewer's properties." | Scenarios > step 8 | preserved as verification |
| Step 3 directive: "Check all properties." | Scenarios > step 13 | preserved (intentionally broad — see Unresolved ambiguities) |
| Step 3 main thing 1: "Switching tables" | Scenarios > step 9 | preserved as verification |
| Step 3 main thing 2: "Check-boxes in selection" | Scenarios > step 10 | preserved as verification |
| Step 3 main thing 3: "Increasing and decreasing the amount of chosen Values" | Scenarios > step 11 | preserved as verification |
| Step 3 main thing 4: "Style (color) changes" | Scenarios > step 12 | preserved as verification |
| Trailing metadata block (`order: 28`, `datasets: [System:DemoFiles/demog.csv]`) | Body trailer (preserved verbatim) + Notes bullet | preserved |

## Decisions

- **Why this `target_layer` (`playwright`):** chosen because a sibling
  `radar-spec.ts` already exists at the playwright layer
  (`public/packages/UsageAnalysis/files/TestTrack/Charts/radar-spec.ts`,
  `existing-test-index.yaml` line 32099, `layer: playwright`,
  `helpers_called: [spec-login]`). Chain rev 1 `output_plan` for
  radar.md proposed `uitests-package` but explicitly deferred the
  decision to Migrator: _"Decision deferred to Step 2."_ Picking
  `playwright` aligns the migrated scenario with the existing
  sibling spec (rather than orphaning it across two layers).
  **chain-vs-frontmatter divergence:** chain proposed
  `uitests-package`; migration honors `playwright` per chain's own
  deferral clause + sibling-spec evidence.
- **Why this `coverage_type` (`smoke`):** chain rev 1
  `dependency_graph.radar.md.pyramid_layer: ui-smoke` and
  `ui_coverage_plan.smoke_scenario: radar.md` — radar.md is the
  elected smoke owner for the entire Charts section (Add Viewer +
  Property-Panel-Gear surface). Smoke shape preferred over deep
  property-grid validation per chain rationale.
- **Why this `strategy` (simple, no chain):** chain rev 1
  `output_plan.radar.md.strategy: simple`,
  `dependency_graph.radar.md.depends_on: []`. Independent
  scenario; opens its own datasets in-session.
- **Sibling tests consulted:**
  - `public/packages/UsageAnalysis/files/TestTrack/Charts/radar-spec.ts` (playwright; uses `loginToDatagrok` + `softStep` + `specTestOptions` from `spec-login.ts`; mirrors original step structure 1→2→3 with `tv.addViewer('Radar')`).
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading.md` (migrated sibling — frontmatter pattern + Notes-section authoring style anchor).
- **Helpers reused (registered in `helpers-registry.yaml`):**
  - `softStep` (`spec-login.ts:18`) — Playwright soft-step wrapper.
  - `loginToDatagrok` (`spec-login.ts:26`) — login helper.
  - Implicit: `specTestOptions` from the same module.
  - These are referenced in the scenario `## Notes` for downstream
    Automator, not invented.
- **Helpers candidate (NOT invented; surfaced for registry curator):**
  - Candidate `helpers.charts.openWithRadarViewer(page, datasetPath)` — an
    open-dataset + addViewer('Radar') wrapper would consolidate the
    repeated Steps 1–2/4–5 pattern visible in `radar-spec.ts` lines
    27–43 (earthquakes) and 47–63 (demog). Propose for addition to
    `helpers-registry.yaml` externally; this Migrator does NOT add
    helpers per A1 boundary.
- **Bug library consulted:**
  - `yes — bugs [GROK-18085]` — `bug-library/charts.yaml` rev 1
    parsed; `GROK-18085` (`affects: [charts.radar,
    charts.echart-base.table]`) intersects this scenario's
    `sub_features_covered`. Listed in frontmatter `related_bugs`.
    GROK-19033, github-2954, github-3221, github-3245, github-3412
    do NOT intersect this scenario's coverage (they affect
    timelines / sunburst / tree).
  - **Bug coverage gap (cross-cutting awareness):** GROK-18085
    requires project save → reopen; the original radar.md has no
    such step. Chain rev 1 records this as
    `bug_match_attempts_skipped` for radar.md
    (skip_category: reproduction_unparseable, matched_steps: []).
    This Migrator surfaces the gap in the migrated scenario's
    `## Notes` and in this report; F-BUG-COVERAGE-01 at
    section-complete is the authoritative gate.
- **Decision log queried:**
  - `yes — empty for this feature` — `decision-log.yaml` was
    grepped for `charts` / `radar`. The only matches are
    coverage-map regenerations (line 6550, 6598, 6602) — no
    `failed_attempts` entries scoped to charts or radar. No
    layer-decision history to apply.

## Opt-outs (SCOPE_REDUCTION proposals)

(none)

_No step is opted out for effort. Every original step lands in
the migrated body or in a tracked report section. The Step 3
"Check all properties" instruction is preserved (intentionally
broad — see Unresolved ambiguities), not opted out._

## Deferred items (NOT opt-outs)

- **GROK-18085 reproduction (project save → reopen with
  table-rebound Radar)**: deferred until either (a) the original
  radar.md is amended to include a save/reopen step, or (b) a
  separate `radar-table-rebind-spec.ts` is authored under chain
  Pattern 2 (bug-focused slice). Prerequisite: amendment to
  scenario `.md` body OR chain-analyzer Mode 2 trigger threshold
  re-evaluation if scenario is amended. The Migrator does NOT
  unilaterally introduce save/reopen steps not in the original
  (would violate "no silent additions" — symmetric to the no-
  silent-drops rule). Real prerequisite cited: chain
  `bug_match_attempts_skipped.GROK-18085.skip_category:
  reproduction_unparseable` (radar.md has no save/reopen body).

## Edge cases

- **Edge case (atlas, source_bug GROK-18085):** "Switching a Radar
  viewer's bound table (tableName) mid-session before saving the
  project breaks deserialization on reopen." (atlas
  `edge_cases[0]`, `sub_features: [charts.radar,
  charts.echart-base.table]`).
  - **Decision:** the **table-switch (mid-session)** half of this
    edge case IS exercised by Scenarios > step 9 (switching the
    bound `table` property between `earthquakes` and `demog`). The
    **save → reopen** half is NOT exercised here (see Deferred
    items). Edge case partially preserved + partially deferred;
    flagged for atlas curator awareness that no single Charts
    scenario currently witnesses both halves end-to-end.

- **Implicit edge case (original Step 3 "showCurrentRow / mouseOver
  / multi-row selection"):** the original scenario calls out
  "Check-boxes in selection" as a Main Thing reflected on the
  viewer. The atlas (`charts.radar.show-current-row`) lists 3
  related properties (showCurrentRow, showMouseOverRow,
  showMouseOverRowGroup).
  - **Decision:** preserved as Scenarios > step 10 (multi-row
    selection → verify selected-row lines reflected on viewer).
    Atlas-listed mouseover variants are NOT enumerated as separate
    steps — the original scenario does not specify mouseover; the
    show-current-row family is treated as a single representative
    behavior per the original's level of specificity. No silent
    drop; documented here.

## Unresolved ambiguities

- **Step 3 "Check all properties" (chain rev 1 unresolved
  ambiguity #1):** the original scenario instructs "Check all
  properties" without enumerating which of the 9 atlas-listed
  `charts.radar.*` properties (`charts.radar.title`,
  `charts.radar.min-percentile`, `charts.radar.max-percentile`,
  `charts.radar.show-current-row`, `charts.radar.show-tooltip`,
  `charts.radar.color-column`, `charts.radar.color-palette`,
  `charts.radar.show-min-max`, `charts.radar.legend-visibility`)
  are required vs. optional. The bullet list ("Switching tables /
  Check-boxes in selection / Values count / Style colors") gives 4
  representative items. **Migrator decision:** preserved verbatim
  as Scenarios > step 13 (broad sweep) plus the 4 explicit Main
  Things as Scenarios > steps 9–12. Did NOT narrow the original to
  just the 4 cited (would silently drop the broader instruction);
  did NOT widen to all 9 atlas properties (would invent specificity
  the original does not have). Chain analyzer recommended
  "narrowing to cited items + atlas cross-reference"; Migrator
  declined to narrow per migration-prompt.md "preserve as-is" rule
  for ambiguity. Flag for QA pair review.

- **Original lacks save/reopen step but bug GROK-18085 cites it
  as the reproduction path:** the chain-analyzer flagged this in
  `bug_match_attempts_skipped.GROK-18085`. Migrator's choice was
  whether to (a) add save/reopen steps to the migrated body
  (silent expansion), (b) drop the related_bugs citation, or (c)
  cite the bug + flag the gap. Chose (c): `related_bugs:
  [GROK-18085]` in frontmatter (per user instruction +
  bug-library affects intersection) + Notes-section gap
  acknowledgement + Deferred-items entry for a future
  bug-focused slice. The merit of this choice: keeps the trail
  auditable for Critic F's section-level bug-coverage assessment
  without inventing scenario steps not in the original.

- **Step 3 "Increasing and decreasing the amount of chosen
  Values":** the original phrasing is ambiguous as to whether
  "Values" means the radar's `colorColumnName` (single column) or
  the multi-column "Values" bag that defines axes. Atlas
  `charts.radar` does not have an explicit "Values" property —
  axes are derived from non-color numeric columns; coloring uses
  `colorColumnName`. **Migrator decision:** preserved as
  Scenarios > step 11 referencing "the columns chosen as Values
  (the color column / values bag)" — both interpretations
  acceptable per the migrated step. Flag for downstream Automator
  to pick a single property path (likely
  `viewer.setOptions({colorColumnName: ...})` and varying numeric
  column visibility on the bound table).
