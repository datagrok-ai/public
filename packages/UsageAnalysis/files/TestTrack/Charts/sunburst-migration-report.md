# Migration Report — sunburst.md

## Step mapping

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Step 1 — Open SPGI_v2.csv & demog.csv; Add Sunburst viewer to each. Expected: viewer opens without errors for both files | Setup steps 2–4 + `### Viewer creation smoke` steps 1–2 | preserved (split for clarity — opens moved to Setup, "without errors" assertion preserved as verification in Scenarios) |
| Step 2 — Click Gear icon. Expected: context panel with viewer properties opens | `### Property panel opens` steps 1–2 | preserved |
| Step 3.1 — Switch between SPGI_v2 / demog. Expected: viewer updates without errors | `### Table switching` steps 1–2 | preserved |
| Step 3.2 — Open Select Columns → Choose 2–4 columns → OK. Expected: hierarchy updates | `### Hierarchy configuration via Select Columns dialog` steps 1–3 | preserved |
| Step 3.2 — Reopen dialog → Search for column → Cancel. Expected: no changes applied | `### Hierarchy configuration via Select Columns dialog` steps 4–5 | preserved |
| Step 3.2 — Verify text fully visible / hidden with tooltip on hover | `### Hierarchy configuration via Select Columns dialog` step 6 | preserved as verification |
| Step 3.2 — Ensure structural columns render at correct size | `### Hierarchy configuration via Select Columns dialog` step 7 | preserved as verification |
| Step 3.3 — In demog: Select SEX → Enable Inherit from grid → Apply categorical coloring on grid SEX. Expected: Sunburst reflects grid colors | `### Inherit from grid (demog.csv)` steps 1–4 | preserved |
| Step 3.3 — Change coloring. Expected: viewer updates accordingly | `### Inherit from grid (demog.csv)` steps 5–6 | preserved |
| Step 3.4 — In SPGI_v2: select Core/R101 + Enable Include nulls. Expected: grey segments appear | `### Include nulls (SPGI_v2.csv)` steps 1–3 | preserved |
| Step 3.4 — Disable Include nulls. Expected: grey segments disappear | `### Include nulls (SPGI_v2.csv)` steps 4–5 | preserved |
| Step 4 — Double-click empty space OR context menu → Reset View. Expected: view resets | `### View reset` steps 1–4 | preserved (split for clarity — both reset paths exercised) |
| Step 5 — Click / Ctrl+Click / Ctrl+Shift+Click. Expected: grid rows update accordingly | `### Multi-selection behaviour` steps 1–6 | preserved |
| Step 6 — Open SPGI_v2; Select column with nulls (Sampling Time); Click null segment. Expected: behaves like any other category | `### Select / filter on empty category (SPGI_v2.csv)` steps 1–3 | preserved |
| Step 7 — Configure Sunburst with 3–4 columns on SPGI_v2 → Save the project; Close all → Reopen. Expected: viewer restores correctly | `### Project save and reopen` steps 1–5 | preserved |
| Step 7 — Save the layout and apply it. Expected: viewer retains saved layout | `### Layout save and apply` steps 1–4 | preserved (split into its own scenario block — distinct subsystem from project save) |
| Step 8 — Open the layout from issue #2979. Expected: viewer shows selected columns; Select Columns dialog in sync | `### Old layout compatibility (issue #2979)` steps 1–2 + Notes | preserved (with explicit fixture-provisioning gap flagged in Unresolved ambiguities) |
| Step 9 — In demog, configure Sunburst + apply internal AND panel filters. Expected: filters combine correctly; viewer shows intersection | `### Collaborative filtering (demog.csv)` steps 1–4 | preserved |

## Decisions

- **Why this `target_layer`:** chose `playwright` because the scenario
  requires save → close-all → reopen across a navigation boundary
  (`### Project save and reopen`) plus UI menu interaction for Save
  Layout / Apply Layout (`### Layout save and apply`), plus context-
  menu UI driving (`### View reset`). Chain `output_plan` for
  `sunburst.md` proposes `target_layer: playwright` for the same
  reasons (chain `charts.yaml` rev 1 line 209-217). Sibling
  `sunburst-spec.ts` already exists at the playwright layer per
  `existing-test-index.yaml` (line 32119) — confirms the layer choice.
- **Why this `pyramid_layer`:** kept `integration` per chain Rule 4
  (multi-subsystem co-existence). Original co-mixes viewer + Layouts
  + Projects + Filters + GitHub-#2979 old-layout-format compatibility.
  Did NOT shrink to `ui-smoke` — the smoke slot in the chain is owned
  by `radar.md` (3-step shortest qualifying scenario per Rule 1
  tie-break); see chain notes for `sunburst.md` `pyramid_layer`
  rationale.
- **Why this `coverage_type`:** `regression` — both cited bugs
  (`github-2954`, `github-3412`) carry `regression-risk` /
  `fixed-then-test-coverage-needed` status; this scenario is the
  primary mid-tier safety net for those Sunburst regressions.
- **Why this `strategy`:** `simple` per chain `output_plan` — single
  scenario, no chaining, no matrix, no fixture sharing
  (`fixtures_extracted: []`).
- **Sibling tests consulted:**
  - `public/packages/UsageAnalysis/files/TestTrack/Charts/sunburst-spec.ts`
    — already exists at playwright layer (currently `test.skip`); its
    feature coverage `[charts.sunburst]` confirmed.
  - `public/packages/UsageAnalysis/files/TestTrack/Charts/radar-spec.ts`
    — sibling viewer Charts spec at playwright layer
    (`existing-test-index.yaml` line 32099).
  - `public/packages/UsageAnalysis/files/TestTrack/Charts/tree-spec.ts`
    — sibling viewer Charts spec at playwright layer
    (`existing-test-index.yaml` line 32134).
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/complex-integration.md`
    — sibling integration-pyramid scenario; borrowed Setup / Notes
    structure (Authenticate → environment dependencies → cleanup), and
    explicit "no UI smoke consolidation here" framing.
- **Helpers reused:**
  - `findViewer` (registered in `helpers-registry.yaml :: grok_test_layer`).
  - `isViewerPresent` (registered).
  - `uploadProject` (registered) — candidate for
    `### Project save and reopen` Step 2.
  - No new helper invented. Downstream Automator may evaluate adding
    a Sunburst-specific `selectHierarchyColumns(viewer, columns)`
    helper if Select Columns dialog driving is repeated across specs;
    candidate noted here, NOT added to the registry by Migrator.
- **Bug library consulted:** yes — bugs `[github-2954, github-3412]`
  populate `related_bugs` in frontmatter. Both `affects` lists
  intersect `sub_features_covered` of this scenario
  (`charts.sunburst`, `charts.sunburst.inherit-from-grid`).
  Cross-cutting candidates from chain YAML
  `bug_focused_candidates: []` — no per-scenario chain candidate to
  cite. Chain `bug_match_attempts_skipped` records flagged both
  github-2954 (single-scenario, date-column repro path not exercised
  here — coverage gap) and github-3412 (single-scenario,
  Scatterplot side absent — cross-viewer pairing gap). Both gaps
  surfaced in chain YAML for Critic F downstream review; this
  Migrator does not invent the missing repro steps (date-column
  selection or Scatterplot pairing) because the chain analyzer
  explicitly classified them as gaps requiring scenario authoring,
  not migration.
- **Decision log queried:** yes — empty (`failed_attempts: []` at
  decision-log.yaml line 3286). No prior failed-attempt records for
  feature `charts`.
- **UI delegation:** none. This scenario does NOT delegate any UI
  flow to another spec — it owns its full UI surface
  (`ui_coverage_delegated_to: null` per chain). Sunburst-specialty
  flows (Select Columns dialog, Reset View context-menu, multi-select
  gestures, Save/Apply Layout) are not witnessed by `radar.md`'s
  smoke (which covers only Add Viewer + Gear/Property Panel) — so no
  SCOPE_REDUCTION cite required.

## Opt-outs (SCOPE_REDUCTION proposals)

(none)

_No step is opted out for effort. Every step from the original landed
in the migrated body or in a tracking section below._

## Deferred items (NOT opt-outs)

(none)

_All original steps are exercised in the migrated scenario. The
fixture provisioning for `### Old layout compatibility (issue #2979)`
is NOT a deferral — the step itself is preserved; only the fixture
sourcing is unresolved (see Unresolved ambiguities)._

## Edge cases

- **Date / unsupported column types in Select Columns dialog
  (`github-2954`):** the bug describes selecting a date column in the
  hierarchy and crashing the viewer with `NullError` in
  `tree-utils.ts`. Original Step 3.2 says "Choose 2–4 columns" without
  enumerating types — the date-column edge case is implicit in the
  generic instruction. Decision: **moved to atlas** as
  `edge_cases` candidate (atlas curator will action) — non-string
  column types (date, datetime, numerical) selected as Sunburst
  hierarchy columns must either be filtered by the selector or
  handled gracefully by the viewer. Chain analyzer already classified
  this as a coverage gap (`bug_match_attempts_skipped` for
  github-2954) — the migrated scenario does NOT invent the
  date-column step because the chain explicitly punted it to a future
  scenario authoring round.
- **Cross-viewer color-state pollution (`github-3412`):** the bug
  describes Sunburst + Scatterplot sharing a color column (e.g.
  Stereo Category) and Scatterplot's color editor breaking
  Sunburst's render. Original sunburst.md is single-viewer — no
  Scatterplot pairing. Decision: **moved to atlas** as `edge_cases`
  candidate (atlas curator will action). Chain
  `bug_match_attempts_skipped` for github-3412 already surfaced the
  cross-viewer pairing as a gap; Migrator does not invent the
  Scatterplot pairing step.
- **Sampling Time as a column with nulls
  (Step 6 example):** original cites Sampling Time as the canonical
  null-bearing column for the empty-category click test. Preserved
  as the explicit example in `### Select / filter on empty
  category (SPGI_v2.csv)` step 1 (parenthetical), with the
  understanding that any null-bearing column in SPGI_v2 satisfies
  the assertion if Sampling Time is renamed in a future fixture.
- **`onClick` action default (Select vs Filter):** atlas
  `charts.sunburst.on-click` notes the property defaults to
  `Select`. Original Step 9 ("apply both internal and panel
  filters") implies switching `onClick` to `Filter`, OR using a
  segment context-menu filter action. Decision: **preserved as
  step** in `### Collaborative filtering (demog.csv)` step 2
  (parenthetical clarification covers both paths); downstream
  Automator picks one driver path.
- **Layout save shortcut:** original Step 7 says "Save the layout"
  without specifying mechanism. Per persistent memory note
  ([feedback_no_ctrlS_for_layouts]), Ctrl+S saves a project, NOT a
  layout. Decision: **preserved as step** in `### Layout save and
  apply` step 1 with explicit "via Layouts menu / Save action — not
  Ctrl+S" guidance.
- **Auto-share dialog after Save Project:** original Step 7 does NOT
  mention the share dialog; sibling `complex-integration.md` Setup
  notes "Cancel auto-share". Preserved as a parenthetical instruction
  in `### Project save and reopen` step 2 ("Cancel any auto-share
  dialog that opens").

## Unresolved ambiguities

- **Step 8 — Old layout fixture sourcing (issue #2979):** the
  original step references "the layout from issue #2979" as if it
  were a known artifact, but the layout file is NOT stored in the
  repo or `System:DemoFiles`, and no `grok.dapi.layouts` provisioning
  is described. This was flagged by chain analyzer as
  `unresolved_ambiguities[1]` for `sunburst.md` — recommended
  resolution path: (a) commit a saved layout to repo as a fixture,
  (b) provision via `grok.dapi.layouts` pre-test, or (c) drop Step
  8 with documented reasoning. Migrator preserved Step 8 verbatim
  in `### Old layout compatibility (issue #2979)` and added an
  inline NOTE marking the provisioning gap. Downstream Automator
  must resolve before generating spec code; in the interim, the
  scenario block can be `test.skip`'d at the spec layer with the
  fixture as the cited blocker. Migrator did NOT invent the
  fixture path.
- **Old layout file location / ticket commentary:** the linked
  GitHub issue [#2979](https://github.com/datagrok-ai/public/issues/2979)
  is the bug that the layout was crafted to reproduce; the layout
  itself may be attached to the ticket OR derived from a
  pre-fix-version saved project. Migrator did not browse the ticket
  to confirm. Flag for QA pair review during fixture provisioning.
- **`Inherit from grid` toggle visibility on SPGI_v2 vs demog:**
  original Step 3.3 specifies the inherit-from-grid demo runs on
  demog (categorical coloring on SEX) and Step 3.4 (Include nulls)
  runs on SPGI_v2. Migrator preserved this split exactly. Downstream
  Automator should not consolidate the two onto a single table — the
  bug behaviour intended for surface (`charts.sunburst.inherit-from-
  grid` on demog SEX) requires the demog table specifically.
- **Helper candidate `selectHierarchyColumns(viewer, columns)`:**
  the Select Columns dialog appears in `### Hierarchy configuration`,
  `### Inherit from grid`, `### Include nulls`, and
  `### Project save and reopen` (4 of 11 scenario blocks). A helper
  would reduce duplication, but this is not a Migrator decision —
  flagged for `helpers-registry.yaml` curator. Migrator did not
  invent the helper name in the migrated body; only documented as
  a candidate here.
