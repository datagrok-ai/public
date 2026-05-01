# Migration Report — opening.md

## Step mapping

The original is 5 numbered steps + trailing JSON metadata
`{ "order": 3 }`. All 5 numbered steps are preserved; the parameterised
nature (over the 18 projects from `uploading.md`) is made explicit in
the migrated body via a project-name table.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Go to Browse > Dashboards" | Scenarios > "Open each ..." step 1 | preserved |
| 2. "Find projects from the previous step (Uploading)" | Scenarios > step 2 + project-name table (18 rows) | preserved-with-resolution. The phrase "projects from the previous step (Uploading)" is resolved per `scenario-chains/projects.yaml` rev 2 to the 18 explicit names produced by `uploading.md`. The parameter list is enumerated as an 18-row table in the migrated body so `chained_tests` strategy has a concrete input set. |
| 3. "Click each project. On Context Panel, check all attributes added (sharing, description, name, picture)" | Scenarios > steps 3-4 (single-click + verify-attributes) | preserved-with-axis-clarity-split. Original's single step 3 is split into migrated step 3 (single-click) and migrated step 4 (verify Context Panel attributes — Sharing, Description, Name, Picture, each as a sub-bullet). Mechanical D-STEP-02 auditability: each verification is explicit. |
| 4. "Double-click the project" | Scenarios > steps 5-6 (double-click + verify-on-open) | preserved-with-verification-split. Original's single step 4 (double-click) is followed by an implicit verification (the project must open successfully); migrated body splits this into step 5 (double-click) + step 6 (verify open: tables/viewers render, no console errors) for D-STEP-02 explicitness. |
| 5. "Close All" | Scenarios > step 7 | preserved as cleanup |
| Trailing JSON `{ "order": 3 }` | (dropped from body) | metadata-not-step (per orchestrator chain analysis convention; original `order: 3` captured in `scenario-chains/projects.yaml` rev 2 `order_from_files`) |

No original numbered step is silently dropped. The renumbering (5
original → 7 migrated) is axis-clarity expansion: the original's
multi-axis steps 3 and 4 are split into atomic verifications. No
semantic content is added or removed.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.opening.md.target_layer = playwright`. The scenario
  exercises Browse > Dashboards UI navigation, Context Panel
  rendering verification, and the project open lifecycle —
  persistence-bearing and UI-driven, ruling out `api-contract`.
  Existing `opening-spec.ts` is already Playwright; alignment keeps
  the section consistent.
- **Why this `priority`:** chose `regression` per the canonical
  A-STRUCT-MECH-06 enum (`smoke | regression | edge | perf`). The
  scenario is matrix-coverage (18 projects × 7-step verification
  flow) of the open + Context-Panel-render lifecycle — NOT a single
  golden-path smoke. Per-cycle Invariant 1 honored: Migrator
  SKILL.md:222 (`p0|p1|p2|p3`) NOT consulted.
- **Why this `strategy`:** `chained_tests` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.opening.md.strategy = chained_tests`. The 18 cases
  share one fixture (`uploading-18-projects`) and run identical
  workflows parameterised by project name — the canonical
  `chained_tests` shape (`describe.each` over a project list).
- **D-STRUCT-02 compliance:** all 18 paths preserved as a
  parameter table. No reduction applied at the .md level. The
  existing `opening-spec.ts` reduces to 1 project
  (`AutoTest-Opening-<timestamp>` per its line 18) — that reduction
  is an Automator-stage decision, NOT a Migrator-stage opt-out.
  The migrated `.md` remains the full 18-path source of truth for
  this scenario's intent.
- **Sibling tests consulted (READ-ONLY per Invariant 2 of the
  per-cycle override):**
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/opening-spec.ts` —
    confirmed Playwright; uses `softStep`, `evalJs`, `closeAll`,
    `loginToDatagrok` per section convention; covers Cases 1-5 with
    a single project. Read-only — not modified.
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading-spec.ts` —
    sibling at the playwright layer; confirms same `Date.now()`-suffix
    naming convention used across all 4 spec files in the section.
- **Helpers reused:** none invented. Existing `closeAll(page)` and
  `evalJs(page, ...)` patterns from sibling specs are the convention
  the Automator will follow at spec time. Migrator does not invoke
  helpers in the `.md` body — prose like "Close All" is the
  convention for migrated TestTrack scenarios.
- **Bug library consulted:** yes — `bug-library/projects.yaml`
  rev 2. No bug intersects this scenario's flow:
  - GROK-19750 (save-copy with Link), GROK-19212 (rename + datasync),
    GROK-19103 (join + save), GROK-19403 (share with un-shared dep),
    GROK-18345 (Spaces datasync share), GROK-19728 (view-and-use
    failure) — none exercise the plain open + Context-Panel-render
    flow. `related_bugs: []` in frontmatter is the correct call.
- **Decision log queried:** yes — `decision-log.yaml` rev 7 read.
  No prior `migration_decisions`, `layer_decisions`, or
  `manual_only` entry applies retroactively to `opening.md`
  content. The relevant per-cycle entries
  (`mig-2026-04-30-priority-enum-drift-discovery`,
  `mig-2026-04-30-uploading-migration`,
  `mig-2026-04-30-share-project-migration`,
  `atlas-2026-04-30-add-projects-shell-share-via-context-menu`,
  `atlas-2026-04-30-cascade-share-via-context-menu`) are honored
  by writing canonical priority + by NOT including
  `projects.shell.share-via-context-menu` in `sub_features_covered`
  (Invariant 3 of the override correctly does NOT apply — this
  scenario does not exercise the right-click Share dialog; its
  step-3 "Sharing" attribute on Context Panel is rendered metadata,
  not a share-flow operation).
- **Per-cycle override invariants (all 3) status:**
  - Invariant 1 (priority enum source-of-truth): honored —
    `priority: regression` is canonical per A-STRUCT-MECH-06.
  - Invariant 2 (existing -spec.ts / -api.ts READ-ONLY):
    SATISFIED — `opening-spec.ts` was read for sibling-test
    convention but NOT modified. Migrator's only writes are
    `opening.md` (overwrote source) and `opening-migration-report.md`
    (new file).
  - Invariant 3 (atlas-aware sub_features_covered for share):
    correctly does NOT apply — the scenario does not exercise the
    right-click Share dialog. Step 4's "Sharing" attribute is
    Context Panel rendering of pre-existing share metadata, not a
    new share operation.

## Opt-outs (SCOPE_REDUCTION proposals)

(none) — every numbered step of the original is preserved; the
parameterisation over 18 projects is preserved per D-STRUCT-02.

## Deferred items (NOT opt-outs)

- **`uploading-18-projects` fixture availability.** This scenario
  consumes 18 projects produced by `uploading.md`. Migration-stage
  cannot guarantee fixture presence at spec-run time; Automator
  stage owns the `beforeAll` fixture build (replay
  `uploading.md`'s save flow OR reuse from a prior chain run).
  Real prerequisite, not effort.
- **Cases 4-6 environment dependency on SPGIs Space (transitive
  via `uploading.md`).** Projects `Test_Case4_Sync`,
  `Test_Case4_NoSync`, `Test_Case5_Sync`, `Test_Case5_NoSync`,
  `Test_Case6_Sync`, `Test_Case6_NoSync` all require the SPGIs
  Space to exist on the test server (per `uploading.md` migration
  report's deferred items). On environments without SPGIs, those
  6 of 18 paths in `opening.md`'s parameter table will fail at
  step 2 (project not findable in Dashboards).
- **Project-name collision risk in CI.** All 18 project names
  are deterministic (`Test_Case<N>_Sync` / `Test_Case<N>_NoSync`).
  Concurrent CI runs of `uploading.md` would collide on these
  names, and `opening.md` would observe whichever instance
  remained. The existing `opening-spec.ts` works around this with
  a `Date.now()`-suffixed name (`AutoTest-Opening-<timestamp>`)
  but at the cost of using a fresh single project rather than the
  18 from `uploading.md`. Automator stage at spec time decides
  whether to (a) parameterise project names with timestamps, OR
  (b) serialize chain runs to avoid collision, OR (c) accept
  reduction to 1 fresh project.

## Edge cases

The original lists no explicit edge cases. Implicit edge cases
derivable from the steps:

- **Project with empty Description.** Step 4 (Verify Context
  Panel) lists "description" as an attribute. Projects from
  `uploading.md` save flow do not set a description; the
  Description attribute should render as empty/placeholder, not
  as a render error.
- **Project with empty Sharing list.** Same axis as above —
  `uploading.md` does not share its 18 projects, so the Sharing
  list will be empty (or show "Not shared"). The render must
  succeed; an empty list is not an error.
- **Project Picture rendering on first open.** Step 4 lists
  "picture" as a Context Panel attribute. `uploading.md`'s
  saved projects do not have explicit pictures set — a generated
  default thumbnail (or dashboard preview) is expected. The
  render must succeed even when no picture asset is bound.
- **Tables-and-viewers re-render on open (step 6).** Each of
  the 18 projects has 2+ tables (linked or joined), and Cases 7-9
  also have derivative tables (joins, pivots, aggregates).
  Step 6 verifies all such tables/viewers render on open without
  console errors. Implicit edge: a project with derivative
  tables that fails to re-render its derivatives is a defect.
  Verification is preserved as "no errors in the browser console
  (F12)".

(none additional)

## Unresolved ambiguities

- **"Projects from the previous step (Uploading)" — explicit list
  vs. dynamic.** The original phrasing implies the test should
  use whichever projects `uploading.md` happens to produce. The
  migrated body enumerates 18 explicit names per
  `scenario-chains/projects.yaml` rev 2's `produces` field. If
  `uploading.md` ever changes its naming convention or set of
  cases, this enumeration will need to be regenerated. Flag for
  retro: should the migrated body reference the chain's `produces`
  field by reference rather than by enumerated copy?
- **Order of iteration over 18 projects.** The original does not
  specify an order. Migrated body lists them in the natural
  (Case 1 Sync, Case 1 NoSync, Case 2 Sync, ..., Case 9 NoSync)
  order — chronological-by-case-number. If parallel iteration
  is preferred for spec-time speed, the order is irrelevant; if
  serial iteration is required, the chosen order is reasonable.
  Flag for Automator.
- **"Sharing" attribute rendering — what is checked exactly?**
  Step 4 lists "sharing" as one of the Context Panel attributes
  to verify. The original is silent on what specifically about
  the Sharing attribute is verified — presence of the section?
  Render without error? Empty list shown as "Not shared"? The
  migrated body says "list of users/groups with whom the project
  is shared (may be empty for solo-owner projects)". Automator
  decides spec-time assertion specificity.
- **Single-click vs. double-click on Dashboards card.** The
  original mixes single-click (step 3, to select for Context
  Panel inspection) and double-click (step 4, to open the
  project). The migrated body preserves both as distinct steps
  (3 and 5). Risk: if the Dashboards card UI changes such that
  single-click also opens the project (eliminating the
  inspect-only state), step 3 would no longer leave the project
  in a selected-but-not-opened state. Flag for Automator: verify
  the single-click-selects vs. double-click-opens distinction is
  still the UI behavior at spec time.
- **Order field in source.** Original `{ "order": 3 }` is the
  Test Track display ordering. With `upload-project.md` and
  `uploading.md` both at `order: 1` (collision per chain
  analysis rev 2 `unresolved_ambiguities`) and `share-project.md`
  at `order: 2`, this `opening.md` at `order: 3` is consistent
  with the dependency graph. Captured in chain artifact, not a
  Migrator concern.
