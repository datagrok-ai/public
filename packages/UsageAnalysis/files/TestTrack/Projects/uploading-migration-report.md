# Migration Report — uploading.md

## Step mapping

The original is a 9-case matrix with each case running twice (Data Sync
ON, Data Sync OFF) — 18 logical paths total. Step numbering in the
original is per-case (each Case starts at step 1). The migrated body
preserves all 9 cases under a single `## Scenarios` H2 with each Case
as an `### Case N` H3, and preserves both Sync ON and Sync OFF variants
within each Case body.

| Original location                    | Migrated location                   | Decision                                                                                                                                                                                                          |
|--------------------------------------|-------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `#### Matrix: Table source combinations` (header + 9-row table) | `## Scenarios > ### Matrix: Table source combinations` (table) | preserved (matrix table moved under the `## Scenarios` H2; row content unchanged)                                                                                                                                 |
| Case 1, original steps 1–6            | Case 1, migrated steps 1–6          | preserved (Browse > Files navigation, dataset double-click, Link Tables, verify linking)                                                                                                                          |
| Case 1, original step 7 ("Save with Data Sync") | Case 1, migrated step 7 ("Save with Data Sync ON") | preserved (renamed to "Save with Data Sync ON" for axis-clarity vs. step 10's OFF)                                                                                                                                |
| Case 1, original steps 8–9            | Case 1, migrated steps 8–9          | preserved (close, reopen, verify on reopen)                                                                                                                                                                       |
| Case 1, original step 10 ("Save without Data Sync") | Case 1, migrated step 10 ("Save with Data Sync OFF") | preserved (renamed for axis-clarity; "Repeat steps 1–6" idiom retained)                                                                                                                                           |
| Case 1, original steps 11–12          | Case 1, migrated steps 11–12        | preserved                                                                                                                                                                                                         |
| Case 2, original steps 1–11           | Case 2, migrated steps 1–12         | preserved (the original collapses verify+steps into 11 steps; migrated separates Sync-ON-reopen-verify from Sync-OFF-reopen-verify into distinct steps for D-STEP-02 explicitness — 12 numbered steps in migrated) |
| Case 3, original steps 1–10           | Case 3, migrated steps 1–12         | preserved (same Sync-ON / Sync-OFF separation as Case 2)                                                                                                                                                          |
| Case 4, original steps 1–10           | Case 4, migrated steps 1–12         | preserved                                                                                                                                                                                                         |
| Case 5, original steps 1–8            | Case 5, migrated steps 1–10         | preserved                                                                                                                                                                                                         |
| Case 6, original steps 1–8            | Case 6, migrated steps 1–10         | preserved                                                                                                                                                                                                         |
| Case 7, original steps 1–12           | Case 7, migrated steps 1–14         | preserved (original step 7 "Repeat step 6 for join types Outer, Left, Right" preserved as migrated step 7; original steps 8–12 expanded to migrated 8–14 to give each verify/save/reopen its own number)            |
| Case 8, original steps 1–10           | Case 8, migrated steps 1–12         | preserved                                                                                                                                                                                                         |
| Case 9, original steps 1–10           | Case 9, migrated steps 1–12         | preserved                                                                                                                                                                                                         |
| Trailing JSON metadata `{ "order": 1 }` | (dropped from body)                | metadata-not-step (per orchestrator chain analysis convention; original `order: 1` is captured in `scenario-chains/projects.yaml` rev 2 `order_from_files`, not in the migrated body)                              |

No original numbered step is silently dropped; every step has an
explicit migrated counterpart. The renumbering in Cases 2–9 (e.g. 11
original steps → 12 migrated steps) is an axis-clarity expansion: the
original's "Save without Data Sync" sub-flow was expressed as a single
multi-line step in some cases; the migrated expresses each save +
close + reopen + verify cycle as its own numbered step so D-STEP-02
("every Expected result is preserved as a verification step") is
mechanically auditable. No semantic content is added or removed.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 2 `output_plan.uploading.md.target_layer = playwright`,
  consistent with the existing sibling spec
  `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading-spec.ts`
  (already at the playwright layer, covering Cases 1, 3, 4, 9 with
  Sync ON only). Aligning the layer keeps the section's house style
  consistent and avoids cross-layer churn for the remaining 14 paths.
- **Why this `priority`:** chose `regression` per the canonical
  A-STRUCT-MECH-06 enum (`smoke | regression | edge | perf`) — the
  scenario is matrix-coverage of the upload+save+reopen lifecycle
  across 9 source combinations and 2 sync states, NOT a single
  golden-path smoke. NOTE: per the per-cycle priority-enum override
  (decision-log entry `mig-2026-04-30-priority-enum-drift-discovery`),
  the Migrator schema at `grok-migrate-scenario/SKILL.md:222`
  (`priority: p0|p1|p2|p3`) is in drift and was NOT consulted for
  this field; A-STRUCT-MECH-06 at `01-architecture-details.md:988-990`
  was used instead.
- **Why this `strategy`:** `data_driven` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.uploading.md.strategy = data_driven`. The 9×2 matrix
  is the canonical data-driven shape — same step skeleton parameterized
  over 9 table-source combinations and 2 sync states.
- **D-STRUCT-02 compliance:** all 18 paths preserved (9 cases × 2
  sync variants). No reduction applied at the .md level. The existing
  `uploading-spec.ts` covers 4 of 9 cases at Sync ON only (4 of 18
  paths) — that reduction is an Automator-stage decision, not a
  Migrator-stage opt-out. The migrated .md remains the full 18-path
  source of truth.
- **Sibling tests consulted:**
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading-spec.ts` —
    confirms playwright is the convention; Cases 1, 3, 4, 9 are
    already speced (Sync ON only) using a local `saveProject(page, name)`
    helper, `evalJs(page, ...)` for JS-API replay of table-open, and
    `softStep(...)` for soft-error accumulation.
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/upload-project-spec.ts` —
    not present; only `upload-project.md` is the upload golden-path
    sibling (no spec file yet).
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading-run.md` —
    2026-03-09 run summary noting 8 of 14 step-equivalents passed,
    Cases 4–8 skipped due to missing Spaces on the release server.
- **Helpers reused:**
  - At the playwright layer: `saveProject(page, name)` — local helper
    in `uploading-spec.ts:14-23`. Migrator does not invent helpers;
    the body of the migrated .md uses prose ("**File** > **Save
    Project**") rather than a helper invocation, as is convention for
    migrated test-track scenarios.
  - At the grok_test_layer: `uploadProject(projectName, tableInfo, view, df)`
    at `public/packages/UITests/src/gui/gui-utils.ts:100` — exists
    but is NOT playwright-compatible and was NOT chosen as the target
    layer for this scenario. Same flag as in upload-project.md
    migration: a Playwright-layer counterpart helper is a candidate
    addition for `helpers-registry.yaml :: playwright_layer` if Step 2
    (Automator) decides multiple cases of this scenario should share a
    factored helper.
  - **Candidate helper (not invented; flagged for registry):**
    `helpers.playwright.projects.saveAndReopen(page, name, syncOn)` —
    a Playwright-layer counterpart that drives the Save Project
    dialog with a Data Sync toggle, saves under the given name,
    closes, and reopens; would be reused by all 18 paths in this
    scenario plus the 8 paths still pending in scenarios 3, 5, 6.
    Propose to Andrew for addition to `helpers-registry.yaml ::
    playwright_layer` if Step 2 of the orchestrator chain produces
    ≥4 paths that would reuse it (this scenario alone produces 18).
- **Bug library consulted:** yes — `bug-library/projects.yaml` rev 2
  was read. Two bugs intersect this scenario's coverage areas:
  - `GROK-19103` — Join result silently saved as separate project
    that fails to open. Case 7 produces 4 join results saved within
    the same project; not the exact reproduction of the bug (which
    saves the join as a *separate* project), but the upload+save
    flow on a workspace containing joined tables is the same code
    path. Listed in `related_bugs` for traceability.
  - `GROK-18345` — Project built on a Spaces dataset and saved with
    datasync, then shared, fails for the recipient. Cases 4–6 cover
    Spaces+datasync but do NOT exercise the share step; the
    upload+save+reopen portion is the same code path up to the share.
    Listed in `related_bugs` for traceability.
  Other bugs (`GROK-19750`, `GROK-19212`, `GROK-19403`, `GROK-19728`)
  do not intersect this scenario's flows: GROK-19750 is save-copy
  (not save), GROK-19212 is rename-then-save (not save-without-rename),
  GROK-19403 is share-with-unshared-deps (not in this scenario's
  scope), GROK-19728 is view-and-use-failure-state (sharing flow,
  not present here).
- **Decision log queried:** yes — `decision-log.yaml` rev 4 read.
  No prior `migration_decisions`, `layer_decisions`, or `manual_only`
  entry applies retroactively to `uploading.md` content. The
  `mig-2026-04-30-priority-enum-drift-discovery` entry (from scenario 1
  closeout, pending paste) is the relevant cross-reference and was
  honored by writing `priority: regression` per A-STRUCT-MECH-06,
  not `p0..p3` per drifted Migrator SKILL.md.

## Opt-outs (SCOPE_REDUCTION proposals)

(none) — every numbered step of every case is preserved at the same
target layer; no SCOPE_REDUCTION is proposed for the .md migration.
The fact that `uploading-spec.ts` only covers Cases 1, 3, 4, 9 with
Sync ON is a property of the existing spec, not a reduction made by
this migration.

## Deferred items (NOT opt-outs)

- **Cases 4–6 environment dependency on Spaces.** Original
  `uploading-run.md` (2026-03-09) recorded SKIP for Cases 4–8 with
  "No Spaces set up on release server". This is a real environment
  prerequisite, not a content gap in the scenario. Deferred for
  Automator/Validator stage to either (a) configure the SPGIs Space
  in CI before running, or (b) gate Cases 4–6 behind a feature flag
  matching server config. Not a migration-stage decision.
- **`order: 1` collision with `upload-project.md`.** Both
  `upload-project.md` and `uploading.md` carried `{ "order": 1 }` in
  their original trailing JSON metadata. The collision is recorded
  in `scenario-chains/projects.yaml` rev 2 `unresolved_ambiguities`
  and was resolved in chain analysis by alphabetical tie-break
  (upload-project.md before uploading.md). Migrator does not modify
  the chain artifact; this remains a downstream concern if Test
  Track ever surfaces ordering UI based on the original metadata.

## Edge cases

The original lists no explicit edge cases beyond the matrix axes
themselves. Implicit edge cases derivable from the matrix:

- **Pivot Table preserve-on-reopen** (Case 8) — the pivot derivative
  is workspace-only (not synced); reopen must reconstruct it from
  saved project state. Migrated body's Case 8 step 9 verifies "both
  the source table and the pivot table are loaded".
- **Aggregate Rows preserve-on-reopen** (Case 9) — same shape as
  Case 8. Migrated body's Case 9 step 9 verifies the aggregated
  table loads alongside the source.
- **4 join-types preserve-on-reopen** (Case 7) — Inner / Outer /
  Left / Right joins each produce a separate result table, all 5
  tables (source + 4 joins) must reload. Migrated body's Case 7
  step 11 verifies "all tables (source + 4 joined) are loaded".
- **Data Sync toggle persistence** (every case) — the toggle's
  state must be preserved as project metadata so reopen restores
  the same sync behavior. Implicit; verified indirectly by the
  reopen-and-verify-no-console-errors steps.

Feature-wide edge cases for `projects.api.save` and
`projects.api.files.sync` documented in `bug-library/projects.yaml`
rev 2 (GROK-19103, GROK-18345) are listed in this scenario's
`related_bugs` for traceability but are not the exact reproductions
exercised here — they have their own (or near-future) dedicated
scenarios.

(none additional)

## Unresolved ambiguities

- **"Selection to filter" linking-direction asymmetry.** The
  original specifies `Selection to filter` as the link type for all
  Link-Tables cases but only verifies one direction (selection in
  Table 1 → filter in Table 2). The reverse direction (selection in
  Table 2 → filter in Table 1) is not explicitly verified. Migrated
  body preserves the original's one-direction verification.
  Flagged for Gate B if reverse-direction asymmetry produces a
  bug in production.
- **`SPGI_v2_infinity.csv` source ambiguity (Case 1).** The original
  step 2 references `SPGI_v2_infinity.csv` from `System:Demo`, but
  this file is not present in the standard demo bundle as of
  scan-coverage rev 1. Migrated preserves the reference as-is; if
  the file is unavailable in CI, Case 1 fails at step 2. Flag for
  Automator to either substitute a present file or surface this as
  a Gate B environment issue.
- **Project-naming collision risk in CI.** All 18 saved projects
  use deterministic names (`Test_Case<N>_Sync` / `Test_Case<N>_NoSync`).
  Concurrent CI runs would collide. The existing `uploading-spec.ts`
  uses a `Date.now()` suffix (`AutoTest-Upload-Case<N>-<timestamp>`)
  to avoid this — a divergence from the original that the migrated
  preserves the original of. Flag for Automator: replace deterministic
  names with timestamp-suffixed names at spec time, OR clean up
  before each run.
- **"Add to workspace" semantics for Pivot/Aggregate (Cases 8, 9).**
  The original assumes "Add to workspace" is a button on the Pivot
  Table / Aggregate Rows dialog that produces a derivative table in
  a new tab. UI-tier verification of this behavior was not in the
  original; migrated preserves "**Verify:** the pivot table opens as
  a separate table in a new tab" but does not specify the precise
  selector or DOM hook. Flag for Automator: cross-reference
  `grok-browser/references/` for pivot/aggregate dialog patterns at
  spec time.
- **Original step renumbering across Sync-ON / Sync-OFF blocks.**
  Cases 2, 3, 4 in the original use steps 7–9 for Sync ON and
  steps 10–12 for Sync OFF; Cases 5, 6 use steps 5–7 for Sync ON
  and steps 8–10 for Sync OFF. Migrated body normalizes all cases
  to a consistent numbering pattern (Sync ON: save → close+reopen
  → verify; Sync OFF: save → close+reopen → verify) for D-STEP-02
  mechanical auditability. No content change; only numbering
  hygiene.
