---
feature: projects
sub_features_covered: [projects.view.browse, projects.api.search, projects.shell.open]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/opening.md
migration_date: 2026-04-30
migration_report: opening-migration-report.md
ui_companion: opening-ui.md
related_bugs: []
---

# Opening

For each of the 18 saved projects produced by `uploading.md`
(`Test_Case<N>_Sync` and `Test_Case<N>_NoSync` for N in 1..9), navigate
to Browse > Dashboards, locate the project, verify its Context Panel
attributes render correctly, then double-click to open and verify the
project loads. Close all views between iterations.

This scenario depends on `uploading.md` having produced its 18 projects
(consumes the `uploading-18-projects` fixture per scenario-chains rev 2).

## Setup

1. **Fixture prerequisite:** all 18 projects produced by `uploading.md`
   must exist on the server with their canonical names
   (`Test_Case1_Sync`, `Test_Case1_NoSync`, …, `Test_Case9_Sync`,
   `Test_Case9_NoSync`). In `chained_tests` strategy, the Automator
   stage will set this up via a `beforeAll` block that either
   (a) reuses the projects produced by a prior `uploading.md` test
   run in the same chain, OR (b) builds the fixture on demand via
   `js-api-replay` of the relevant `uploading.md` cases.
2. Each test case in the chain (one per project name) starts from a
   clean Datagrok session — `beforeEach` calls `grok.shell.closeAll()`.

## Scenarios

### Open each Uploading project and verify Context Panel + open

The following workflow is parameterised over the 18 project names from
the `uploading-18-projects` fixture. Implemented as a `describe.each`
over the project list at the Playwright layer (per `chained_tests`
strategy in `scenario-chains/projects.yaml` rev 2).

For each `<projectName>` in the 18-project list:

1. Go to **Browse** > **Dashboards**.
2. Locate `<projectName>` in the Dashboards listing (search bar or
   scroll/filter).
3. Click `<projectName>` to select it (single-click).
4. **Verify on Context Panel — Name attribute renders:**
   - Name — `<projectName>` rendered correctly.

   > **UI-only Step 4 sub-bullets moved to `opening-ui.md`** (preserving
   > original sub-bullet labels for cross-reference): Sharing,
   > Description, Picture — render-quality verifications requiring
   > human QA judgment per D3 bucket-b classification.
5. Double-click `<projectName>` to open it.
6. **Verify on open:** the project loads, its tables and viewers
   render in the workspace, no errors in the browser console (F12).
7. **Close All** open views and viewers (cleanup before next
   iteration).

### Project name list (18 paths)

Per `chained_tests` strategy, the parameter list is:

| # | Project name |
|---|--------------|
| 1 | `Test_Case1_Sync` |
| 2 | `Test_Case1_NoSync` |
| 3 | `Test_Case2_Sync` |
| 4 | `Test_Case2_NoSync` |
| 5 | `Test_Case3_Sync` |
| 6 | `Test_Case3_NoSync` |
| 7 | `Test_Case4_Sync` |
| 8 | `Test_Case4_NoSync` |
| 9 | `Test_Case5_Sync` |
| 10 | `Test_Case5_NoSync` |
| 11 | `Test_Case6_Sync` |
| 12 | `Test_Case6_NoSync` |
| 13 | `Test_Case7_Sync` |
| 14 | `Test_Case7_NoSync` |
| 15 | `Test_Case8_Sync` |
| 16 | `Test_Case8_NoSync` |
| 17 | `Test_Case9_Sync` |
| 18 | `Test_Case9_NoSync` |

## Notes

- **Original `order: 3`** — runs after `upload-project.md` (`order: 1`)
  and `uploading.md` (`order: 1`). Captured in
  `scenario-chains/projects.yaml` rev 2 `dependency_graph`
  (`opening.md.depends_on: [uploading.md]`).
- **"Projects from the previous step (Uploading)"** in the original
  source resolves to the 18 named projects produced by `uploading.md`
  per chain-analysis rev 2. The migrated body enumerates them
  explicitly so the parameter list is auditable and the
  `chained_tests` strategy has a concrete table of test inputs.
- **"Sharing" Context Panel attribute** in step 4 is the rendered
  metadata about who the project is shared with — NOT the right-click
  Share dialog operation. The Sharing tab/section will be empty for
  solo-owner projects (no recipients yet); a non-empty list appears
  only after `share-project.md` or similar has run. The migrated body
  preserves the original's verification of the rendering, not of any
  share-completion semantics. Invariant 3 of the per-cycle override
  (atlas-aware sub_features_covered for share) does NOT apply here —
  this scenario does not exercise the Share dialog.
- **Cleanup responsibility:** `closeAll` between iterations and at
  the end of the chain. The 18 projects themselves are NOT deleted
  by this scenario — `deleting.md` (scenario 9, `must_run_last`)
  owns terminal cleanup of all projects produced upstream.
- **Sibling spec convention:** existing `opening-spec.ts` exists at
  `public/packages/UsageAnalysis/files/TestTrack/Projects/opening-spec.ts`
  but covers only ONE project (`AutoTest-Opening-<timestamp>`) rather
  than iterating over the 18 from `uploading.md`. The existing spec's
  one-project reduction is an Automator-stage decision; this migrated
  `.md` preserves the full 18-path chain per D-STRUCT-02. The existing
  spec is **READ-ONLY for Migrator** in this cycle (Invariant 2 of the
  per-cycle override) and was not touched.
