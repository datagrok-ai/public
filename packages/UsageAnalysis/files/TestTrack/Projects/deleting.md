---
feature: projects
sub_features_covered: [projects.api.search, projects.api.delete]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/deleting.md
migration_date: 2026-04-30
migration_report: deleting-migration-report.md
related_bugs: []
---

# Deleting

Terminal cleanup of every project produced by preceding scenarios in
this chain. For each project produced upstream, locate it in Browse >
Dashboards, trigger Delete via the right-click context menu (or via
the Context Panel name-dropdown alternative), confirm in the dialog,
and verify the project disappears from Dashboards.

This scenario is **must_run_last: true** per
`scenario-chains/projects.yaml` rev 2 — placing it earlier in the
chain would invalidate the rest of the run by deleting projects that
downstream scenarios consume. Strategy is `end_to_end_fixtures` per
`output_plan.deleting.md.strategy = end_to_end_fixtures`: the chain's
upstream scenarios (1-8) collectively produce the fixture set
consumed here for cleanup.

## Setup

1. **Fixture prerequisite (from upstream chain):** every project
   produced by preceding scenarios must exist on the server at
   the start of this scenario. The complete consumed set per chain
   rev 2 `dependency_graph.deleting.md.consumes`:
   - `demog` (from `upload-project.md`)
   - 18 `Test_Case<N>_Sync` / `Test_Case<N>_NoSync` for N in 1..9
     (from `uploading.md`)
   - shared `demog` variants (from `share-project.md` — same
     `demog` project after sharing; not a new project)
   - 3 copy-clone-customizations variants:
     `demog-link`, `demog-clone`,
     `demog-personal-view-customizations` (from
     `projects-copy-clone.md`)
   - up to 4 projects from `complex.md` (steps 2/5/6/8; exact names
     depend on Automator's spec-time naming; e.g.
     `Complex_Test_Project_1`, `Complex_Test_Project_1_copy_NoSync`,
     `Complex_Test_Project_1_copy_Sync`,
     `Complex_Test_Project_1_renamed_Sync`)
   - 1 project from `custom-creation-scripts.md`
     (`Custom_Creation_Script_Test`)

   Total: ~26-30 projects depending on upstream chain configuration
   and which optional variants ran. The Automator's `beforeAll`
   block is responsible for collecting the actual produced-project
   list at runtime (e.g. via
   `grok.dapi.projects.filter('author = currentUser AND createdOn > chainStart').list()`)
   rather than relying on a hardcoded list, since exact names depend
   on chain-state-at-runtime.
2. **Authentication:** session must be authenticated as the project
   owner (delete authority required).
3. **No additional environment dependencies** beyond what upstream
   scenarios already required.

## Scenarios

### Delete each upstream-produced project and verify removal

The following workflow is parameterised over the consumed fixture
list. Implemented per `end_to_end_fixtures` strategy at the
Playwright layer (per `output_plan.deleting.md.strategy =
end_to_end_fixtures` in `scenario-chains/projects.yaml` rev 2).

For each `<projectName>` in the upstream-produced project list:

1. Go to **Browse** > **Dashboards** and locate `<projectName>`.
   Use the search bar or scroll/filter as needed.
2. Trigger the **Delete project** action. **Either**:
   - **Option A (right-click context menu):** Right-click
     `<projectName>` in the Dashboards listing → select **Delete
     project** from the context menu. A confirmation dialog opens.
   - **Option B (Context Panel dropdown):** Click `<projectName>`
     to select; on the Context Panel, click the dropdown next to
     the project name → select **Delete**. A confirmation dialog
     opens.
3. In the confirmation dialog, click **DELETE** to confirm. (The
   DELETE button has no `name=` attribute; locate by visible text
   per `grok-browser/references/projects.md` line 153.)
4. **Verify on Dashboards:** the project no longer appears in
   Browse > Dashboards (search/scroll confirms absence). Also
   confirm via JS API:
   `(await grok.dapi.projects.filter('name = "<projectName>"').list()).length === 0`.

## Notes

- **Must-run-last invariant** (line 1 of source: "Note: make sure
  this is the last test case"). The chain orchestrator
  (`scenario-chains/projects.yaml` rev 2) enforces this via
  `dependency_graph.deleting.md.must_run_last: true`. Any chain
  rearrangement that places `deleting.md` before another scenario
  is invalid by construction — the orchestrator must reject it.
- **Original `order: 6`** — captured in
  `scenario-chains/projects.yaml` rev 2 `order_from_files` for
  TestTrack-display purposes. NOT the chain execution order
  (`must_run_last: true` overrides numeric ordering for chain
  semantics).
- **No bug intersections.** Per
  `decision-log.yaml :: migration_decisions` entry
  `mig-2026-04-29-bug-removed`, the historical github-3752 bug
  ("Deleting a project fails with FK violation on
  view_layouts_columns / table_columns") was removed from
  `bug-library/projects.yaml` (rev 2) per 2026-04-29 ("bad
  exemplar"). The atlas critical_path
  `delete-with-view-layouts` that originally provided regression
  coverage for this bug was also deleted in that cleanup. As a
  result, this scenario carries `related_bugs: []` — no
  bug-library entries intersect plain delete-flow as of rev 2.
- **Per-cycle Invariant 3 (atlas-aware sub_features for share)
  does NOT apply.** No Share dialog operation in this scenario.
  `projects.shell.share-via-context-menu` correctly OMITTED from
  `sub_features_covered`.
- **Per-cycle Invariant 2 (existing -spec.ts / -api.ts READ-ONLY)
  applies and is satisfied.** An existing `deleting-spec.ts`
  exists at
  `public/packages/UsageAnalysis/files/TestTrack/Projects/deleting-spec.ts`
  (read-only sibling-test-convention reference for this
  migration; not modified). Its content reflects pre-migration
  deleting.md semantics; spec regeneration to match the migrated
  body is Step 7 of Phase 1 plan, owned by Automator (not
  Migrator).
- **Cleanup completeness as the chain's terminal contract.**
  After this scenario completes, the test environment should be
  free of every project produced by the chain. If a downstream
  retro reports leftover projects, this scenario's parameter list
  was incomplete (Automator's `beforeAll` filter under-collected
  produced projects). Surfaced as a deferred item in the
  migration report.
