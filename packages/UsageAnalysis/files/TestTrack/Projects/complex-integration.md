---
feature: projects
sub_features_covered:
  - projects.upload
  - projects.api.save
  - projects.api.files.sync
  - projects.api.namespaces
  - projects.add-relation
  - projects.add-link
target_layer: playwright
coverage_type: regression
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
migration_report: complex-integration-migration-report.md
related_bugs: []
---

# Complex — Multi-source integration

Source-agnostic op (Wave 2B). The only intentionally-long Wave 2B
spec — value is multi-source co-existence verification. Decomposition
of `complex.md` Step 1 in full + Step 2 (Save with Sync ON + reopen)
per Plan line 170.

**8 sources** opened in one workspace + Pivot + Aggregate + Join + Clone
+ Save Sync ON + reopen. Tests that the platform correctly handles
projects with heterogeneous source coexistence.

UI coverage delegated to `projects-ui-smoke.md`.

## Setup

1. Authenticate as test user.
2. Project name: `integration-test-${Date.now()}`.
3. **Source provisioning** — all required prerequisites are created
   by the spec itself (no external package or DB env dependency):
   - Files: `System:AppData/Chem/tests/spgi-100.csv` and
     `System:DemoFiles/demog.csv` (file presence is not a
     prerequisite per project-wide convention).
   - Space: inline-created `test-projects-demo-integration-*`
     populated with `demog.csv` (Spaces prelude pattern; subject to
     the platform-level Spaces blocker — defensive env-skip per
     `projects-lifecycle-spaces-spec.ts` pattern).
   - DB table: `System:Datagrok / public.groups` opened ad-hoc via
     `helpers/openers.ts:openTableFromDbTable`.
   - Saved query: provisioned in-test via
     `helpers/openers.ts:provisionSystemDatagrokQuery({sql:
     SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE})`.
   - Dataframe-output script: provisioned in-test via
     `helpers/openers.ts:provisionDataframeScript`.
4. Cleanup: delete project; delete inline Space; invoke
   `provisionedQuery.cleanup()` and
   `deleteProvisionedScript(provisionedScript.scriptId)`.

## Scenarios

### Main flow — multi-source integration

0. **Provision in-test prerequisites.** Create the saved query and
   the dataframe-output script via the helpers listed in Setup
   point 3.
1. **Open 8 tables from heterogeneous sources.**
   - File share: `System:AppData/Chem/tests/spgi-100.csv` (via
     `helpers/openers.ts:openTableFromFile`).
   - Space: `test-projects-demo-integration-* > demog.csv`
     (prelude-created; defer / env-skip if Spaces blocker active).
   - DB table: `System:Datagrok / public.groups` opened via
     `openTableFromDbTable` (mirrors `Browse > Databases` double-
     click).
   - Query: run the provisioned saved query via
     `openTableFromDbQuery(page, provisionedQuery.queryNqName)`.
   - Script: run the provisioned script via
     `openTableFromScript(page, provisionedScript.resolvedNqName)`.
   - Pivot: configure Pivot Table on the demog table > **Add**
     (use `addAggregateToWorkspace({via: 'pivot-viewer'})`).
   - Aggregate: configure Aggregate Rows on the script-output
     table > **Add** (use `addAggregateToWorkspace({via: 'menu'})`).
   - Join: Join the spgi-100 table with the DB-table on a
     common-name key column (or any plausible join condition);
     produce as new tab.
   - Clone: right-click any open table > **Clone** to produce
     an independent copy.
   - Verify `grok.shell.tables.length >= 8` after all 8 sources
     are open.
2. **Save project with Data Sync ON (the integration save).**
   Save Project, name from Setup, Data Sync **ON** for every
   applicable table (some tables — derived Pivot/Aggregate/Join,
   Clone — do not have Data Sync configurable at that level;
   their parents do). Click **OK**. Cancel auto-share.
3. **Verify save persisted with all 8 tables.** Verify
   `(await grok.dapi.projects.find(<id>)).relations.length` is
   the expected count (depends on parent vs derived counting
   semantics; assert `>= 5` for parent tables at minimum:
   spgi, demog, db-table, query-result, script-result).
4. **Close all views and reopen the integration project.** Close
   via `grok.shell.closeAll()`. Reopen via
   `Browse > Dashboards`. Verify all parent tables load (5+);
   verify derivations (Pivot, Aggregate, Join, Clone) restore
   correctly. **Multi-source co-existence assertion:** all 8
   tabs are present in the workspace; no missing-source errors.
5. **Cleanup.** Delete the project. Delete the inline Space.
   Invoke `provisionedQuery.cleanup()` and
   `deleteProvisionedScript(provisionedScript.scriptId)`.

### Expected results

- 8 heterogeneous sources can coexist in one project.
- Save Project with Data Sync ON persists all sources +
  derivations.
- Reopen restores all 8 tables / tabs without errors.
- No source-class interaction conflicts (e.g. Spaces +
  System:Datagrok + Files all coexist; Pivot + Aggregate + Join
  all derive correctly; Clone is independent).

## Notes

- **Self-contained source provisioning.** The Query and Script
  prerequisites are created and deleted within the spec; the DB
  table is on the built-in `System:Datagrok` connection. No
  Samples package dependency.
- **Origin: `complex.md` Step 1 in full + Step 2 (Wave 2B
  re-framing per Plan line 170).** Plan explicitly notes:
  "The only intentionally-long spec — value is multi-source
  co-existence."
- **Why this is NOT bug-focused.** None of the GROK bugs
  target multi-source co-existence specifically. Each bug
  is source-specific or op-specific. This scenario covers
  the gap: regressions in multi-source coordination would
  be caught here proactively.
- **`projects.add-link` in sub_features.** Cloning
  exercises `addLink` (vs `addRelation` for fresh sources).
- **Long runtime acknowledged.** Plan explicitly accepts
  this is the longest Wave 2B spec (~10-15 min in CI).
  Value is in multi-source coverage; signal/noise improves
  vs running 8 separate single-source specs.
- **UI coverage delegated.** No UI surface beyond Save
  Project. Pivot/Aggregate/Join/Clone UI surfaces are owned
  by `uploading.md` source-matrix Cases 8/9.
- **Spaces blocker.** Cases involving Spaces (Step 1
  sub-bullet 2) defensively env-skip via the same pattern as
  `projects-lifecycle-spaces-spec.ts` until the platform-level
  Spaces bug is fixed.
- **No `related_bugs`.** Pure proactive multi-source
  coverage.
- **Self-cleaning.** Step 5 deletes project + Space + provisioned
  query + provisioned script.
