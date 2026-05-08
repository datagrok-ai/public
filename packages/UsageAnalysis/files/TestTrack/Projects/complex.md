---
feature: projects
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.shell.open, projects.add-relation, projects.api.namespaces, projects.shell.share-via-context-menu]
target_layer: playwright
coverage_type: regression
pyramid_layer: bug-focused
ui_coverage_responsibility:
  - context-menu-rename-project
  - context-menu-rename-query
  - context-menu-rename-script
  - pcmdShareProject
  - share-dialog-permissions-editor
  - logout-login-as-second-user
  - data-sync-refresh-verification
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
migration_report: complex-migration-report.md
ui_companion: complex-ui.md
realized_as:
  - complex-derived-tables-spec.ts
  - complex-rename-spec.ts
  - complex-share-second-user-spec.ts
related_bugs: [GROK-19212, GROK-19103, GROK-19403, GROK-18345, GROK-19728, github-3550]
---

# Complex

End-to-end mega-scenario exercising the full lifecycle of a project
across many data sources (file share, Space, database, query, script,
pivot, aggregate, join), Data Sync mode toggles, drag-and-drop
augmentation, table rename inside a project, entity rename (Project /
Query / Script), entity move (to file share, then to Space), share with
a second user at two access levels (View-and-Use, Full), and re-auth
verification as that second user.

This scenario is **complex-standalone** per
`scenario-chains/projects.yaml` rev 3: no cross-file fixtures consumed;
all required state is created inline. The "second user" identity is a
LITERAL placeholder per decision-log entry
`mig-2026-04-29-fixture-placeholder` (eventual fixture account TBD; not
bound to any named real account).

This scenario is **the central cross-cutting bug exemplar** for the
Projects chain. Per chain YAML rev 3 `pyramid_layer: bug-focused` (Rule
3 — 6 bugs walked here): Steps 7-9 (rename) + Step 11 (Data Sync
verify) walk the GROK-19212 reproduction path (table-rename-then-
reopen reference resolution under datasync); Step 1 sub-bullet "Join
two tables" + Step 2 (Save with Sync) walks GROK-19103's repro (join
result silently saved as separate project that fails to open); Step 1
"Table from Space" + Step 2 (Save with Sync) + Step 12 (Share with
second user) + Step 13 (open as second user) walks GROK-18345's repro
(recipient cannot open shared project that uses Spaces dataset saved
with data sync); Step 9 sub-bullet "Query — rename" + Step 11 (open
and verify) walks github-3550's repro (external-entity rename
invalidation, sister of GROK-19212); Step 12 (Share at View-and-Use)
+ Step 13 (recipient open) intersects GROK-19403 (recipient cannot
open shared project when underlying script not also shared) at the
share + recipient-open touchpoint; Step 12 (Share at View-and-Use) +
Step 13 (recipient access under failure state) intersects GROK-19728
(view-and-use users edit creation script under failure state) at the
view-and-use + failure-state surface. **GROK-18345 is the primary
motivating example for chain-analyzer-prompt.md Edit 5 Trigger 2 (single-
scenario non-adjacent steps): the bug spans Step 1 + Step 2 + Step 12 +
Step 13 — non-adjacent matched steps within this single scenario.**

This scenario was **already split** at Wave 1b into 3 satellite specs
per decision-log `wave-1b-complex-split-b70-followup` (2026-05-01):

- `complex-derived-tables-spec.ts` — Step 1 multi-source open + Join
  sub-bullet + Step 2 Save with Sync; targets GROK-19103.
- `complex-rename-spec.ts` — Steps 7-8 table rename + Step 11 reopen-
  verify; targets GROK-19212; partial github-3550.
- `complex-share-second-user-spec.ts` — Step 12 share at View-and-Use +
  Full; targets GROK-18345 (Step 13 deferred pending Helper 3
  `helpers.playwright.session.logoutAndLoginAs` registration).

The `complex-move` sub-spec was **dropped** per decision-log
`wave-1b-complex-split-b70-followup` and reinforced by
`b2-2026-05-03-drag-drop-ui-only-reclassification` — both documented
paths (drag-drop and right-click `Move to`) are unautomatable in
current state. The original 13-step .md remains the canonical scope
reference; realized-state coverage is tracked in scenario-chains rev 3
`output_plan.complex.md.realized_as`.

## Setup

1. **"Second user" identity placeholder** — step 12 shares the project
   with "another user" at View-and-Use AND Full access levels; step 13
   logs in as that user. Per decision-log entry
   `mig-2026-04-29-fixture-placeholder`, the literal string "second
   user" is the placeholder identity in this scenario; eventual fixture
   account binding (a specific email + credentials, OR a per-run
   timestamp-suffixed test account) is TBD and is the Automator
   stage's responsibility.
2. **Re-auth pattern dependency** — step 13 ("Log in as the second
   user") requires a logout + login-as-different-credentials + verify
   + (optionally) login-back-as-original sequence. Per decision-log
   `b14-2026-04-30-re-auth-pattern-applied`, the re-auth pattern was
   added to `.claude/skills/grok-browser/references/projects.md` (UI
   path documented; JS API path TODO awaiting external confirmation).
   Helper candidate `helpers.playwright.session.logoutAndLoginAs`
   surfaced; not yet registered (Helper 3 from helpers-batch-1
   backlog). Current realized coverage in
   `complex-share-second-user-spec.ts` defers Step 13 pending the
   helper.
3. **Source provisioning** — the scenario uses these table sources:
   `System:AppData/Chem/tests/spgi-100.csv` (file share), an inline-
   created `test-projects-demo` Space populated with `demog.csv`
   from `System:DemoFiles` via JS API (see Setup point 4), an
   in-test-provisioned saved query on `System:Datagrok` (created via
   `helpers/openers.ts:provisionSystemDatagrokQuery({sql:
   SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE})`), an in-test-provisioned
   dataframe-output script (created via
   `helpers/openers.ts:provisionDataframeScript`), and an ad-hoc DB
   table on `System:Datagrok` (`public.groups` — opened via
   `openTableFromDbTable`). All sources are self-contained; no
   external package (Samples) or env-provisioned DB connection is
   required.
4. **Spaces prelude / postlude** — before step 1's Spaces sub-bullet,
   create the Space via JS API:
   ```js
   const space = await grok.dapi.spaces.createRoot('test-projects-demo');
   const client = await grok.dapi.spaces.spaceClient(space.id);
   const file = (await grok.dapi.files.list('System:DemoFiles', false, 'demog.csv'))[0];
   await client.addEntity(file, /*link=*/true);
   ```
   After the scenario completes (post step 13), delete the Space:
   ```js
   await grok.dapi.spaces.delete(space);
   ```

## Scenarios

### Mega-flow: open many sources, save, evolve, rename, move, share, re-auth verify

0. **Provision in-test prerequisites.** Before opening sources, the
   spec creates: (a) a saved query on `System:Datagrok` via
   `helpers/openers.ts:provisionSystemDatagrokQuery`, (b) a
   dataframe-output JS script via
   `helpers/openers.ts:provisionDataframeScript`. Both are namespaced
   under the test user's login and cleaned up in `finally`.
1. **Open tables from different sources** — open at least one table
   each from the following 8 sources, ensuring all open tables
   coexist in the workspace:
   - Table from a file share (`System:AppData` > `Chem` > `tests` >
     `spgi-100.csv`).
   - Table from a Space (`Browse` > `Spaces` > `test-projects-demo` >
     `demog.csv` — Space created via JS API per Setup point 4).
   - Table from a database — open `public.groups` on
     `System:Datagrok` via
     `helpers/openers.ts:openTableFromDbTable(page, {connectionNqName:
     SYSTEM_DATAGROK_NQNAME, schemaName: 'public', tableName:
     'groups'})` (mirrors `Browse` > `Databases` double-click
     semantics).
   - Table from a Query — run the provisioned query via
     `helpers/openers.ts:openTableFromDbQuery(page,
     provisionedQuery.queryNqName)`.
   - Table generated by a Script — run the provisioned script via
     `helpers/openers.ts:openTableFromScript(page,
     provisionedScript.resolvedNqName)`.
   - Add tables via **Pivot Table** > **Add** (configure rows / cols
     / values on an existing source table, click **Add** to add as
     new tab) AND via **Data** > **Aggregate Rows** > **Add** (group
     + aggregate on another source table, click **Add**) — see
     `helpers/openers.ts:addAggregateToWorkspace`.
   - Join two tables — **Data** > **Join Tables** with one DB-source
     table + one query-result table; pick a key column on each;
     produce the join result as a new tab.
   - Clone a table opened from a script — right-click the script-
     output table > **Clone** (or equivalent UI clone action) to
     produce an independent copy in a new tab.
2. **Save all opened tables as a project with Data Sync enabled.**
   Trigger **File** > **Save Project**. In the Save Project dialog,
   ensure the **Data Sync** toggle is **ON** for every applicable
   table. Name the project (e.g. `Complex_Test_Project_1`). Click
   **OK**. Wait for the save to complete (verify project appears in
   `Browse` > `Dashboards`, OR the underlying `POST /projects`
   succeeds). Cancel the auto-opened Share dialog if it appears.
3. **Open tables from 4 additional sources** — open new tables not
   already in the project from each:
   - Spaces (e.g. `test-projects-demo` > re-open `demog.csv` for a
     fresh table tab; auto-disambiguated as `demog (2)`).
   - Files (e.g. `System:AppData` > `Chem` > `tests` > re-open
     `spgi-100.csv` for a fresh table tab).
   - Query result — re-run the provisioned query for a fresh result
     OR provision a second saved query (e.g. with
     `SYSTEM_DATAGROK_QUERIES.GROUPS_RELATIONS`).
   - DB table — open a different table on `System:Datagrok` (e.g.
     `groups_relations`) via `openTableFromDbTable`.
4. **Add newly opened tables to the opened project (drag-and-drop in
   the Dashboards).**

   > **UI-only step moved to `complex-ui.md`** (preserving original
   > step number 4 for cross-reference). Reason: drag-drop event
   > registration not automatable through current Playwright + JS API
   > mechanisms (3 mechanisms tried; see decision-log
   > `b2-2026-05-03-drag-drop-ui-only-reclassification`). Section
   > excluded from `-spec.ts` generation.
5. **Open the project and save a copy without Data Sync for all
   tables.** Open the project from step 2 (now augmented with 4
   tables from step 4). Trigger **Save Copy** (or **File** > **Save
   Project As**). In the Save dialog, turn the **Data Sync** toggle
   **OFF** for every table. Name the copy (e.g.
   `Complex_Test_Project_1_copy_NoSync`). **OK**. Cancel any
   auto-share dialog.
6. **Open the copied project and save it again with Data Sync
   enabled** *(result: two projects, both with Data Sync ON)*.
   Open `Complex_Test_Project_1_copy_NoSync` from step 5. Trigger
   **File** > **Save Project** (re-save under the same name OR
   under a new name like `Complex_Test_Project_1_copy_Sync`). In
   the Save dialog, turn Data Sync **ON** for every table. **OK**.
   At this point: the original project from step 2 has Data Sync
   ON; the copy from step 5 (re-saved here) has Data Sync ON. Two
   projects, both Sync ON.
7. **Open the last saved project and rename all tables inside it.**
   Open the project from step 6 (the re-saved copy with Sync ON).
   For each table in the project, rename it via the table's tab
   context menu OR via `grok.shell.tables[i].name = '<new>'` (or
   the equivalent rename UI). Verify the renamed tables persist.
8. **Save a copy of the project with Data Sync enabled.** Trigger
   **Save Copy** on the project from step 7 (with renamed tables).
   In the Save dialog, Data Sync **ON**. Name the copy (e.g.
   `Complex_Test_Project_1_renamed_Sync`). **OK**.
9. **Rename the following entities** (entities outside the project,
   referenced BY the project):
   - **Project** — rename the project from step 8 (e.g. add a `_v2`
     suffix). Trigger via `Browse` > `Dashboards` > right-click >
     **Rename**, OR via Context Panel.
   - **Query** — rename the provisioned query used as a data source
     (e.g. `<resolvedName>` → `<resolvedName>_renamed`) via JS API:
     ```js
     const q = await grok.dapi.queries.find(provisionedQuery.queryId);
     q.name = `${provisionedQuery.resolvedName}_renamed`;
     await grok.dapi.queries.save(q);
     ```
     The test owns the query (it's namespaced under the test user's
     login) — rename always succeeds. No env-permission caveat.
   - **Script** — rename the provisioned script used as a data
     source via JS API:
     ```js
     const s = await grok.dapi.scripts.find(provisionedScript.scriptId);
     s.name = `${provisionedScript.resolvedName}_renamed`;
     await grok.dapi.scripts.save(s);
     ```
     The test owns the script — rename always succeeds.
10. **Move the following entities to any file share, then to any
    Space.**

    > **UI-only step moved to `complex-ui.md`** (preserving original
    > step number 10 for cross-reference). Reason: both documented
    > paths are unavailable for automation — drag-drop registration
    > mechanism not accessible from Playwright (same blocker as
    > Step 4) AND the right-click `Move to` menu option does not
    > exist in the current UI (verified 2026-05-03 against
    > dev.datagrok.ai). See decision-log
    > `b2-2026-05-03-drag-drop-ui-only-reclassification`. Section
    > excluded from `-spec.ts` generation.
11. **Open the moved project and previously saved project and
    verify:**
    - All tables are available (no missing-table errors on open).
    - Table relationships are preserved (linked tables still
      filter / select-propagate correctly).
    - Data Sync updates correctly (force a Data Sync refresh OR
      re-open after a server-side data change; verify the project's
      tables reflect the latest data).
12. **Configure project sharing for another user** — right-click the
    moved project in `Browse` > `Dashboards`, select **Share**, fill
    "second user" as the recipient (literal placeholder per
    decision-log `mig-2026-04-29-fixture-placeholder`):
    - **View-and-Use access** — grant view-and-use level (no edit
      permission).
    - **Full access** — additionally OR alternatively grant Full
      (edit) level. The original is silent on whether both levels
      are granted simultaneously to one second user OR to two
      different second users; preserved as ambiguity (see migration
      report).
13. **Log in as the second user and open the shared project.**
    Logout from the current session; log in as the "second user"
    placeholder identity (eventual fixture account TBD per Setup
    point 1; re-auth pattern dependency per Setup point 2);
    navigate to `Browse` > `Dashboards`; locate the shared project
    (it should appear in the second user's accessible projects,
    possibly under "Shared with me" or a similar grouping); open
    it.

### Expected results

After completing all 13 steps, verify:

- **All entities function correctly after renaming and moving.**
  Renamed Project / Query / Script remain accessible at their new
  names; the project's references resolve correctly post-rename;
  moved entities are accessible at their new namespace location.
- **Data Sync updates as expected.** Both Sync-ON projects (from
  steps 2 and 6/8) refresh their data correctly when underlying
  data sources change.
- **Users with different access levels see the project and its
  contents according to their permissions.** The "second user"
  with View-and-Use access can OPEN the project but cannot edit;
  with Full access can edit; both can see the project's tables and
  viewers per their permissions.

## Notes

- **Original `order: 8`** — runs after all preceding scenarios per
  `scenario-chains/projects.yaml` rev 3 `order_from_files`. The
  chain ordering interleaves complex.md after deleting.md
  (`order: 7` per source) but `must_run_last: false` for complex.md
  (only `deleting.md` is the chain's terminal `must_run_last: true`).
- **Realized as 3 satellite specs (Wave 1b split — preserved).** Per
  decision-log `wave-1b-complex-split-b70-followup` (2026-05-01),
  this single .md has been realized as three satellite specs:
  `complex-derived-tables-spec.ts` (GROK-19103),
  `complex-rename-spec.ts` (GROK-19212; partial github-3550),
  `complex-share-second-user-spec.ts` (GROK-18345 share-side; Step 13
  deferred). The `complex-move` sub-spec was dropped per
  `b2-2026-05-03-drag-drop-ui-only-reclassification` and is NOT to be
  re-merged. This .md's frontmatter intentionally describes the full
  canonical 13-step scope; the realized coverage is the 3-spec union
  documented in `output_plan.complex.md.realized_as`.
- **Pyramid layer: bug-focused (chain rev 3 Rule 3).** 6 bugs walked
  by this scenario: GROK-19212 (rename), GROK-19103 (join+save),
  GROK-19403 (share with un-shared deps), GROK-18345 (Spaces+sync+
  share), GROK-19728 (view-and-use under failure state), github-3550
  (external Query/Script rename invalidation). All 6 listed in
  `related_bugs` per chain YAML rev 3.
- **GROK-18345 single-scenario non-adjacent (Trigger 2 motivating
  example).** GROK-18345 spans **Step 1** (Spaces sub-bullet),
  **Step 2** (Save with Sync ON), **Step 12** (Share with second
  user), and **Step 13** (recipient opens shared project) — exactly
  the 4-component repro path of the bug, all within this single
  scenario. This is THE primary motivating example for chain-
  analyzer-prompt.md Edit 5 Trigger 2 (single-scenario non-adjacent
  steps). Spec lines explicitly cite complex.md case in Edit 5
  rationale.
- **Cross-cutting candidate specs (chain YAML rev 3 bug_focused_candidates).**
  Per migration-prompt.md "Cross-cutting bug citations from chain
  YAML": this scenario is a span in the following candidates:
  - GROK-19212: cross-cutting candidate spec —
    `projects-grok-19212-spec.ts` (chain rev 3 spans uploading.md:
    Step 7 + complex.md:Step 7).
  - GROK-19103: cross-cutting candidate spec —
    `projects-grok-19103-spec.ts` (chain rev 3 spans upload-project.md:
    Step 1 + uploading.md:Step 6 + projects-copy-clone.md:Step 4 +
    complex.md:Step 1).
  - GROK-19403: cross-cutting candidate spec —
    `projects-grok-19403-spec.ts` (chain rev 3 spans share-project.md:
    Step 4 + complex.md:Step 12 + projects-copy-clone.md:Step 5).
  - GROK-18345: cross-cutting candidate spec —
    `projects-grok-18345-spec.ts` (chain rev 3 spans uploading.md:
    Step 7 + complex.md:Step 12).
  - GROK-19728: cross-cutting candidate spec —
    `projects-grok-19728-spec.ts` (chain rev 3 spans complex.md:
    Step 12 — single-scenario; emitted under Trigger 1 with
    lifecycle-api.md inclusion via semantic step match per chain
    rev 3 rationale).
  - github-3550: cross-cutting candidate spec —
    `projects-github-3550-spec.ts` (chain rev 3 spans complex.md:
    Step 9).
  Per-scenario migration alone does not capture these bugs'
  cross-cutting invariants; F-BUG-COVERAGE-01 at section-complete
  is the authoritative gate.
- **"Second user" placeholder.** Per
  `mig-2026-04-29-fixture-placeholder`, the literal string "second
  user" is used as the recipient identity in step 12 and as the
  re-auth target identity in step 13. Eventual binding to a specific
  test account is TBD; the placeholder is intentional. Decision-log
  cross-reference recorded; not flagged as a Direct-answer gap.
- **Re-auth pattern dependency.** Per decision-log
  `b14-2026-04-30-re-auth-pattern-applied`, the re-auth pattern was
  added to `.claude/skills/grok-browser/references/projects.md`
  (UI path documented; JS API path TODO). Helper candidate
  `helpers.playwright.session.logoutAndLoginAs` not yet registered
  (Helper 3 backlog). Step 13 in `complex-share-second-user-spec.ts`
  is currently deferred pending registration.
- **Invariant 3 (atlas-aware sub_features for share) APPLIES.**
  Step 12 exercises the right-click Share dialog
  (`pcmdShareProject` per atlas rev 11 sub_feature
  `projects.shell.share-via-context-menu`). The sub_feature is
  correctly INCLUDED in `sub_features_covered`.
- **UI coverage owned (chain rev 3 ui_coverage_responsibility, 7
  flows).** This scenario is the chain's witness for: rename
  context-menu surface (`context-menu-rename-project`,
  `context-menu-rename-query`, `context-menu-rename-script` —
  Step 9), share dialog with PermissionsEditor at two access levels
  (`pcmdShareProject`, `share-dialog-permissions-editor` — Step 12),
  re-auth flow (`logout-login-as-second-user` — Step 13),
  data-sync refresh verification (`data-sync-refresh-verification` —
  Step 11). `ui_coverage_delegated_to: null` — these flows are not
  delegated to any sibling scenario in the chain.
- **Cleanup responsibility.** This scenario produces multiple saved
  projects (from steps 2, 5, 6, 8) and possibly leaves entities
  renamed/moved. Terminal cleanup is the `deleting.md` (scenario 9,
  must_run_last) responsibility for the produced projects;
  rename/move undoing is NOT covered by `deleting.md` and may need
  Automator's `afterAll` to restore entity names/locations OR
  isolated-environment teardown. Surfaced as a deferred item in
  the migration report.
- **Sibling spec convention.** The Wave 1b split realized this .md
  as 3 satellite specs (see "Realized as 3 satellite specs" note
  above). The `complex-run.md` companion exists from a prior MCP
  reproduction run but has no spec counterpart. `complex-ui.md`
  carries the UI-only Step 4 + Step 10 content per
  `b2-2026-05-03-drag-drop-ui-only-reclassification`.
