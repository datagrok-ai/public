---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [rename-dependent-entity-reopen, share-with-unshared-deps, share-spaces-datasync, view-and-use-failure-state, derive-then-save-inside-project]
realizes: []
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
source_text_fixes: []
candidate_helpers:
  - helpers.playwright.projects.dragDropOntoDashboard
  - helpers.playwright.session.logoutAndLoginAs
  - helpers.playwright.projects.move
  - helpers.playwright.projects.rename
unresolved_ambiguities:
  - another-user-second-user-single-identity-vs-two
  - order-of-access-level-grants-in-step-12
  - table-from-a-script-semantics-in-step-1
  - pivot-table-add-and-aggregate-rows-add-ui-controls
  - move-to-any-file-share-then-to-any-space-in-step-10
  - rename-ui-access-points-in-step-9
  - verify-data-sync-updates-correctly-in-step-11
  - grok-19728-in-chain-rev-3-emission-rationale
scope_reductions: []
ui_companion: complex-ui.md
realized_as:
  - complex-derived-tables-spec.ts
  - complex-rename-spec.ts
  - complex-share-second-user-spec.ts
related_bugs: [GROK-19212, GROK-19103, GROK-19403, GROK-18345, GROK-19728, github-3550]
---

# Complex — Multi-source project lifecycle mega-scenario

End-to-end scenario exercising the full lifecycle of a project: opening
tables from eight different sources (file share, Space, database,
saved query, script, pivot table, aggregate rows, and a table join),
saving with Data Sync on and off, renaming tables and the entities
that produced them (project, query, script), moving entities across
namespaces, sharing with a second user at two access levels, and
verifying access as that second user.

This is the platform's primary regression scenario for the Projects
area — it walks the reproduction path of six known bugs together:
table rename losing its data-sync link on reopen (GROK-19212), a
joined table silently saved as a separate, broken project
(GROK-19103), a shared project failing to open because its
Spaces-sourced table can't resolve under the recipient's identity
(GROK-18345), a shared project failing to open because its underlying
script wasn't also shared (GROK-19403), a view-and-use user able to
edit the creation script while it's in a failure state (GROK-19728),
and an external query rename invalidating a project reference
(github-3550).

Because a single 13-step scenario this size is unwieldy to automate
and debug as one script, it is realized as three separate Playwright
specs — `complex-derived-tables-spec.ts`, `complex-rename-spec.ts`,
and `complex-share-second-user-spec.ts`. The "second user" identity
used for sharing is a literal placeholder — an eventual fixture
account is still to be decided, not bound to any named real account.
Drag-and-drop steps (adding tables to an open project in Step 4;
moving entities between namespaces in Step 10) are not automatable
through Playwright and are covered instead by the manual companion
`complex-ui.md`.

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
    > excluded from `-spec.ts` generation. The JS-API move path is
    > covered separately by `complex-move.md` /
    > `complex-move-spec.ts`.
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

- **Realized as 3 satellite specs; move + drag-drop split out.** This
  single .md is realized as three Playwright specs:
  `complex-derived-tables-spec.ts` (targets GROK-19103),
  `complex-rename-spec.ts` (targets GROK-19212; partially covers
  github-3550), and `complex-share-second-user-spec.ts` (targets
  GROK-18345; Step 13's recipient-login is deferred pending a
  logout/login-as-other-user helper). The Step 10 move across
  namespaces was decomposed into the separate scenario
  `complex-move.md`, realized by `complex-move-spec.ts` via the
  `grok.dapi.projects` JS API; only the UI move paths (drag-drop and
  right-click "Move to") remain unautomatable in the current UI and
  live in the manual companion `complex-ui.md`. This .md remains the
  canonical 13-step scope; the three specs together are its realized
  coverage.
- **Cross-cutting bug coverage lives partly in sibling scenarios.**
  Several of the six bugs this scenario walks are also targeted, in
  part, by dedicated cross-cutting specs elsewhere in the section:
  GROK-19212 also spans a step in `uploading.md`; GROK-19103 also
  spans steps in `projects-ui-smoke.md`, `uploading.md`, and
  `projects-copy-clone.md`; GROK-19403 also spans steps in
  `projects-ui-smoke.md` and `projects-copy-clone.md`; GROK-18345 also
  spans a step in `uploading.md`; github-3550 is anchored only here
  (Step 9's Query rename). GROK-19728 is included via a semantic
  match with `lifecycle-api.md`. This scenario is the one place all
  six converge in a single run.
- **"Second user" placeholder.** The literal string "second user" is
  used as the recipient identity in step 12 and as the re-auth
  target identity in step 13. Eventual binding to a specific test
  account is intentionally left open — a decision for the automation
  stage, not a gap in this scenario.
- **Re-auth pattern dependency.** Step 13 needs a
  logout + login-as-different-user + verify sequence. The UI path
  for this is documented in
  `.claude/skills/grok-browser/references/projects.md`; a JS API
  path is still TODO. Until a reusable login-as-other-user helper is
  registered, Step 13 is deferred in
  `complex-share-second-user-spec.ts`.
- **UI coverage owned here.** This scenario is the chain's witness
  for: the rename context-menu surface (Project / Query / Script,
  Step 9), the share dialog with a permissions editor at two access
  levels (Step 12), the re-auth/login-as-second-user flow (Step 13),
  and data-sync refresh verification (Step 11). No other scenario in
  the section covers these flows.
- **Cleanup responsibility.** This scenario produces multiple saved
  projects (from steps 2, 5, 6, 8) and leaves entities renamed and
  possibly moved. `projects-ui-smoke.md` (the section's terminal cleanup
  scenario) is responsible for deleting the produced projects, but
  NOT for undoing the renames/moves — that's a known gap, flagged for
  the automation stage to handle via after-test teardown or an
  isolated environment.
- **Companion files.** `complex-run.md` is a prior manual test-run
  record with no automated-spec counterpart. `complex-ui.md` carries
  the UI-only Step 4 (drag-drop add) and Step 10 (move) content that
  can't be automated.
