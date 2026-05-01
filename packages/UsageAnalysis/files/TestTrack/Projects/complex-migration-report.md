# Migration Report — complex.md

## Step mapping

The original is 13 numbered steps (Markdown source uses `1.` for all
items; auto-numbered as 1-13 on render) + an "**Expected result:**"
section with 3 bullets + trailing JSON `{ "order": 8 }`. The migrated
body preserves all 13 steps with axis-clarity expansion of
sub-bullets, promotes the "Expected result" bullets to an explicit
"### Expected results" sub-section under `## Scenarios`, and accounts
for the trailing JSON in this table.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Open tables from different sources" + 8 sub-bullets (file share / Space / DB / Query / Script / Pivot+Aggregate / Join / Clone) | Scenarios > step 1 (with all 8 sub-bullets preserved as bulleted list, each clarified with example sources) | preserved-with-clarification. The original's bare sub-bullets are augmented with example sources (e.g. "System:Demo > SPGI_v2.csv" for file share) drawn from `uploading.md`'s known sources for cross-scenario consistency. |
| 2. "Save all opened tables as a project with Data Sync enabled" | Scenarios > step 2 (with explicit OK + verify-on-server + cancel-share-dialog) | preserved-with-axis-clarity-expansion (the implicit "wait for save", "cancel auto-share dialog" actions are made explicit per D-STEP-02). |
| 3. "Open tables from: Spaces / Files / Query result / DB table" | Scenarios > step 3 (with 4 sub-bullets preserved) | preserved. |
| 4. "Add newly opened tables to the opened project (drag'n'drop them in the Dashboards)" | Scenarios > step 4 (with explicit drag-drop UI mention + Link/Clone/Move/Copy choice + Context Panel verification) | preserved-with-axis-clarity-expansion. |
| 5. "Open the project and save a copy without Data Sync for all tables" | Scenarios > step 5 | preserved-with-axis-clarity-expansion. |
| 6. "Open the copied project and save it again with Data Sync enabled (result: two projects, both with Data Sync)" | Scenarios > step 6 (with parenthetical preserved as italicized inline note) | preserved. The parenthetical "result: two projects, both with Data Sync" is preserved verbatim because it documents the resulting state after the re-save (the Sync-OFF copy from step 5 is now Sync-ON, plus the original from step 2 was Sync-ON). |
| 7. "Open the last saved project and rename all tables inside it" | Scenarios > step 7 (with rename mechanism options: tab context menu OR JS API) | preserved-with-clarification. |
| 8. "Save a copy of the project with Data Sync enabled" | Scenarios > step 8 (renamed-tables copy, Data Sync ON) | preserved. |
| 9. "Rename the following entities: Project / Query / Script" | Scenarios > step 9 (with 3 sub-bullets preserved + rename UI options for each) | preserved-with-clarification. |
| 10. "Move the following entities to any file share, then to any Space: Script / Query / Project" | Scenarios > step 10 (with 3 sub-bullets preserved + each-entity-moved-twice clarification) | preserved-with-clarification. The original "to file share, then to any Space" is per-entity (each moved twice) — clarified as "each entity is moved twice". |
| 11. "Open the moved project and previously saved project and verify: All tables / relationships / Data Sync" | Scenarios > step 11 (with 3 sub-bullets verifications preserved) | preserved. |
| 12. "Configure project sharing for another user with View and Use access / Full access" | Scenarios > step 12 (with right-click Share dialog mechanism + 2 access-level sub-bullets + ambiguity note for "both levels OR one each") | preserved-with-axis-clarity-expansion. The original's "another user" is bound to the literal "second user" placeholder per `mig-2026-04-29-fixture-placeholder` decision-log entry. |
| 13. "Log in as the second user and open the shared project" | Scenarios > step 13 (with logout + login-as-second-user + navigation + open) | preserved-with-axis-clarity-expansion (logout step is implicit in original; explicit in migrated for D-STEP-02). |
| "**Expected result:**" + 3 bullets (entities function after rename/move; Data Sync updates; users see per permissions) | Scenarios > "### Expected results" sub-section with 3 bullets preserved verbatim | preserved as explicit verification block at the end of `## Scenarios`. D-STEP-02 axis-clarity. |
| Trailing JSON `{ "order": 8 }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/projects.yaml` rev 2 `order_from_files`) |

No original numbered step or expected-result bullet is silently
dropped.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.complex.md.target_layer = playwright`. The scenario
  exercises full UI flows (drag-drop in Dashboards, right-click
  context menus, Share dialog, multi-session logout-login), which
  rule out `api-contract`. `manual-only` would be a degradation —
  the scenario is automatable in principle, only the test-length
  concern (likely SCOPE_REDUCTION at Critic Gate A) suggests it
  may be split into smaller specs at Automator stage.
- **Why this `priority`:** chose `regression` per A-STRUCT-MECH-06
  enum. The scenario is regression-class: it covers GROK-19212 +
  github-3550 (rename + datasync), GROK-19103 (join + save), and
  GROK-18345 (Spaces datasync share) — all regression-prone bugs
  with explicit reproduction paths. Per-cycle Invariant 1 honored.
- **Why this `strategy`:** `simple` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.complex.md.strategy = simple`. Despite the 13-step
  length, the scenario is a single linear narrative (no matrix
  axes, no chained tests over multiple fixtures, no cross-fixture
  composite). The Migrator schema's `simple` enum value matches:
  "single standalone scenario, no matrix axes, no chained tests".
- **D-MERIT-01 compliance: NO opt-outs at Migrator stage.** The
  chain notes flag this scenario as a "Likely SCOPE_REDUCTION
  candidate at Critic Gate A — too long for a single test, may
  need split". However, "too long for one test" is an EFFORT axis,
  NOT a technical dependency. D-MERIT-01 requires opt-outs to cite
  real technical dependencies, never effort. Therefore: Migrator
  preserves all 13 steps faithfully + the Expected results block,
  and PROPOSES (Section 1 of the per-scenario REPORT) a split
  along the rename (steps 7-9) / move (step 10) / share-and-
  second-user (steps 12-13) boundaries as a Critic-Gate-A
  candidate. The Test Designer / Automator stage at downstream
  Critic Gate A may apply the split with full SCOPE_REDUCTION
  justification at THAT stage; not at this one.
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - No `complex-spec.ts` exists. Spec generation is Step 7.
  - `complex-run.md` exists from a prior MCP reproduction run
    (read-only inspection); informs which sub-flows are tractable
    in MCP (drag-and-drop in Dashboards is historically hard for
    automation; tree-drag-drop-actions atlas critical_path
    references this).
  - Adjacent specs (uploading-spec.ts, opening-spec.ts) confirm
    the section's `loginToDatagrok` / `softStep` / `evalJs` /
    `closeAll` convention; read-only inspection.
- **Helpers consulted / candidates:**
  - Reuse from prior scenarios:
    - `helpers.playwright.projects.shareProjectViaContextMenu`
      (from scenario 3) for step 12
    - `helpers.playwright.projects.saveAndReopen` (from scenario 2)
      for steps 2, 5, 6, 8
    - `helpers.playwright.projects.saveCopyWithMode` (from
      scenario 6) for step 5 (Save Copy without Sync) and step 6
      (re-save with Sync)
  - **NEW candidates** specific to this scenario:
    - `helpers.playwright.projects.dragDropOntoDashboard(page,
      sourceTable, targetProject, action: 'link'|'clone'|'move'|
      'copy')` — for step 4's drag-drop into Dashboards.
    - `helpers.playwright.session.logoutAndLoginAs(page,
      credentials)` — for step 13's re-auth. **CRITICAL**: this is
      currently UNDOCUMENTED in
      `grok-browser/references/projects.md` (only a stub comment
      exists at line 138). Flagged as a B14 reference-file-
      addition candidate in Section 3 of the per-scenario REPORT.
- **Bug library consulted:** yes — `bug-library/projects.yaml`
  rev 2. Four bugs intersect this scenario's flows:
  - **GROK-19212** — "Projects fail to open with 'Could not
    resolve table' after a referenced table is renamed". Steps
    7-9 rename tables INSIDE the project + Project/Query/Script
    entities; step 11 reopens and verifies. Direct intersection.
  - **GROK-19103** — "Join result silently saved as a separate
    project that later fails to open". Step 1 sub-bullet "Join
    two tables" + step 2 (Save with Sync) — direct intersection
    of the join + save code path.
  - **GROK-18345** — "Recipient cannot open shared project that
    uses a Spaces dataset saved with data sync". Step 1 sub-
    bullet "Table from Space" + step 2 (Save with Sync) + step 12
    (Share with second user) + step 13 (open as second user) —
    direct intersection of all three components of the bug's
    reproduction.
  - **github-3550** — "External-entity rename invalidation (sister
    of GROK-19212, but for queries instead of tables)". Step 9
    sub-bullet "Query — rename" + step 11 (open and verify) —
    direct intersection.
  - GROK-19403 (share with un-shared deps), GROK-19728 (view-and-
    use failure state), GROK-19750 (save-copy with link drops
    viewers): NOT in `related_bugs` because:
    - GROK-19403 requires the project to depend on un-shared
      scripts/queries; this scenario shares the project with
      sufficient dep co-sharing assumed (or fails if not — but
      the bug's specific reproduction is not exercised).
    - GROK-19728 requires INDUCED failure state on the project
      (broken creation script); this scenario does not induce
      failure.
    - GROK-19750 requires Save-Copy-with-Link mode specifically;
      this scenario's step 5 is Save Copy WITHOUT specifying Link
      mode (just toggling Sync OFF). Not the bug's repro.
- **Decision log queried:** yes — `decision-log.yaml` rev 10 read.
  Key cross-references:
  - `mig-2026-04-29-fixture-placeholder` directly applies — "second
    user" placeholder used in steps 12-13 per Olena's decision.
    Not re-asked.
  - `atlas-2026-04-30-add-projects-shell-share-via-context-menu` +
    `atlas-2026-04-30-cascade-share-via-context-menu` — Invariant
    3 applies; sub_feature INCLUDED in `sub_features_covered`.
  - No other prior decision contradicts this scenario's content.
- **Per-cycle override invariants (all 3) status:**
  - Invariant 1 (priority enum source-of-truth): honored —
    `priority: regression` is canonical per A-STRUCT-MECH-06.
  - Invariant 2 (existing -spec.ts / -api.ts READ-ONLY):
    SATISFIED — no `complex-spec.ts` exists; trivially satisfied.
    Sibling specs read for convention only and not modified.
  - Invariant 3 (atlas-aware sub_features_covered for share):
    APPLIES — step 12 exercises right-click Share dialog.
    `projects.shell.share-via-context-menu` INCLUDED in
    `sub_features_covered`.

## Opt-outs (SCOPE_REDUCTION proposals)

(none) at Migrator stage. Per D-MERIT-01, opt-outs require a real
technical dependency. The "too long for one test" concern is effort,
not a technical dependency, so it does NOT justify a Migrator-stage
SCOPE_REDUCTION. Migrator preserves all 13 steps + expected results
faithfully.

A SPLIT proposal IS surfaced for the downstream Critic Gate A (Test
Designer / Automator stage), where the test-length axis is in scope:

- **Proposed split candidate** (Section 1 of per-scenario REPORT):
  - **Split A: setup-and-save** (steps 1-6) — open many sources, save
    with Sync, augment, save copy without Sync, re-save with Sync.
  - **Split B: rename-and-move** (steps 7-11) — rename tables /
    entities, move entities, verify reopen. Atlas critical_path
    `rename-dependent-entity-reopen` (p1) covers this; the split B
    spec could BE that critical_path's regression test.
  - **Split C: share-and-second-user** (steps 12-13 + Expected
    result bullets 3) — share at two access levels, log in as second
    user, verify access. Atlas critical_paths
    `share-spaces-datasync` (p1, GROK-18345) and
    `view-and-use-failure-state` (p2, GROK-19728) partially cover
    this; the split C spec could combine both regression tests.

Critic Gate A may apply this split or accept the single-spec form
based on Automator's spec-time considerations.

## Deferred items (NOT opt-outs)

- **"Second user" identity placeholder.** Per
  `mig-2026-04-29-fixture-placeholder` — eventual binding to a
  specific test account is TBD. Automator stage owns spec-time
  decision (per-run timestamp account, pre-created pool, fixed
  test-account pair). Real fixture dependency, not effort.
- **Re-auth pattern dependency.**
  `grok-browser/references/projects.md` line 138 has only a stub
  comment. Automator's spec-time logout + login-as-different-user
  + (optionally) login-back implementation depends on the pattern
  being documented. Surfaced as B14 reference-file-addition
  candidate in Section 3 of per-scenario REPORT.
- **Environment dependencies** (Postgres, Spaces, file share with
  CSVs, registered Query, Script). Real prerequisites; Automator
  decides per-environment availability or feature-gating.
- **Drag-and-drop in Dashboards (step 4) MCP feasibility.** Per
  `tree-drag-drop-actions` atlas critical_path notes, drag-and-drop
  via DOM events is historically hard to automate reliably. The
  `complex-run.md` companion (read-only inspection) may indicate
  whether prior MCP runs succeeded on this step. Automator may
  fall back to a JS-API-equivalent (e.g.
  `grok.dapi.projects.addRelation(...)`) at spec time IF the
  drag-and-drop UI proves too flaky.
- **Cleanup of renamed/moved entities + multiple saved projects.**
  `deleting.md` (scenario 9, must_run_last) cleans up produced
  projects but does NOT undo entity renames or moves. Automator's
  `afterAll` may need to restore entity names / locations, OR rely
  on isolated test-environment teardown. Real chain dependency for
  the projects-cleanup half; environmental for the rename/move
  undo half.
- **Step 12 ambiguity: "another user" — single user with both
  access levels, OR two users at different levels?** The original
  is silent. Migrated body preserves the ambiguity ("both levels
  simultaneously to one user OR to two different second users").
  Automator decides at spec time.

## Edge cases

The original lists no explicit edge cases beyond the Expected result
bullets. Implicit edge cases derivable:

- **GROK-19212 regression in step 11.** The rename of tables in
  step 7 + entity rename in step 9 + reopen in step 11 = exactly
  the bug's reproduction sequence. If step 11 fails to find a
  table by its (post-rename) name, GROK-19212 has regressed.
- **github-3550 regression in step 11 (queries axis).** Same as
  above but for the renamed Query in step 9.
- **GROK-19103 regression in step 2 / step 6 / step 8.** The Join
  table from step 1 is included in the Save in step 2 (and copies
  in 5, 6, 8). If on reopen the join result silently splits into a
  separate project that fails to open, GROK-19103 has regressed.
- **GROK-18345 regression in step 13.** The second user opens a
  shared project that includes a Spaces-sourced table from step 1.
  If the second user gets a permission failure or null on the
  Spaces table, GROK-18345 has regressed.
- **Drag-and-drop relation type semantics in step 4.** The drag-
  drop offers Link / Clone / Move / Copy choices. The migrated
  body picks Link as default; behavior of Move and Copy on the
  drag-drop is implicit (Move would relocate the source table to
  the project's namespace; Copy would create independent copies).
  Implicit; verifications preserved.
- **Project rename in step 9 — accessibility post-rename.** The
  project from step 8 is renamed in step 9 sub-bullet 1. Step 10
  refers to "the moved project" (the renamed-and-moved one);
  step 11 expects to "open the moved project". The renamed name
  must be the open-target, not the original. Implicit; preserved.

(none additional)

## Unresolved ambiguities

- **"Another user" / "second user" — single identity vs. two?**
  Step 12 grants two access levels (View-and-Use AND Full); step 13
  logs in as "the second user". Original is silent on whether one
  user has both levels (most likely interpretation: the second user
  has Full access, and View-and-Use is the lower-bound; or grants
  are sequential testing). Migrated body preserves the ambiguity
  with an explicit ambiguity note.
- **Order of access-level grants in step 12.** "View and Use access /
  Full access" — apply both? In what order? Migrated body says
  "view-and-use level (no edit permission)" then "additionally OR
  alternatively grant Full". Automator decides spec-time order.
- **"Table from a Script" semantics in step 1.** A Script that
  produces a table on execution — but what does "Clone a table
  opened from a script" (step 1 last sub-bullet) mean? Clone the
  script-output table itself? Or clone-as-relation in the project?
  Migrated body uses "right-click > Clone (or equivalent UI clone
  action) to produce an independent copy". Automator verifies UI
  behavior at spec time.
- **"Pivot Table > Add" and "Aggregate Rows > Add" UI controls.**
  Same axis as scenario 2's `uploading.md` Cases 8-9 — original
  prose preserved verbatim; Automator cross-references
  `grok-browser/references/` for selectors at spec time.
- **"Drag'n'drop in the Dashboards" mechanism.** Step 4's drag-and-
  drop target — is it the Dashboards card view tile, OR the
  project tree node in the Browse panel, OR the open ProjectView
  in the workspace? Most natural interpretation: drag onto the
  Dashboards card. Automator decides per UI inspection.
- **"Move ... to any file share, then to any Space" in step 10.**
  "any" implies the test-author can choose; for reproducibility
  the Automator should pick a deterministic file share + Space
  pair. Real spec-time decision.
- **Rename UI access points in step 9.** The original lists three
  entities to rename (Project, Query, Script). The rename UI for
  each may differ (Browse > Dashboards > right-click for Project;
  Queries pane right-click for Query; Scripts manage UI for
  Script). Migrated body documents the most likely access points;
  Automator verifies each at spec time.
- **"Verify Data Sync updates correctly" in step 11.** The
  verification mechanism — force a refresh? wait for an
  auto-update interval? introduce a server-side data change and
  observe? — is silent in the original. Migrated body says "force
  a Data Sync refresh OR re-open after a server-side data change".
  Automator decides spec-time mechanism.
