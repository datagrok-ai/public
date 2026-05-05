# Migration Report — complex.md (re-migration rev 3)

This is a **re-migration** of `complex.md` per chain YAML revision 3
(`scenario-chains/projects.yaml` rev 3 — 2026-05-04). It overwrites
the prior 2026-04-30 migration (decision-log
`mig-2026-04-30-complex-migration`) and applies the rev-3 schema
additions:

- `pyramid_layer: bug-focused` (Rule 3 — 6 bugs walked).
- `ui_coverage_responsibility:` 7 flows (rename × 3, share × 2, re-
  auth × 1, data-sync × 1).
- `ui_coverage_delegated_to: null`.
- `related_bugs:` expanded to 6 — adds GROK-19403 and GROK-19728 to
  the prior 4 (GROK-19212, GROK-19103, GROK-18345, github-3550) per
  chain rev 3 `bug_focused_candidates[]` spans referencing
  `complex.md:Step 12`.
- Split-spec preservation: Wave 1b 3-spec realization
  (`complex-derived-tables-spec.ts`, `complex-rename-spec.ts`,
  `complex-share-second-user-spec.ts`) honored per decision-log
  `wave-1b-complex-split-b70-followup` (2026-05-01) and reinforced
  by `b2-2026-05-03-drag-drop-ui-only-reclassification` (move sub-
  spec dropped, NOT to be re-merged).

## Step mapping

The original is 13 numbered steps + an "**Expected result:**" section
with 3 bullets + trailing JSON `{ "order": 8 }`. The migrated body
preserves all 13 steps with axis-clarity expansion of sub-bullets,
promotes the Expected result bullets to an explicit
`### Expected results` sub-section under `## Scenarios`, and accounts
for the trailing JSON in this table. Steps 4 and 10 are flagged as
UI-only and content moved to `complex-ui.md` per
`b2-2026-05-03-drag-drop-ui-only-reclassification` (preserving the
original step numbers in body for cross-reference).

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Open tables from different sources" + 8 sub-bullets (file share / Space / DB / Query / Script / Pivot+Aggregate / Join / Clone) | Scenarios > step 1 (with all 8 sub-bullets preserved as bulleted list, each clarified with example sources) | preserved (split for clarity) |
| 2. "Save all opened tables as a project with Data Sync enabled" | Scenarios > step 2 (with explicit OK + verify-on-server + cancel-share-dialog) | preserved as verification |
| 3. "Open tables from: Spaces / Files / Query result / DB table" | Scenarios > step 3 (with 4 sub-bullets preserved) | preserved |
| 4. "Add newly opened tables to the opened project (drag'n'drop them in the Dashboards)" | Scenarios > step 4 reference + body content moved to complex-ui.md | preserved (moved to atlas/ui-companion — see complex-ui.md; UI-only per b2-2026-05-03-drag-drop-ui-only-reclassification) |
| 5. "Open the project and save a copy without Data Sync for all tables" | Scenarios > step 5 | preserved |
| 6. "Open the copied project and save it again with Data Sync enabled (result: two projects, both with Data Sync)" | Scenarios > step 6 (with parenthetical preserved as italicized inline note) | preserved |
| 7. "Open the last saved project and rename all tables inside it" | Scenarios > step 7 (with rename mechanism options: tab context menu OR JS API) | preserved (split for clarity) |
| 8. "Save a copy of the project with Data Sync enabled" | Scenarios > step 8 (renamed-tables copy, Data Sync ON) | preserved |
| 9. "Rename the following entities: Project / Query / Script" | Scenarios > step 9 (with 3 sub-bullets preserved + rename UI options for each) | preserved (split for clarity) |
| 10. "Move the following entities to any file share, then to any Space: Script / Query / Project" | Scenarios > step 10 reference + body content moved to complex-ui.md | preserved (moved to atlas/ui-companion — see complex-ui.md; UI-only per b2-2026-05-03-drag-drop-ui-only-reclassification) |
| 11. "Open the moved project and previously saved project and verify: All tables / relationships / Data Sync" | Scenarios > step 11 (with 3 sub-bullets verifications preserved) | preserved as verification |
| 12. "Configure project sharing for another user with View and Use access / Full access" | Scenarios > step 12 (with right-click Share dialog mechanism + 2 access-level sub-bullets + ambiguity note for "both levels OR one each") | preserved (split for clarity) |
| 13. "Log in as the second user and open the shared project" | Scenarios > step 13 (with logout + login-as-second-user + navigation + open) | preserved as verification |
| "**Expected result:**" + 3 bullets (entities function after rename/move; Data Sync updates; users see per permissions) | Scenarios > "### Expected results" sub-section with 3 bullets preserved verbatim | preserved as verification |
| Trailing JSON `{ "order": 8 }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/projects.yaml` rev 3 `order_from_files`) |

No original numbered step or expected-result bullet is silently
dropped.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 3
  `output_plan.complex.md.target_layer = playwright`. The scenario
  exercises full UI flows (right-click context menus, Share dialog
  with PermissionsEditor, multi-session logout-login), which rule
  out `api-contract`. Drag-drop sub-flows (Steps 4 and 10) moved
  to `complex-ui.md` per
  `b2-2026-05-03-drag-drop-ui-only-reclassification`.
- **Why this `coverage_type`:** chose `regression` per A-STRUCT-MECH-06
  enum. The scenario is regression-class: it covers GROK-19212 +
  github-3550 (rename + datasync), GROK-19103 (join + save),
  GROK-18345 (Spaces datasync share), GROK-19403 (share with un-
  shared deps), and GROK-19728 (view-and-use under failure state)
  — all regression-prone bugs with explicit reproduction paths.
- **Why this `pyramid_layer`:** chose `bug-focused` per chain rev 3
  Rule 3 (related_bugs frontmatter has 6 entries; multiple Steps
  walk distinct bug repros — Steps 7-9 + 11 → GROK-19212; Step 1 +
  Step 2 → GROK-19103; Step 1 + Step 2 + Step 12 + Step 13 →
  GROK-18345 [non-adjacent, Trigger 2 motivating]; Step 9 → github-
  3550; Step 12 + Step 13 → GROK-19403; Step 12 → GROK-19728).
  Discriminator passes for all 6: each bug fails this scenario before
  fix. Bug-focused dominates over Rule 4 multi-source per chain
  rev 3 priority 3 → 4.
- **Why this `strategy`:** `split (Wave 1b)` per chain rev 3
  `output_plan.complex.md.strategy`, with `realized_as` listing the
  3 satellite specs. Single-spec coverage of all 13 steps was
  rejected at Critic E re-run as too heavy; the 3-spec split was
  executed at Wave 1b (decision-log
  `wave-1b-complex-split-b70-followup`) and is preserved here. The
  `complex-move` sub-spec was dropped per
  `b2-2026-05-03-drag-drop-ui-only-reclassification` and is NOT
  to be re-merged.
- **UI coverage delegation status:** `ui_coverage_delegated_to: null`.
  All 7 UI flows in `ui_coverage_responsibility` are owned here, not
  delegated. UI flows: `context-menu-rename-project`,
  `context-menu-rename-query`, `context-menu-rename-script`
  (Step 9); `pcmdShareProject`, `share-dialog-permissions-editor`
  (Step 12); `logout-login-as-second-user` (Step 13);
  `data-sync-refresh-verification` (Step 11).
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `complex-derived-tables-spec.ts` (Wave 1b — GROK-19103) — read-
    only inspection; honors split decision, do NOT re-merge.
  - `complex-rename-spec.ts` (Wave 1b — GROK-19212; partial github-
    3550) — read-only inspection.
  - `complex-share-second-user-spec.ts` (Wave 1b — GROK-18345 share-
    side; Step 13 deferred) — read-only inspection.
  - `complex-run.md` (prior MCP reproduction run) — read-only
    inspection.
  - `complex-ui.md` (UI-companion for Steps 4 and 10) — read-only
    inspection per `b2-2026-05-03-drag-drop-ui-only-reclassification`.
  - Adjacent specs (uploading-spec.ts, opening-spec.ts) — confirmed
    section convention (`loginToDatagrok`, `softStep`, `evalJs`,
    `closeAll`); read-only.
- **Helpers consulted / candidates:**
  - **Reused (registered in helpers-registry.yaml):**
    - `helpers.playwright.projects.shareProjectViaContextMenu`
      (Step 12).
    - `helpers.playwright.projects.saveAndReopen` (Steps 2, 5, 6, 8).
    - `helpers.playwright.projects.saveCopyWithMode` (Steps 5, 6, 8).
  - **Candidate helpers (NOT yet in registry — flagged for
    addition):**
    - `helpers.playwright.projects.dragDropOntoDashboard` (Step 4 —
      flagged in `mig-2026-04-30-complex-migration`; UI-only per
      `b2-2026-05-03-drag-drop-ui-only-reclassification`; helper
      candidate parked pending Playwright drag-drop support).
    - `helpers.playwright.session.logoutAndLoginAs` (Step 13 — Helper
      3 from helpers-batch-1 backlog; blocking spec realization of
      Step 13 in `complex-share-second-user-spec.ts`; surfaced by
      `wave-1b-complex-split-b70-followup`).
    - `helpers.playwright.projects.move` (Step 10 — surfaced by
      Wave 1b complex-move sub-spec drop analysis;
      `complex-move-spec.ts` deferred pending registration).
    - `helpers.playwright.projects.rename` (Step 9 Project rename —
      surfaced by Wave 1b complex-rename round-1 hypothesis; JS API
      `p.name = '<new>'; dapi.save(p)` failed to propagate on dev;
      either UI helper or platform investigation needed).
- **Bug library consulted:** yes — `bug-library/projects.yaml`
  rev 2. Six bugs intersect this scenario's flows (per chain rev 3
  `bug_focused_candidates[]` spans referencing complex.md):
  - **GROK-19212** — Steps 7-9 (rename) + Step 11 (reopen-verify);
    cross-cutting candidate spec `projects-grok-19212-spec.ts`
    (spans uploading.md:Step 7 + complex.md:Step 7).
  - **GROK-19103** — Step 1 sub-bullet "Join two tables" + Step 2
    (Save with Sync); cross-cutting candidate spec
    `projects-grok-19103-spec.ts` (4 affecting scenarios — upload-
    project, uploading, projects-copy-clone, complex).
  - **GROK-19403** — Step 12 (Share at View-and-Use) + Step 13
    (recipient open) intersects share-with-un-shared-deps repro;
    cross-cutting candidate spec `projects-grok-19403-spec.ts`
    (spans share-project.md:Step 4 + complex.md:Step 12 + projects-
    copy-clone.md:Step 5).
  - **GROK-18345** — Step 1 (Spaces) + Step 2 (Sync ON) + Step 12
    (Share) + Step 13 (recipient open); **single-scenario non-
    adjacent (Trigger 2 motivating example)**; cross-cutting
    candidate spec `projects-grok-18345-spec.ts` (spans uploading.md:
    Step 7 + complex.md:Step 12).
  - **GROK-19728** — Step 12 (View-and-Use grant) + Step 13
    (recipient access under failure state) intersects view-and-use
    failure-state repro; cross-cutting candidate spec
    `projects-grok-19728-spec.ts` (chain rev 3 emits under Trigger 1
    with lifecycle-api.md inclusion via semantic step match;
    flagged for review per chain rev 3 rationale).
  - **github-3550** — Step 9 sub-bullet "Query/Script — rename" +
    Step 11 (open and verify); cross-cutting candidate spec
    `projects-github-3550-spec.ts` (chain rev 3 spans complex.md:
    Step 9; emits under Trigger 1 cross-scenario).
  - GROK-19750 (Save-Copy-with-Link drops viewers): NOT in
    `related_bugs` — Step 5 toggles Sync OFF on Save Copy, not
    "Link" mode specifically; not the bug's repro. Owned by
    `projects-copy-clone.md` and `upload-project.md` per chain rev 3.
- **Decision log queried:** yes — `decision-log.yaml` read; key
  cross-references applied (D10 active query):
  - `mig-2026-04-29-fixture-placeholder` — "second user" placeholder
    preserved as literal string in steps 12-13.
  - `mig-2026-04-30-complex-migration` — prior 2026-04-30 migration;
    this re-migration overwrites it per chain rev 3 schema.
  - `b14-2026-04-30-re-auth-pattern-applied` — re-auth pattern
    documented in grok-browser/references/projects.md (UI path; JS
    API TODO); Helper 3 candidate flagged.
  - `wave-1b-complex-split-b70-followup` — 3-spec split realized;
    move sub-spec dropped; preserved here per task constraint
    "split decision honored (don't re-merge into single spec)".
  - `b2-2026-05-03-drag-drop-ui-only-reclassification` — Steps 4
    and 10 reclassified UI-only; body moved to complex-ui.md;
    preserved here.
  - `atlas-2026-04-30-add-projects-shell-share-via-context-menu` +
    `atlas-2026-04-30-cascade-share-via-context-menu` — Invariant 3
    applies; sub_feature INCLUDED in `sub_features_covered`.
- **Cross-cutting bug citations (chain rev 3 bug_focused_candidates).**
  Per `migration-prompt.md` "Cross-cutting bug citations from chain
  YAML": all 6 bugs whose `spans` reference complex.md are cited
  in Notes section of migrated body — see "Cross-cutting candidate
  specs" note. RECOMMENDED, not mandatory; F-BUG-COVERAGE-01 at
  section-complete is authoritative.
- **Per-cycle override invariants status:**
  - Invariant 1 (coverage_type enum source-of-truth): honored —
    `coverage_type: regression` is canonical.
  - Invariant 2 (existing -spec.ts / -api.ts READ-ONLY): SATISFIED —
    3 satellite specs (`complex-derived-tables-spec.ts`,
    `complex-rename-spec.ts`, `complex-share-second-user-spec.ts`)
    inspected read-only; not modified per task constraint.
  - Invariant 3 (atlas-aware sub_features_covered for share):
    APPLIES — step 12 exercises right-click Share dialog.
    `projects.shell.share-via-context-menu` INCLUDED.

## Opt-outs (SCOPE_REDUCTION proposals)

- **`complex-move` sub-spec dropped** (covers Step 10 — move entities
  to file share / Space). Technical dependency:
  `b2-2026-05-03-drag-drop-ui-only-reclassification` — drag-drop event
  registration not automatable through current Playwright + JS API
  mechanisms (3 mechanisms tried: synthetic DragEvent, raw mouse
  events, CDP Input.dispatchMouseEvent); right-click `Move to` menu
  option does not exist in current UI (verified 2026-05-03 against
  dev.datagrok.ai). Reduce-scope path chosen over inlining fragile
  selector code. Tracked for `helpers.playwright.projects.move` helper
  registration follow-up. UI coverage of the underlying flow lives
  in `complex-ui.md` (Step 10 body preserved there).
- **Step 13 (recipient-side open) deferred in
  `complex-share-second-user-spec.ts`.** Technical dependency:
  `helpers.playwright.session.logoutAndLoginAs` (Helper 3 from
  helpers-batch-1 backlog) NOT yet registered. The Step 12 share-
  side grant via JS API IS exercised (when env permits — currently
  skips with FK violation note on dev). Recipient-side reopen
  verification (the actual GROK-18345 + GROK-19403 + GROK-19728
  recipient-side bug reproduction) requires the re-auth helper.
- **Step 4 drag-drop content moved to `complex-ui.md`.** Technical
  dependency: same as `complex-move` above —
  `b2-2026-05-03-drag-drop-ui-only-reclassification`. UI coverage
  lives in `complex-ui.md` (Step 4 body preserved there). Section
  excluded from `-spec.ts` generation.

## Deferred items (NOT opt-outs)

- **"Second user" identity placeholder.** Per
  `mig-2026-04-29-fixture-placeholder` — eventual binding to a
  specific test account is TBD. Automator stage owns spec-time
  decision (per-run timestamp account, pre-created pool, fixed
  test-account pair). Real fixture dependency, not effort.
- **Re-auth pattern JS API path TODO.** Per
  `b14-2026-04-30-re-auth-pattern-applied`, the UI re-auth path is
  documented in `grok-browser/references/projects.md`; JS API path is
  flagged TODO awaiting external confirmation of whether
  `grok.dapi` exposes a session re-auth API. Real reference-file
  dependency.
- **Environment dependencies** (Postgres, Spaces, file share with
  CSVs, registered Query, Script). Real prerequisites; Automator
  decides per-environment availability or feature-gating.
- **Step 9 Project-entity rename JS API investigation.** Per
  `wave-1b-complex-split-b70-followup` Wave 1b round-1 finding,
  `p.name = '<new>'; await dapi.save(p)` failed to propagate on dev
  (renamed project not findable via dapi.filter at new name). Either
  documented JS API approach is broken (file as bug) or alternative
  approach needed. Real platform-investigation prerequisite for the
  Step 9 Project-entity rename surface in
  `complex-rename-spec.ts`.
- **Cleanup of renamed/moved entities + multiple saved projects.**
  `deleting.md` (must_run_last) cleans up produced projects but does
  NOT undo entity renames or moves. Automator's `afterAll` may need
  to restore entity names / locations, OR rely on isolated test-
  environment teardown. Real chain dependency for the projects-
  cleanup half; environmental for the rename/move undo half.
- **Step 12 ambiguity: "another user" — single user with both
  access levels, OR two users at different levels?** The original
  is silent. Migrated body preserves the ambiguity ("both levels
  simultaneously to one user OR to two different second users").
  Real source-text ambiguity; Automator decides at spec time.
- **Step 13 second user account provisioning.** Auto-creating the
  second user account at test-time depends on either a fixture
  builder (helpers candidate) OR a pre-provisioned test account
  pool. Real fixture dependency.

## Edge cases

The original lists no explicit edge cases beyond the Expected result
bullets. Implicit edge cases derivable:

- **GROK-19212 regression in step 11.** The rename of tables in
  step 7 + entity rename in step 9 + reopen in step 11 = exactly
  the bug's reproduction sequence. If step 11 fails to find a
  table by its (post-rename) name, GROK-19212 has regressed.
  PRESERVED as scenario step (Step 11 verification).
- **github-3550 regression in step 11 (queries axis).** Same as
  above but for the renamed Query in step 9. PRESERVED as scenario
  step (Step 11 verification); query/script axis partial coverage
  in `complex-rename-spec.ts` per Wave 1b outcome.
- **GROK-19103 regression in step 2 / step 6 / step 8.** The Join
  table from step 1 is included in the Save in step 2 (and copies
  in 5, 6, 8). If on reopen the join result silently splits into a
  separate project that fails to open, GROK-19103 has regressed.
  PRESERVED as scenario step. Realized in
  `complex-derived-tables-spec.ts`.
- **GROK-18345 regression at Step 13 (Trigger 2 motivating
  example).** The second user opens a shared project that includes
  a Spaces-sourced table from step 1, saved with Sync ON in step 2,
  shared in step 12. If the second user gets a permission failure
  or null on the Spaces table, GROK-18345 has regressed. PRESERVED
  as scenario step (Step 13 verification); recipient-side coverage
  deferred in `complex-share-second-user-spec.ts` pending Helper 3.
- **GROK-19403 regression at Step 13 recipient open.** If the
  second user (View-and-Use access) opens the shared project and
  the underlying Script/Query was not also shared, recipient sees
  silent null instead of explicit failure. PRESERVED as scenario
  step (Step 13 verification); recipient-side coverage deferred.
- **GROK-19728 regression at Step 12 + Step 13 (View-and-Use under
  failure state).** If the project enters failure state (e.g.
  broken creation script via the renamed entity in step 9), the
  View-and-Use recipient must NOT be able to edit the creation
  script. PRESERVED as scenario step (Expected results bullet 3);
  recipient-side coverage deferred.
- **Drag-and-drop relation type semantics in step 4.** The drag-
  drop offers Link / Clone / Move / Copy choices. Body content
  moved to `complex-ui.md` per
  `b2-2026-05-03-drag-drop-ui-only-reclassification`; UI-companion
  preserves the choice surface.
- **Project rename in step 9 — accessibility post-rename.** The
  project from step 8 is renamed in step 9 sub-bullet 1. Step 11
  expects to "open the moved project". The renamed name must be
  the open-target. PRESERVED as scenario step; investigation
  deferred per Wave 1b round-1 finding.

## Unresolved ambiguities

- **"Another user" / "second user" — single identity vs. two?**
  Step 12 grants two access levels (View-and-Use AND Full); step 13
  logs in as "the second user". Original is silent on whether one
  user has both levels (most likely interpretation: one user with
  Full access, View-and-Use as the lower-bound; OR sequential
  test-time switching). Migrated body preserves the ambiguity with
  an explicit ambiguity note.
- **Order of access-level grants in step 12.** "View and Use access
  / Full access" — apply both? In what order? Migrated body says
  "view-and-use level (no edit permission)" then "additionally OR
  alternatively grant Full". Automator decides spec-time order;
  current Wave 1b realization in `complex-share-second-user-spec.ts`
  applies both sequentially via JS API.
- **"Table from a Script" semantics in step 1.** A Script that
  produces a table on execution — what does "Clone a table opened
  from a script" (step 1 last sub-bullet) mean? Clone the script-
  output table itself? Or clone-as-relation in the project?
  Migrated body uses "right-click > Clone (or equivalent UI clone
  action) to produce an independent copy". Automator verifies UI
  behavior at spec time.
- **"Pivot Table > Add" and "Aggregate Rows > Add" UI controls.**
  Same axis as scenario 2's `uploading.md` Cases 8-9 — original
  prose preserved verbatim; Automator cross-references
  `grok-browser/references/` for selectors at spec time.
- **"Move ... to any file share, then to any Space" in step 10.**
  "any" implies the test-author can choose; for reproducibility
  the Automator should pick a deterministic file share + Space
  pair. Body content moved to `complex-ui.md`; spec-time decision
  deferred there.
- **Rename UI access points in step 9.** The original lists three
  entities to rename (Project, Query, Script). The rename UI for
  each may differ. Migrated body documents the most likely access
  points; Automator verifies each at spec time. Wave 1b round-1
  finding: JS API project-rename setter not propagating on dev —
  flagged for platform investigation.
- **"Verify Data Sync updates correctly" in step 11.** The
  verification mechanism — force a refresh? wait for an
  auto-update interval? introduce a server-side data change and
  observe? — is silent in the original. Migrated body says "force
  a Data Sync refresh OR re-open after a server-side data change".
  Automator decides spec-time mechanism. UI flow
  `data-sync-refresh-verification` owned by this scenario per
  chain rev 3.
- **GROK-19728 in chain rev 3 emission rationale.** Chain rev 3
  notes: "stricter intersection yields exactly 1 affecting scenario
  (lifecycle-api.md) which lacks the failure-state repro flow;
  emitting under Trigger 1 with complex.md inclusion via semantic
  step match for cross-cutting invariant of 'view-and-use cannot
  edit creation script under failure state' — flagged for review."
  Migrated frontmatter INCLUDES GROK-19728 in `related_bugs` per
  chain rev 3 task constraint; the semantic-match rationale is
  preserved here for audit.
