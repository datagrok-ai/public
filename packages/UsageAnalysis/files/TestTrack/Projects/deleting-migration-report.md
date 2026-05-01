# Migration Report — deleting.md

## Step mapping

The original is 4 numbered steps + a leading `>Note: make sure this
is the last test case` blockquote + trailing JSON `{ "order": 6 }`.
The migrated body preserves all 4 steps + the must-run-last semantics
(now captured by `must_run_last: true` in chain rev 2 dependency_graph
+ explicit Notes-section reference) and the trailing JSON.

| Original element | Migrated location | Decision |
|------------------|-------------------|----------|
| `>Note: make sure this is the last test case` (leading blockquote) | Notes section (`Must-run-last invariant`) + chain reference | preserved-as-invariant. The note is no longer a literal step but a chain-level invariant captured in `scenario-chains/projects.yaml` rev 2 `dependency_graph.deleting.md.must_run_last: true`. Migrated body cites this in the Notes section explicitly. |
| 1. "Find the projects from the previous steps" | Scenarios > step 1 (Browse + locate project) | preserved-with-resolution. The phrase "projects from the previous steps" is resolved per chain rev 2 `consumes` field to an explicit upstream-produced project enumeration in the Setup section. |
| 2. "Right-click the project and select Delete project from the context menu. A dialog opens. Optionally, delete it using the drop-down menu next to the project name on the Context Panel" | Scenarios > step 2 (with two explicit Option A/B sub-bullets) | preserved-with-axis-clarity-expansion. The original's "Optionally..." alternative is preserved as Option B (Context Panel dropdown), with Option A being the right-click context menu. Both achieve the same delete trigger. |
| 3. "Click DELETE in the confirmation dialog" | Scenarios > step 3 (with selector hint citation) | preserved-with-clarification. Added cross-reference to `grok-browser/references/projects.md` line 153 for the DELETE button's no-`name=` attribute selector convention. |
| 4. "Check that the project has been deleted and is no longer present" | Scenarios > step 4 (Browse verification + JS API double-check) | preserved-with-axis-clarity-expansion. Added explicit JS API verification (`grok.dapi.projects.filter(...).list().length === 0`) as a complementary check to the visual Dashboards-listing verification. |
| Trailing JSON `{ "order": 6 }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/projects.yaml` rev 2 `order_from_files`). Note: the chain's `must_run_last: true` overrides the numeric `order: 6` for chain execution semantics — `order: 6` is purely TestTrack-display ordering. |

No original numbered step or note is silently dropped. The
must-run-last invariant is now captured TWO ways (chain artifact +
migrated body Notes), making it discoverable from either side.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.deleting.md.target_layer = playwright`. The scenario
  exercises Browse > Dashboards UI navigation, right-click context
  menus, confirmation dialog interaction, and post-delete UI
  verification. `api-contract` cannot exercise the UI flow (and
  `grok.dapi.projects.delete()` is a JS API path that bypasses the
  UI delete-confirmation, which is what the test needs to verify).
- **Why this `priority`:** chose `regression` per A-STRUCT-MECH-06
  enum (`smoke | regression | edge | perf`). The scenario is multi-
  target delete verification (~26-30 projects depending on chain
  state) — NOT a single golden-path smoke. Per-cycle Invariant 1
  honored.
- **Why this `strategy`:** `end_to_end_fixtures` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.deleting.md.strategy = end_to_end_fixtures`. The
  scenario consumes the cumulative output of all 8 preceding
  scenarios as its fixture set. This is the canonical
  `end_to_end_fixtures` shape: a per-test consumed fixture
  composed of upstream side-effects.
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - **`deleting-spec.ts` EXISTS** at
    `public/packages/UsageAnalysis/files/TestTrack/Projects/deleting-spec.ts`.
    Read-only inspection (per Invariant 2 of the per-cycle
    override): the existing spec uses a single hardcoded test
    project (similar pattern to `opening-spec.ts`'s 1-project
    reduction) rather than iterating over the full 26-30
    upstream-produced project list. The existing spec's
    reduction is an Automator-stage decision; this migration
    preserves the full upstream-fixture-consumer intent at the
    `.md` level.
  - Adjacent specs (`uploading-spec.ts`, `opening-spec.ts`,
    `browser-spec.ts`) all use the section's
    `loginToDatagrok` / `softStep` / `evalJs` / `closeAll`
    convention; read-only inspection confirms.
- **Helpers consulted / candidates:**
  - Reuse from prior scenarios: `helpers.playwright.projects.saveAndReopen`
    (from scenario 2) is NOT applicable here (no save+reopen).
    `helpers.playwright.projects.openProjectFromDashboards` (from
    scenario 4) is partially applicable (the find-in-Dashboards
    half), but this scenario doesn't open the project — it
    deletes it.
  - **NEW candidates** specific to this scenario:
    - `helpers.playwright.projects.deleteProjectViaContextMenu(page, projectName, mode: 'right-click'|'context-panel-dropdown')`
      — drives steps 1-3 (locate + right-click OR Context Panel
      dropdown + DELETE confirmation). The `mode` parameter
      preserves the original's two-option flexibility.
    - `helpers.playwright.projects.collectChainProducedProjects(page, options: { sinceTimestamp: number, ownedBy: string })`
      — collect the produced-project list at runtime via JS API
      filter for the `beforeAll` block. Replaces a hardcoded
      project list (which would be stale on chain config
      changes).
- **Bug library consulted:** yes —
  `bug-library/projects.yaml` rev 2. **No bug intersects plain
  delete-flow as of rev 2.** The historical github-3752 bug
  (FK violation on delete) was REMOVED from bug-library per
  `decision-log.yaml :: migration_decisions :: mig-2026-04-29-bug-removed`
  ("bad exemplar per Olena 2026-04-29"). The associated atlas
  critical_path `delete-with-view-layouts` was also deleted in
  that cleanup (per
  `mig-2026-04-29-atlas-cascade-completion`). As a result,
  `related_bugs: []` is the correct call — no bug-library entry
  is in scope for this scenario.
- **Decision log queried:** yes —
  `decision-log.yaml` rev 13 read. Cross-references:
  - `mig-2026-04-29-bug-removed` directly applies (github-3752
    removed; this scenario's `related_bugs: []` is the downstream
    consequence). Honored.
  - `mig-2026-04-29-atlas-cascade-completion` directly applies
    (cascade swept github-3752 references from atlas; this
    scenario's coverage doesn't reference the deleted critical_path).
  - No other prior decision contradicts this scenario's content.
- **Per-cycle override invariants (all 3) status:**
  - Invariant 1 (priority enum source-of-truth): honored —
    `priority: regression` is canonical per A-STRUCT-MECH-06.
  - Invariant 2 (existing -spec.ts / -api.ts READ-ONLY):
    SATISFIED — existing `deleting-spec.ts` was read for
    sibling-test convention but NOT modified. Verified via git
    status (only `.md` files modified/untracked; spec files
    untouched).
  - Invariant 3 (atlas-aware sub_features_covered for share):
    correctly does NOT apply — no Share dialog operation in this
    scenario. `projects.shell.share-via-context-menu` correctly
    OMITTED.

## Opt-outs (SCOPE_REDUCTION proposals)

(none) — every original numbered step + the must-run-last note +
the trailing JSON is preserved. The "Optionally..." alternative
in step 2 is preserved as Option B (not a SCOPE_REDUCTION; both
options are viable trigger paths).

## Deferred items (NOT opt-outs)

- **Upstream-produced project list completeness.** The Automator's
  `beforeAll` filter must collect ALL projects produced by chain
  scenarios 1-8, not a hardcoded subset. If the filter
  under-collects, this scenario's terminal cleanup will leave
  leftover projects in the test environment, breaking the chain's
  termination contract. Real prerequisite, not effort. The
  migrated body proposes a runtime-collection approach
  (`grok.dapi.projects.filter('author = currentUser AND createdOn > chainStart').list()`)
  but the exact filter expression depends on Automator's spec
  context.
- **Order of deletion within the chain.** The original is silent
  on whether deletions can happen in parallel OR must be
  sequential. Sequential is safer for verification (each delete
  → verify pair is atomic), but parallel may be needed for
  test-runtime efficiency on large fixture sets. Automator
  decides at spec time.
- **DELETE button selector.** Per `grok-browser/references/projects.md`
  line 153, the DELETE button has no `name=` attribute and must
  be located by visible text. Stable as long as button text doesn't
  change in the platform. Automator may want to add a name=
  attribute upstream (out of Migrator scope). Real test-stability
  prerequisite.
- **JS API double-check vs UI-only verification.** Step 4's added
  JS API double-check (`grok.dapi.projects.filter(...).list()`) is
  a complementary verification — useful when UI listing has
  caching/refresh latency. Automator may use only one of the two
  if the platform's UI consistency is sufficient. Either is valid;
  both is most robust.
- **Cleanup of side-effects beyond projects.** Some chain scenarios
  produced side-effects beyond the project list:
  - `complex.md` step 9 renames Project / Query / Script
  - `complex.md` step 10 moves Script / Query / Project to file
    share / Space
  - `share-project.md` + `complex.md` step 12 auto-create user
    accounts via email-invite
  - `custom-creation-scripts.md` step 4 mutates
    `System:DemoFiles/chem` filesystem
  These are NOT cleaned up by this scenario (which only deletes
  projects). The renamed/moved entities, auto-created users, and
  filesystem mutations require separate cleanup at Automator's
  `afterAll` OR isolated-environment teardown. Surfaced as the
  primary deferred item for chain-completeness.

## Edge cases

The original lists no explicit edge cases. Implicit edge cases
derivable:

- **Project with active sharing relations.** Some chain projects
  (from `share-project.md`, `complex.md` step 12) are shared with
  recipients (registered users + auto-created email users). The
  delete operation MUST cascade or refuse cleanly — the platform's
  delete semantic for shared projects is what's exercised.
- **Project consumed by another open chain step.** If the chain
  is run with parallel scenarios (not the Phase 1 plan, but
  hypothetically), deleting a project that another scenario is
  reading could cause the reader to fail. Sequential chain
  execution avoids this.
- **Project with linked vs. cloned tables.** From
  `projects-copy-clone.md`, the `<original>-link` variant shares
  data with the `<original>`. Deleting the link variant should
  NOT affect the original; deleting the original SHOULD cascade
  to the link variant (or the link variant becomes broken on
  reopen). Implicit; the test's coverage is delete-success
  verification, not link-semantic verification.
- **Project with bound creation script (from
  `custom-creation-scripts.md`).** Deleting the project should
  also handle the bound creation-script reference cleanly (no
  dangling reference). The historical github-3752 bug (now
  removed from bug-library) was about FK-violation in this
  general area; current state is bug-free per rev 2 cleanup.
- **DELETE confirmation dialog cancellation.** Step 3 says "Click
  DELETE in the confirmation dialog" — the original does NOT
  exercise the Cancel path. Implicit edge: cancelling the dialog
  should NOT delete the project. NOT covered in this scenario;
  flag for retro if a separate cancel-test is needed.

(none additional)

## Unresolved ambiguities

- **"Projects from the previous steps" — runtime collection
  vs. enumeration.** Migrated body Setup section provides BOTH
  an enumerated upstream-produced list (from chain rev 2
  `consumes`) AND a runtime-collection JS API filter approach.
  Automator decides which to use; the runtime approach is more
  robust to chain-config changes.
- **Sequential vs parallel delete iteration.** See Deferred
  items.
- **DELETE confirmation dialog timing.** The original implies a
  click-and-wait pattern; Automator may need explicit waits for
  the dialog to render before clicking DELETE, and for the
  Dashboards refresh to reflect the deletion before step 4
  verification.
- **Verify-absence depth.** Step 4 says "no longer present" —
  visual absence in Dashboards listing? Or absence under any
  filter? Or absence in `grok.dapi.projects.list()` as well? The
  migrated body uses dual verification (visual + JS API filter)
  for robustness; Automator may relax.
- **Cleanup of non-project side-effects.** See Deferred items
  (cleanup of side-effects beyond projects).
- **Scope of "the projects" in step 1.** Strict interpretation:
  ALL projects produced by chain. Loose interpretation: only
  projects whose names appear in upstream chain output. The
  migrated body's runtime-collection approach handles both
  interpretations naturally (filter by author + createdOn since
  chain start).
