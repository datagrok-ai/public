# Migration Report — custom-creation-scripts.md

## Step mapping

The original is 5 numbered steps (Markdown source uses `1.` for all
items; auto-numbered as 1-5 on render) + an inline JS code block
embedded in step 1 + bullet-list verifications in step 5. **No
trailing JSON `order` metadata** (unlike sibling scenarios in this
section). The migrated body preserves all 5 steps + the embedded
script + all step-5 verifications.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Run the script:" + inline JS code block | Scenarios > step 1 (with code block preserved + execution-mechanism note added) | preserved-with-axis-clarity-expansion. Original's bare "Run the script" lacks the execution mechanism. Migrated body documents UI options ("Tools > Scripting > paste + Run, OR `grok.functions.eval(...)` programmatically") + adds explicit verification of script output (dataframe produced, table opens). The JS code block is preserved verbatim character-for-character; only the indentation is normalized to 3-space (consistent within the migrated body, vs the source's mixed 0/1/4-space indentation which was a Markdown rendering artifact). |
| 2. "Add some viewers and save the project with data sync enabled" | Scenarios > step 2 (with explicit Save Project dialog mechanism + project name + creation-script-binding semantic note) | preserved-with-axis-clarity-expansion. The semantic note "Data Sync ON binds the script as the project's creation script" makes the test intent explicit (the test verifies creation-script rebuild behavior on reopen, which requires the script to BE the creation script — Data Sync ON is the binding mechanism). |
| 3. "Close all" | Scenarios > step 3 (with workspace-clear verification) | preserved. |
| 4. "Update any CSV file in the 'System:DemoFiles/chem' folder (add file with numeric suffix or rename any file (add numeric suffix) since the script sorts by numeric suffix)" | Scenarios > step 4 (with two explicit options: add OR rename + warning about destructive mutation) | preserved-with-axis-clarity-expansion. Original's parenthetical "(add file with numeric suffix or rename any file)" expanded to two explicit option bullets. Added explicit WARNING about destructive shared-state mutation. |
| 5. "Open the saved project - verify, that:" + 2 verification bullets | Scenarios > step 5 (with verification bullets preserved) | preserved. The implicit "Open the saved project" action + the 2 explicit verification bullets ("No errors occur", "The most recently created or modified file in the `chem` folder is loaded") all preserved. |
| (no trailing JSON) | (n/a) | metadata-not-step note: source has NO `{ "order": N }` JSON. Chain analysis rev 2 notes this and assigns alphabetical ordering. Migrated body's `## Notes` section documents this absence. |

No original numbered step or verification bullet is silently
dropped. The inline JS script in step 1 is preserved verbatim.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.custom-creation-scripts.md.target_layer = playwright`.
  The scenario exercises Save Project dialog (UI), close-and-reopen
  lifecycle (UI session), AND filesystem mutation (server-side via
  JS API). Mixed concern — `api-contract` cannot exercise the
  reopen UI flow; the Playwright layer handles both UI navigation
  and JS-API-driven filesystem mutation cleanly.
- **Why this `priority`:** chose `regression` per A-STRUCT-MECH-06
  enum (`smoke | regression | edge | perf`). The creation-script-
  rebuild-on-reopen behavior is regression-prone (the rebuild
  mechanism is server-side and depends on script-execution
  semantics that have historically been fragile). Per-cycle
  Invariant 1 honored.
- **Why this `strategy`:** `simple` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.custom-creation-scripts.md.strategy = simple`. The
  scenario is a single linear narrative (no matrix axes, no
  chained tests, no cross-fixture composite). Migrator schema's
  `simple` enum matches.
- **D-MERIT-01 compliance: NO Migrator-stage SCOPE_REDUCTION on
  destructive-mutation grounds.** The pre-known context flags
  this as "likely SCOPE_REDUCTION ... destructive shared-state
  mutation may not pass Critic Gate D without per-test working-
  copy of `System:DemoFiles/chem`". However: dropping step 4 from
  the migrated .md would BREAK the test's intent. The dynamic-
  data-resolution verification IS the test — without step 4's
  mutation between save and reopen, step 5's verification ("the
  most recently created or modified file is loaded") has nothing
  to verify against. Therefore, Migrator preserves all 5 steps
  faithfully + flags spec-time mitigations for the destructive
  mutation in `## Deferred items`. The "shared-state mutation"
  concern is a real technical dependency for SPEC-TIME mitigation
  (mitigated at Automator stage), NOT a Migrator-stage SCOPE_
  REDUCTION justification. Critic Gate D may pass this scenario
  even with the destructive operation present, because the
  mutation is intrinsic to test intent.
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - No `custom-creation-scripts-spec.ts` exists. Spec generation
    is Step 7.
  - Adjacent specs (`uploading-spec.ts`, `opening-spec.ts`,
    `browser-spec.ts`, `deleting-spec.ts`) confirm the section's
    `loginToDatagrok` / `softStep` / `evalJs` / `closeAll`
    convention; read-only inspection.
  - No sibling spec exercises filesystem mutation in the
    `System:DemoFiles/chem` path; this scenario is unique in that
    regard.
- **Helpers consulted / candidates:**
  - Reuse from prior scenarios: `helpers.playwright.projects.saveAndReopen`
    (from scenario 2) for steps 2-3 + step 5's reopen.
  - **NEW candidates** specific to this scenario:
    - `helpers.playwright.projects.runInlineScript(page, scriptText)` —
      execute a JS creation script via the platform's scripting UI,
      return the resulting dataframe (or its identifier) for
      downstream addViewers + Save Project chain.
    - `helpers.playwright.files.mutateChemFolder(page,
      action: 'add'|'rename', suffix: number)` — perform the
      destructive mutation in step 4 with explicit suffix control;
      complementary helper `restoreChemFolder(page, baseline)`
      MUST exist for `afterAll` cleanup.
- **Bug library consulted:** yes — `bug-library/projects.yaml`
  rev 2. No bug directly intersects this scenario's flow:
  - GROK-19403 (share with un-shared deps) — not exercised; no
    sharing in this scenario.
  - GROK-19728 (view-and-use failure state) — not exercised; no
    induced failure state, no view-and-use sharing.
  - GROK-19212 (rename + datasync, fails on reopen) — tangentially
    related (the script does NOT rename a referenced table; it
    selects a different CSV file based on filesystem state), but
    NOT the bug's reproduction path. The scenario's flow is
    creation-script-rebuild, not table-rename.
  - GROK-18345 (Spaces datasync share) — not exercised; no
    Spaces, no share.
  - GROK-19103 (join + save) — not exercised; no join.
  - GROK-19750 (save-copy with link) — not exercised; no save-copy.
  - github-3550 (external rename invalidation, queries axis) —
    not exercised; the script is inline (not a saved Query
    entity), and no entity rename happens.
  `related_bugs: []` is the correct call.
- **Decision log queried:** yes — `decision-log.yaml` rev 12 read.
  Cross-references:
  - `mig-2026-04-29-source-text-correction` does NOT apply
    (this scenario does not have "with layout" or "personal view
    customizations" content).
  - `mig-2026-04-29-fixture-placeholder` does NOT apply (no
    "second user" placeholder).
  - `atlas-2026-04-30-add-projects-shell-share-via-context-menu` /
    `atlas-2026-04-30-cascade-share-via-context-menu` do NOT
    apply (no Share dialog).
  - `b14-2026-04-30-re-auth-pattern-applied` does NOT apply (no
    re-auth).
  No prior decision contradicts this scenario's content.
- **Per-cycle override invariants (all 3) status:**
  - Invariant 1 (priority enum source-of-truth): honored —
    `priority: regression` is canonical per A-STRUCT-MECH-06.
  - Invariant 2 (existing -spec.ts / -api.ts READ-ONLY):
    SATISFIED — no `custom-creation-scripts-spec.ts` exists;
    trivially satisfied.
  - Invariant 3 (atlas-aware sub_features_covered for share):
    correctly does NOT apply — no Share dialog operation in this
    scenario. `projects.shell.share-via-context-menu` correctly
    OMITTED.

## Opt-outs (SCOPE_REDUCTION proposals)

(none) at Migrator stage. Per D-MERIT-01, opt-outs require a real
technical dependency, not effort. The destructive shared-state
mutation IS a real technical dependency, but the right level of
mitigation is SPEC-TIME (Automator owns per-test working-copy or
isolation), NOT a Migrator-stage drop. Dropping step 4 would
break test intent.

A spec-time mitigation IS proposed for Automator (see Deferred
items below); Critic Gate A may apply per-test isolation as a
spec-design decision (which IS Gate A's scope), without a
Migrator-stage opt-out being needed.

## Deferred items (NOT opt-outs)

- **Destructive shared-state mutation mitigation.** Step 4
  modifies `System:DemoFiles/chem` which is shared filesystem
  state. Real technical dependency; spec-time mitigation options:
  - **Option A: per-test working-copy folder.** Before step 1,
    copy `System:DemoFiles/chem` contents to a per-test temp
    folder (e.g. `System:Tmp/chem-test-<timestamp>`); rewrite
    the inline script's `chem` path to point at the temp folder;
    perform mutation in step 4 against the temp folder; clean
    up temp folder in `afterAll`. Test intent fully preserved;
    no shared-state collision.
  - **Option B: isolated test-environment teardown.** Run this
    scenario in an isolated test environment where
    `System:DemoFiles/chem` can be freely mutated; teardown the
    environment between runs. More heavyweight; appropriate if
    multiple destructive scenarios share the environment.
  - **Option C: test serialization + restore.** Run this
    scenario serially (no concurrent runs touching `chem`);
    capture `System:DemoFiles/chem` baseline before step 4;
    restore baseline in `afterAll`. Lightest weight but doesn't
    handle interleaved test failures cleanly.
  Automator picks one at spec time. Real prerequisite, not
  effort.
- **`System:DemoFiles/chem` folder pre-existing content.** The
  test assumes the folder contains at least one CSV with a
  numeric suffix. Real environmental prerequisite; Automator
  validates in `beforeAll` OR fixture-builds the folder content.
- **Inline-script execution UI mechanism.** Step 1's "Run the
  script" can be triggered via the platform's scripting UI OR
  programmatically via `grok.functions.eval(...)` /
  `grok.dapi.scripts.run(...)`. Migrated body documents both
  options; Automator picks per spec convenience. The Critic Gate
  D pass through doesn't depend on this choice.
- **Cleanup of saved project.** The produced project
  `Custom_Creation_Script_Test` is cleaned up by `deleting.md`
  (scenario 9, must_run_last). Real chain dependency.
- **Cleanup of filesystem mutation.** The CSV added or renamed
  in step 4 is NOT cleaned up by `deleting.md`. Automator's
  `afterAll` (or per-test environment teardown) must restore
  `System:DemoFiles/chem` to a known state. Tied to the
  mitigation chosen (Option A: drop temp folder; Option B:
  teardown env; Option C: explicit restore).

## Edge cases

The original lists no explicit edge cases. Implicit edge cases
derivable:

- **Empty `chem` folder edge.** The script throws "No CSV files
  found" if the folder is empty. Test must ensure non-empty
  folder before step 1.
- **No-numeric-suffix files in folder.** The `getNumberSuffix`
  function returns -1 for files without a numeric suffix. If ALL
  files in the folder lack numeric suffixes, the script picks
  the file with `-1` (alphabetically first per the sort
  comparison), which may be ambiguous. Implicit edge: the test
  assumes at least one file with a numeric suffix exists.
- **Suffix collision.** Two files with the same numeric suffix —
  the script picks one based on sort stability (which may not be
  deterministic across JS engines). Implicit edge: the test
  assumes unique suffixes.
- **Mutation timing — between save and reopen vs. before save.**
  Step 4 happens AFTER step 3 (Close All) and BEFORE step 5
  (reopen). The script SHOULD re-run on reopen and pick up the
  new file. If the project's creation-script binding caches the
  step-1 result instead of re-running on reopen, the test will
  fail (expected: new file loaded; actual: original step-1 file
  loaded). This IS the bug-detection axis the test exists to
  cover.
- **Filesystem-update propagation latency.** If the file system
  has caching/CDN latency, step 5's reopen may not yet see the
  step 4 mutation. Implicit edge: depends on environment.

(none additional)

## Unresolved ambiguities

- **"Update any CSV file" semantics in step 4.** The original
  offers two options ("add file with numeric suffix or rename
  any file (add numeric suffix)") — both achieve the same test
  goal (make a different file the highest-suffixed). Migrated
  body preserves both options as explicit bullets; Automator
  picks at spec time.
- **Specific file naming for the mutation.** The original is
  silent on what suffix to use. A safe choice is a number HIGHER
  than any existing suffix (e.g. 999 if existing range is 1-99).
  Migrated body uses 999 as the example. Automator may use a
  per-run timestamp if collision-resistance is needed.
- **Inline script storage location.** The script is inline in
  the `.md` step 1 — does it become a saved script entity (e.g.
  in the user's namespace), OR is it a one-shot execution that
  produces a dataframe but is not persistable as a Script
  entity? The "Save Project with Data Sync" in step 2 implies
  the script becomes the project's creation script (so it's
  bound to the project, not a standalone Script entity).
  Migrated body assumes the binding-as-creation-script
  interpretation; Automator verifies at spec time.
- **Step 1 "Run the script" in workspace context.** Does the
  script need to be run while a particular table is active, OR
  is it run in a fresh session? Original is silent. Migrated
  body assumes fresh session (no table prerequisite); Automator
  verifies.
- **Step 5 verification specificity for "no errors occur".** The
  original says "No errors occur" — what specifically is
  checked? Console errors? Error balloons? Failure dialogs?
  Migrated body lists "console errors (F12)" + "no error balloon"
  as explicit verifications; Automator may add deeper checks.
- **Atlas sub_feature naming.** The closest atlas match for
  creation-script-rebuild-on-reopen is
  `projects.url-params.rebuild-creation-script` — but this
  scenario does NOT use URL params. A non-URL-param sub_feature
  for creation-script flows (e.g.
  `projects.creation-script.rebuild-on-reopen`) would be a more
  precise match. Flagged as Section 1 candidate atlas addition;
  not applied autonomously.
