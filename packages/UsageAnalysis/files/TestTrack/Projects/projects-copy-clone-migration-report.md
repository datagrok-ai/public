# Migration Report — projects-copy-clone.md

## Step mapping

The original is 5 numbered steps + nested 4 sub-bullets in step 4 +
trailing JSON metadata `{ "order": 5 }`. The migrated body preserves
all original steps and sub-bullets, applies one MANDATORY addition
(GROK-19750 invariant assertion), and one OPTIONAL harmonization
(close-and-reopen verification on the personal-view-customizations
sub-bullet, for symmetry with sibling sub-bullets).

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Go to Browse and check the preview of all created projects" | Scenarios > step 1 (Browse + preview verification) | preserved |
| 2. "Share created projects" | Scenarios > step 2 (right-click Share, recipients, verify Sharing tab) | preserved-with-axis-clarity-expansion. Original's bare "Share created projects" expanded to explicit right-click Share + recipient fill + Sharing tab verification per D-STEP-02. Same Share-dialog mechanism as `share-project.md`. |
| 3. "Open created projects" | Scenarios > step 3 (open + verify on open) | preserved-with-verification-split. Original's bare "Open created projects" expanded to explicit open + render verification + no-console-errors. |
| 4. "Edit projects, save, and reopen:" + 4 sub-bullets | Scenarios > step 4 (with sub-flows 4a, 4b, 4c, 4d) | preserved-with-mandatory-addition (GROK-19750 invariant in 4b) and harmonization (close-and-reopen in 4d). |
| 4 sub-bullet 1: "add any viewer and save the original project, close all" | Scenarios > step 4 sub-flow **4a** (steps 1-4) | preserved-with-axis-clarity-expansion. Sub-bullet expanded into 4 atomic numbered steps for D-STEP-02. |
| 4 sub-bullet 2: "open the original project, add any viewer and save a copy with the **link**, close all, reopen it. Close all" | Scenarios > step 4 sub-flow **4b** (steps 1-8) | **preserved-with-mandatory-addition**. Sub-bullet preserved as steps 1-6 (open + add viewer + Save Copy Link + Close All + reopen-copy + verify-copy + Close All). **Step 7 added: GROK-19750 invariant assertion** — reopen the ORIGINAL project after the copy roundtrip and verify ORIGINAL's viewers remain intact. **Step 8: Close All.** This addition is the regression-coverage assertion for GROK-19750 (per bug-library rev 2: "Original project opens without viewers after Save Copy with table link mode"). The bug's reproduction is exactly this sub-flow's mechanics; verifying the original's intactness IS the regression test. Per the per-cycle prompt: "MUST add an explicit assertion ... mandatory migration-time addition, not a candidate". |
| 4 sub-bullet 3: "open the original project, add any viewer and save a copy with the **clone**, close all, reopen it. Close all" | Scenarios > step 4 sub-flow **4c** (steps 1-6) | preserved-with-axis-clarity-expansion. Sub-bullet expanded into 6 atomic numbered steps for D-STEP-02. NO mandatory addition for clone mode (GROK-19750 is specific to LINK mode; clone mode does not have the leak-back semantics). |
| 4 sub-bullet 4: "open the original project, add any viewer and save personal view customizations" | Scenarios > step 4 sub-flow **4d** (steps 1-6) | **preserved-with-harmonization-addition**. Sub-bullet preserved as steps 1-4. **Steps 5-6 added: close-and-reopen verification of the variant**, mirroring the sibling sub-flows 4b (steps 5-6) and 4c (steps 5-6). Sub-bullet 4 was the only one missing this verification; harmonization not a behavior change — just symmetry with sibling sub-flows. Per the per-cycle prompt's pre-known context: "Migrator may auto-add the missing close-and-reopen verification ... acceptable to apply autonomously". Auto-applied. |
| 5. "Share the new created projects. Check ability to open the shared projects" | Scenarios > step 5 (re-share + recipient open verification) | preserved-with-axis-clarity-expansion. Original's "Check ability to open the shared projects" expanded to explicit recipient session open + render verification per D-STEP-02. |
| Trailing JSON `{ "order": 5 }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/projects.yaml` rev 2 `order_from_files`) |

No original step or sub-bullet is silently dropped. The two additions
(GROK-19750 mandatory + harmonization) are explicitly documented in
both this report and the migrated body.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.projects-copy-clone.md.target_layer = playwright`. The
  scenario exercises right-click context menus (Share, Save Copy
  modes), Save Project dialog interactions, Browse > Dashboards
  navigation, project reopen lifecycle, and recipient-session share
  verification — all UI-driven. `api-contract` cannot exercise the
  multi-mode Save Copy UI flow.
- **Why this `priority`:** chose `regression` per A-STRUCT-MECH-06
  enum (`smoke | regression | edge | perf`). The scenario is multi-
  mode chained (3 save modes × open/save/reopen lifecycle) plus
  carries the GROK-19750 regression-coverage invariant — clearly
  regression-class, not a single golden-path smoke. Per-cycle
  Invariant 1 honored.
- **Why this `strategy`:** `chained_tests` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.projects-copy-clone.md.strategy = chained_tests`. The
  steps share state across the chain (open → preview → share → 4a →
  4b → 4c → 4d → re-share); no tear-down between steps. The
  `multi-source-saved-projects` fixture is reused across all chained
  steps. The 3 save modes are NOT independent matrix axes (they
  share state in chronological order — sub-flow 4b's GROK-19750
  invariant depends on the original's state being preserved across
  4a's save), so `data_driven` is incorrect. `end_to_end_fixtures`
  would imply per-step fixture rebuild, but here the same fixture is
  carried forward.
- **Sibling tests consulted (READ-ONLY per Invariant 2 of the
  per-cycle override):**
  - No `projects-copy-clone-spec.ts` exists. Spec generation is Step
    7 (Phase 1 plan separation).
  - Adjacent specs (`browser-spec.ts`, `deleting-spec.ts`,
    `opening-spec.ts`, `uploading-spec.ts`) all use the
    `loginToDatagrok`, `softStep`, `evalJs`, `closeAll` convention
    and `Date.now()` suffix naming pattern. Read-only inspection
    confirms; not modified.
  - `projects-copy-clone-run.md` exists from a prior MCP reproduction
    run; its run-log structure may inform the future Step 7 spec
    codegen but is not consumed by this Migrator stage.
- **Helpers consulted / candidate:** the migrated body identifies
  several candidate Playwright-layer helpers, surfaced in Section 3
  of the per-scenario REPORT for B14 propose-only flow:
  - `helpers.playwright.projects.shareProjectViaContextMenu(...)`
    (already candidate from scenario 3, reuse here for steps 2 + 5)
  - `helpers.playwright.projects.saveAndReopen(...)` (already
    candidate from scenario 2, reuse here for sub-flows 4a/4b/4c/4d)
  - `helpers.playwright.projects.saveCopyWithMode(page, name, mode)`
    where mode ∈ `{link, clone, personal-view-customizations}` —
    NEW candidate specific to this scenario; would parameterize
    sub-flows 4b/4c/4d cleanly.
- **Bug library consulted:** yes — `bug-library/projects.yaml`
  rev 2. **GROK-19750** intersects this scenario directly:
  - `affects: [projects.api.save, projects.add-relation]`
  - `reproduction:` (lines 9-17 of bug-library) is exactly the
    sequence of sub-flow 4b: open original → add viewer → save
    project → close+reopen (verifies the viewer persists) → add
    another viewer → Save Copy with **Link** mode → close all →
    open original → original opens WITHOUT its viewer (the bug).
  - Listed in `related_bugs: [GROK-19750]` in frontmatter.
  - The migrated body's sub-flow 4b step 7 ("GROK-19750 INVARIANT")
    is the regression test for this bug.
  Other bugs (GROK-19212, GROK-19103, GROK-19403, GROK-18345,
  GROK-19728) do not intersect this scenario's flows directly:
  GROK-19212 is rename + sync (not Save Copy modes); GROK-19103 is
  join + save (not Save Copy); GROK-19403, GROK-18345, GROK-19728
  are share + recipient-side load failures (this scenario's share
  steps are simple share without unshared deps / Spaces / failure
  states). Steps 2 and 5 share-flow could potentially trigger
  GROK-19403 IF the variants depended on un-shared scripts/queries
  — but the variants here are pure dataset projects without script
  deps, so the bug's reproduction path is not triggered.
- **Decision log queried:** yes — `decision-log.yaml` rev 9 read.
  Three prior entries directly apply:
  - `mig-2026-04-29-source-text-correction`: "with layout" →
    "personal view customizations" terminology was applied to this
    scenario's source step 4 sub-bullet 4. Migrated body preserves
    the corrected wording verbatim.
  - `mig-2026-04-29-fixture-synthesized-inline`: confirms "save
    personal view customizations" is a save-flow flag produced by
    THIS scenario (no external fixture needed). Migrated body
    treats it as a save-mode within sub-flow 4d.
  - `atlas-2026-04-30-add-projects-shell-share-via-context-menu`
    + `atlas-2026-04-30-cascade-share-via-context-menu`: the new
    sub_feature was added to atlas; this scenario's `sub_features_
    covered` correctly INCLUDES `projects.shell.share-via-context-
    menu` per Invariant 3.
- **Per-cycle override invariants (all 3) status:**
  - Invariant 1 (priority enum source-of-truth): honored —
    `priority: regression` is canonical per A-STRUCT-MECH-06.
  - Invariant 2 (existing -spec.ts / -api.ts READ-ONLY):
    SATISFIED — no `projects-copy-clone-spec.ts` exists. Sibling
    specs read for convention only and not modified.
  - Invariant 3 (atlas-aware sub_features_covered for share):
    APPLIES — steps 2 and 5 both exercise the right-click Share
    dialog. `projects.shell.share-via-context-menu` is INCLUDED in
    `sub_features_covered`.

## Opt-outs (SCOPE_REDUCTION proposals)

(none) — every original numbered step + sub-bullet is preserved at
the same target layer. The two additions (GROK-19750 invariant in
4b; harmonization close-and-reopen in 4d) are EXPANSIONS of coverage,
not reductions.

## Deferred items (NOT opt-outs)

- **`multi-source-saved-projects` fixture availability.** This
  scenario consumes saved projects from upstream (`upload-project.md`
  for `demog`, optionally `uploading.md` for `Test_Case<N>_Sync`/
  `_NoSync`). Migration-stage cannot guarantee fixture presence at
  spec-run time; Automator stage owns `beforeAll` fixture build.
- **Recipient identity for steps 2 and 5.** Same axis as
  `share-project.md`: registered user (e.g. `Olena Ahadzhanian`)
  AND/OR email recipient (auto-creates user account; cleanup in
  `afterAll`). Automator decides between hardcoded recipients,
  pre-created user-pool, or per-run timestamp emails.
- **"Add any viewer" semantics across sub-flows.** Each of 4a/4b/4c/
  4d says "add any viewer". Automator decides at spec time which
  specific viewer (scatter, bar, line, histogram) to add per sub-
  flow. The choice does NOT affect the regression test for GROK-19750
  (the invariant is about ANY viewer's preservation, not a specific
  type), but consistency across sub-flows aids assertion clarity.
- **"Personal view customizations" semantics in 4d.** The original
  scenario does not define WHICH customizations to apply (filter?
  sort? column visibility? layout reposition?). The migrated body
  uses "viewer addition OR a personal view customization" as a
  flexible interpretation. Automator picks specific customizations
  at spec time per the `mig-2026-04-29-fixture-synthesized-inline`
  guidance (the customizations are a save-mode flag, not a specific
  fixture state).
- **Side-effect cleanup of variant projects.** Steps 4b/4c/4d each
  produce a NEW variant project. These are CONSUMED by
  `project-url.md` (chain rev 2 dependency); they should NOT be
  cleaned up at `afterAll` of THIS scenario — `deleting.md`
  (scenario 9, must_run_last) owns terminal cleanup. Real chain
  dependency, not effort.

## Edge cases

The original lists no explicit edge cases. Implicit edge cases
derivable from the steps:

- **GROK-19750 regression itself.** Sub-flow 4b step 7 explicitly
  tests this. If the original loses viewers after Save-Copy-with-
  Link, the test FAILS — regression coverage achieved.
- **Save Copy with Link vs. Clone semantic distinction.** Sub-flows
  4b and 4c produce semantically different copies: Link copies
  share data with the original; Clone copies are independent. On
  reopen, the link variant should reflect changes to the original's
  data; clone should not. The migrated body's verification ("tables
  render via the link to the original's data" for 4b vs. "tables
  are independent copies" for 4c) preserves this distinction.
- **Personal view customizations apply but data is shared.** Sub-
  flow 4d's variant has customized view-state but shares data with
  the original (per save-mode semantics). On reopen, customizations
  must apply but data tables must be the original's. Implicit; the
  verification step is preserved.
- **Re-share permission propagation (step 5).** Each variant is
  shared with recipients in step 5. The recipients must have access
  to the variant's underlying data (which, for Link mode, is the
  original's data). If the original is NOT shared with the same
  recipient, the recipient may hit a permission failure when trying
  to access linked data. This is the GROK-19403 reproduction path
  — but only if the variant has un-shared dependencies. Pure
  dataset projects without script/query deps avoid this; the
  scenario's verification step 5 is preserved as-is.
- **Order of variant creation matters for GROK-19750.** Sub-flow 4a
  must complete BEFORE 4b for the GROK-19750 invariant to be
  testable (the original needs at least one viewer before sub-flow
  4b adds another and saves a Link copy). The chained_tests
  ordering enforces this.

(none additional)

## Unresolved ambiguities

- **"All created projects" in step 1.** The original's step 1 says
  "all created projects" — referring to upstream-produced projects.
  The migrated body interprets this as: at minimum the `demog`
  project from `upload-project.md`, optionally the 18 from
  `uploading.md`. Automator decides which subset to iterate over;
  the GROK-19750 invariant is testable against any single one (the
  `demog` project is the minimum sufficient case).
- **"Save Copy" UI control location.** The original assumes a
  "Save Copy" UI control with mode selection (link / clone /
  personal-view-customizations). The exact UI placement (File menu?
  Project context menu? Save dialog dropdown?) is not specified.
  Automator must cross-reference `grok-browser/references/` for
  Save Copy dialog and mode-selection control selectors at spec
  time.
- **"Save with personal view customizations" mode trigger.** The
  source-text correction (`mig-2026-04-29-source-text-correction`)
  established that "save personal view customizations" is a
  save-mode flag, but the exact UI trigger (separate menu item?
  toggle in Save dialog? specific button labeled "Save with
  customizations"?) is not specified. Automator must verify the UI
  control name at spec time.
- **Recipient open verification depth (step 5).** The original
  says "Check ability to open the shared projects" without
  specifying what depth of verification. The migrated body
  expanded this to "recipient session opens + render + no console
  errors" but does NOT verify deep behavior (linked data updates
  propagate? customizations preserved cross-user?). Automator may
  add deeper assertions at spec time.
- **GROK-19750 invariant: which viewers count as "intact"?** The
  added sub-flow 4b step 7 verifies "the original project's viewers
  ... are still present and intact". Strict interpretation: count
  = original-pre-sub-flow-4b + (sub-flow 4b's additional viewer if
  it was added BEFORE the Save Copy was triggered, since save-copy
  preserves the current state). Looser interpretation: at minimum
  the original's pre-existing viewers must persist. Per bug-library
  rev 2 reproduction, the bug's failure mode is the original
  losing the viewer added in sub-flow 4b's step 2 — so the strict
  interpretation is the regression test. Migrated body's wording
  ("viewers ... existing before sub-flow 4b plus the sub-flow 4b
  addition") matches the strict interpretation. Flag for retro:
  is the strict interpretation universally agreed?
