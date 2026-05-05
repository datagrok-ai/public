# Migration Report — projects-copy-clone.md

Re-migration under chain YAML revision 3 (2026-05-04). Supersedes the
prior 2026-04-30 migration. Schema and content changes from rev 2 → 3:

- Frontmatter `pyramid_layer: bug-focused` added (chain rev 3 Edit 2,
  Rule 3 — GROK-19750 regression-coverage discriminator).
- Frontmatter `ui_coverage_responsibility:` populated with 6 flows
  owned by this scenario (chain rev 3 Edit 1; includes the right-
  click Share dialog flows delegated from `share-project.md` per
  `ui_coverage_plan.delegated_scenarios`).
- Frontmatter `ui_coverage_delegated_to: null` (this scenario owns
  its UI surface; delegates to no one).
- Frontmatter `related_bugs: [GROK-19750, GROK-19103, GROK-19403]`
  expanded from rev 2 `[GROK-19750]` to reflect chain rev 3
  `bug_focused_candidates` cross-cutting cover where this scenario
  is a span participant.
- Body Notes updated to surface (a) `pyramid_layer: bug-focused`
  discriminator citation, (b) the 6 owned UI flows + delegation
  rationale, (c) cross-cutting bug citations for GROK-19750 /
  GROK-19103 / GROK-19403 per chain rev 3
  `bug_focused_candidates`.
- All five original numbered steps (and the 4 step-4 sub-bullets)
  preserved verbatim from the rev 2 migration; the 2026-04-30
  GROK-19750 mandatory addition (sub-flow 4b step 7) and the
  harmonization addition (sub-flow 4d steps 5-6) are retained
  unchanged. No silent drops.

## Step mapping

The original is 5 numbered steps + nested 4 sub-bullets in step 4 +
trailing JSON metadata `{ "order": 5 }`. The migrated body preserves
all original steps and sub-bullets, applies one MANDATORY addition
(GROK-19750 invariant assertion, retained from rev 2), and one
OPTIONAL harmonization (close-and-reopen verification on the
personal-view-customizations sub-bullet, for symmetry with sibling
sub-bullets, retained from rev 2).

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Go to Browse and check the preview of all created projects" | Scenarios > step 1 (Browse + preview verification) | preserved |
| 2. "Share created projects" | Scenarios > step 2 (right-click Share, recipients, verify Sharing tab) | preserved (split for clarity). Original's bare "Share created projects" expanded to explicit right-click Share + recipient fill + Sharing tab verification per D-STEP-02. Same Share-dialog mechanism as `share-project.md`; this scenario OWNS the right-click Share dialog UI per chain rev 3 `ui_coverage_plan.delegated_scenarios`. |
| 3. "Open created projects" | Scenarios > step 3 (open + verify on open) | preserved as verification. Original's bare "Open created projects" expanded to explicit open + render verification + no-console-errors. |
| 4. "Edit projects, save, and reopen:" + 4 sub-bullets | Scenarios > step 4 (with sub-flows 4a, 4b, 4c, 4d) | preserved (split for clarity) with mandatory addition (GROK-19750 invariant in 4b) and harmonization (close-and-reopen in 4d). |
| 4 sub-bullet 1: "add any viewer and save the original project, close all" | Scenarios > step 4 sub-flow **4a** (steps 1-4) | preserved (split for clarity). Sub-bullet expanded into 4 atomic numbered steps for D-STEP-02. |
| 4 sub-bullet 2: "open the original project, add any viewer and save a copy with the **link**, close all, reopen it. Close all" | Scenarios > step 4 sub-flow **4b** (steps 1-8) | preserved (split for clarity, mandatory addition). Sub-bullet preserved as steps 1-6 (open + add viewer + Save Copy Link + Close All + reopen-copy + verify-copy + Close All). **Step 7 added: GROK-19750 invariant assertion** — reopen the ORIGINAL project after the copy roundtrip and verify ORIGINAL's viewers remain intact. **Step 8: Close All.** This addition is the regression-coverage assertion for GROK-19750 (per bug-library rev 2 reproduction). Per per-cycle prompt: "MUST add an explicit assertion ... mandatory migration-time addition, not a candidate". |
| 4 sub-bullet 3: "open the original project, add any viewer and save a copy with the **clone**, close all, reopen it. Close all" | Scenarios > step 4 sub-flow **4c** (steps 1-6) | preserved (split for clarity). Sub-bullet expanded into 6 atomic numbered steps for D-STEP-02. NO mandatory addition for clone mode (GROK-19750 is specific to LINK mode). |
| 4 sub-bullet 4: "open the original project, add any viewer and save personal view customizations" | Scenarios > step 4 sub-flow **4d** (steps 1-6) | preserved (split for clarity, harmonization addition). Sub-bullet preserved as steps 1-4. **Steps 5-6 added: close-and-reopen verification of the variant**, mirroring the sibling sub-flows 4b (steps 5-6) and 4c (steps 5-6). Sub-bullet 4 was the only one missing this verification; harmonization is not a behavior change — just symmetry with sibling sub-flows. Per per-cycle prompt's pre-known context: "Migrator may auto-add the missing close-and-reopen verification ... acceptable to apply autonomously". Auto-applied. |
| 5. "Share the new created projects. Check ability to open the shared projects" | Scenarios > step 5 (re-share + recipient open verification) | preserved (split for clarity). Original's "Check ability to open the shared projects" expanded to explicit recipient session open + render verification per D-STEP-02. |
| Trailing JSON `{ "order": 5 }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/projects.yaml` rev 3 `order_from_files`) |

No original step or sub-bullet is silently dropped. The two additions
(GROK-19750 mandatory + harmonization) are explicitly documented in
both this report and the migrated body and are retained verbatim
from the rev 2 migration.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 3
  `output_plan.projects-copy-clone.md.target_layer = playwright`. The
  scenario exercises right-click context menus (Share, Save Copy
  modes), Save Project dialog interactions, Browse > Dashboards
  navigation, project reopen lifecycle, and recipient-session share
  verification — all UI-driven. `api-contract` cannot exercise the
  multi-mode Save Copy UI flow.
- **Why this `coverage_type`:** chose `regression` per the
  `smoke | regression | edge | perf` enum. The scenario is multi-
  mode chained (3 save modes × open/save/reopen lifecycle) plus
  carries the GROK-19750 regression-coverage invariant — clearly
  regression-class, not a single golden-path smoke.
- **Why this `pyramid_layer: bug-focused`:** per chain rev 3
  Rule 3 discriminator. `related_bugs: [GROK-19750, GROK-19103,
  GROK-19403]`; sub-flow 4b step 7 explicitly walks GROK-19750's
  reproduction path (open original → add viewer → Save Copy with
  Link → close all → reopen original → assert viewers intact) and
  IS the regression test. Discriminator test passes: GROK-19750
  fails this scenario before fix. Bug-focused dominates over Rule 4
  multi-source per chain rev 3 heuristic priority 3 → 4.
- **Why this `strategy`:** `chained_tests` per
  `scenario-chains/projects.yaml` rev 3
  `output_plan.projects-copy-clone.md.strategy = chained_tests`. The
  steps share state across the chain (open → preview → share → 4a →
  4b → 4c → 4d → re-share); no tear-down between steps. The
  `multi-source-saved-projects` fixture is reused across all chained
  steps. The 3 save modes are NOT independent matrix axes (sub-flow
  4b's GROK-19750 invariant depends on the original's state being
  preserved across 4a's save), so `data_driven` is incorrect.
- **UI coverage ownership (chain rev 3):** this scenario owns 6 UI
  flows in `ui_coverage_responsibility`:
  1. `save-copy-with-link-dialog` (sub-flow 4b)
  2. `save-copy-with-clone-dialog` (sub-flow 4c)
  3. `save-personal-view-customizations-dialog` (sub-flow 4d)
  4. `pcmdShareProject` (steps 2 + 5; right-click Share dialog —
     **delegated from `share-project.md`** per chain rev 3
     `ui_coverage_plan.delegated_scenarios`)
  5. `share-dialog-recipients` (steps 2 + 5)
  6. `context-panel-sharing-tab` (step 2 verification)
  `ui_coverage_delegated_to: null` — owns its UI surface. The
  delegation citation `share-project.md → projects-copy-clone.md`
  appears in `share-project.md`'s `ui_coverage_delegated_to:` field
  (verified 2026-05-04). UI-only verification flows (preview
  thumbnail rendering quality at step 1; view-state customization
  preservation at sub-flow 4d step 5) are SCOPE_REDUCED to
  `projects-copy-clone-ui.md` — see Opt-outs section.
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `projects-copy-clone-spec.ts` EXISTS at the same path (created
    by a prior cycle). Per-cycle Invariant 2 enforces READ-ONLY on
    existing -spec.ts / -api.ts; Migrator does not modify it.
    `existing-test-index.yaml` lists it with
    `features_covered: [projects.add-link, projects.add-relation,
    projects.api.save, projects.shell.open]` and helpers
    `[spec-login]`.
  - Adjacent specs (`opening-spec.ts`, `uploading-spec.ts`,
    `share-project-spec.ts`, `deleting-spec.ts`) all use the
    `loginToDatagrok`, `softStep`, `evalJs`, `closeAll` convention
    and `Date.now()` suffix naming pattern; consulted READ-ONLY
    for house style.
  - `projects-copy-clone-run.md` exists from a prior MCP
    reproduction run; its run-log structure may inform future spec
    codegen but is not consumed by this Migrator stage.
- **Helpers consulted / candidate:** none of the candidates below
  exist in `helpers-registry.yaml` rev 1 (registered helpers are
  `grok_test_layer` `gui-utils.ts` style — `uploadProject`,
  `readDataframe`, `findViewer`, `waitForElement`, etc.).
  Surfaced as Migrator candidates for downstream addition:
  - `helpers.playwright.projects.shareProjectViaContextMenu(...)`
    (candidate from scenario 3 reuse — steps 2 + 5)
  - `helpers.playwright.projects.saveAndReopen(...)` (candidate from
    scenario 2 reuse — sub-flows 4a/4b/4c/4d)
  - `helpers.playwright.projects.saveCopyWithMode(page, name, mode)`
    where mode ∈ `{link, clone, personal-view-customizations}` —
    NEW candidate; would parameterize sub-flows 4b/4c/4d cleanly.
  Per Migrator helpers discipline (read-only on registry), surfaced
  as candidates only; Migrator does NOT modify the registry.
- **Bug library consulted:** yes — `bug-library/projects.yaml`
  rev 2. Per chain rev 3 `bug_focused_candidates`, this scenario is
  a span participant for:
  - **GROK-19750** (`projects-grok-19750-spec.ts` proposed; spans
    `upload-project.md:Step 1` + `projects-copy-clone.md:Step 4`).
    `affects: [projects.api.save, projects.add-relation]`
    intersects this scenario's `sub_features_covered`. Sub-flow 4b
    step 7 IS the regression assertion. Cross-cutting invariant:
    original-project-viewers-intact post-Save-Copy-with-Link.
    PRIMARY witness within the chain.
  - **GROK-19103** (`projects-grok-19103-spec.ts` proposed; spans
    `upload-project.md:Step 1` + `uploading.md:Step 6` +
    `projects-copy-clone.md:Step 4` + `complex.md:Step 1`).
    `affects: [projects.api.save, projects.add-relation]`
    intersects. Cross-cutting invariant: derivation lands in
    active project, not as stray new project. This scenario covers
    the save-after-edit-active-project surface (step 4); chain
    Trigger 1 cross-scenario emission.
  - **GROK-19403** (`projects-grok-19403-spec.ts` proposed; spans
    `share-project.md:Step 4` + `complex.md:Step 12` +
    `projects-copy-clone.md:Step 5`).
    `affects: [projects.api.get-by-id, projects.add-relation]`
    intersects via `projects.add-relation` and `projects.add-link`.
    Cross-cutting invariant: share project with un-shared script
    dependency → recipient sees explicit failure, not silent null.
    Step 5 share + recipient open is the candidate cover surface.
  - All three listed in `related_bugs:` frontmatter. Per Edit 5
    Design point A, the citation is RECOMMENDED (early visibility);
    F-BUG-COVERAGE-01 at section-complete is the authoritative gate
    for cross-cutting bug coverage at section-level.
  - Other curated bugs (GROK-19212, GROK-18345, GROK-19728,
    github-3550) do not have this scenario in their chain
    `bug_focused_candidates` spans → not in `related_bugs:`.
- **Decision log queried:** yes — `decision-log.yaml` read; entry
  `mig-2026-04-30-projects-copy-clone-migration` honored
  (rev 2 migration outcomes preserved verbatim: GROK-19750 invariant
  at sub-flow 4b step 7, harmonization at sub-flow 4d, all 5 step-
  mapping decisions, all 9 deferred-items decisions). Other relevant
  prior entries:
  - `mig-2026-04-29-source-text-correction`: "with layout" →
    "personal view customizations" terminology preserved verbatim.
  - `mig-2026-04-29-fixture-synthesized-inline`: "save personal
    view customizations" treated as a save-mode flag produced by
    sub-flow 4d (no external fixture).
  - `atlas-2026-04-30-add-projects-shell-share-via-context-menu`
    + `atlas-2026-04-30-cascade-share-via-context-menu`: atlas
    sub_feature `projects.shell.share-via-context-menu` correctly
    INCLUDED in `sub_features_covered`.
  - `mig-2026-04-30-priority-enum-drift-discovery`: scenario
    frontmatter uses `coverage_type: regression` per the
    plan-01-canonical enum (not the SKILL.md rev that briefly
    declared `priority` — drift reverted).
  No `failed_attempts` entries reference this scenario directly.

## Opt-outs (SCOPE_REDUCTION proposals)

- **Preview thumbnails render-quality check (step 1).** SCOPE_REDUCED
  to `projects-copy-clone-ui.md` — UI-only visual judgment. The
  Playwright spec cannot reliably verify thumbnail visual rendering
  without a screenshot-diff fixture (technical dependency on a
  screenshot-diff tooling layer not present in the test
  infrastructure). UI coverage of the preview thumbnails lives in
  `projects-copy-clone-ui.md` (declared in chain rev 3
  `extracted_ui_flows` and the section's `ui_coverage_plan` —
  no scenario currently owns thumbnail render-quality at the
  Playwright layer; flagged for Gate F section-level review with
  the existing `-ui.md` companion as the manual coverage anchor).
- **View-state customization visual preservation (sub-flow 4d step 5).**
  SCOPE_REDUCED to `projects-copy-clone-ui.md` — UI-only visual
  judgment. Confirming filter/sort/column-visibility/layout-position
  preservation requires DOM-level visual comparison that is
  selector-fragile in Playwright (technical dependency on the
  `view-state-snapshot` helper which does not exist in
  `helpers-registry.yaml`). UI coverage gap: the
  `save-personal-view-customizations-dialog` flow IS owned in
  Playwright (variant loads + tables present + no console errors);
  the visual customization-preservation aspect is UI-only and
  delegated to `projects-copy-clone-ui.md`. Flagged for Gate F
  section-level review.

Both opt-outs cite real technical dependencies (D-MERIT-01
satisfied), not effort.

## Deferred items (NOT opt-outs)

- **`multi-source-saved-projects` fixture availability.** This
  scenario consumes saved projects from upstream (`upload-project.md`
  for `demog`, optionally `uploading.md` for `Test_Case<N>_Sync`/
  `_NoSync`). Migration-stage cannot guarantee fixture presence at
  spec-run time; Automator stage owns `beforeAll` fixture build per
  chain rev 3 `output_plan.strategy = chained_tests`.
- **Recipient identity for steps 2 and 5.** Same axis as
  `share-project.md`: registered user (e.g. `Olena Ahadzhanian`)
  AND/OR email recipient (auto-creates user account; cleanup in
  `afterAll`). Automator decides between hardcoded recipients,
  pre-created user-pool, or per-run timestamp emails. Prerequisite:
  recipient-provisioning policy at `helpers-registry.yaml` (helper
  `loginAsUser` / equivalent not yet registered).
- **"Add any viewer" semantics across sub-flows.** Each of 4a/4b/4c/
  4d says "add any viewer". Automator decides at spec time which
  specific viewer (scatter, bar, line, histogram) to add per sub-
  flow. The choice does NOT affect the regression test for GROK-19750
  (the invariant is about ANY viewer's preservation, not a specific
  type), but consistency across sub-flows aids assertion clarity.
  Prerequisite: viewer-add helper standardization in
  `helpers-registry.yaml`.
- **"Personal view customizations" semantics in 4d.** The original
  scenario does not define WHICH customizations to apply (filter?
  sort? column visibility? layout reposition?). The migrated body
  uses "viewer addition OR a personal view customization" as a
  flexible interpretation. Automator picks specific customizations
  at spec time per the `mig-2026-04-29-fixture-synthesized-inline`
  guidance. Prerequisite: customization-axis decision at Automator
  Step 7 — cannot be resolved at Migrator Step 6.
- **Side-effect cleanup of variant projects.** Steps 4b/4c/4d each
  produce a NEW variant project. These are CONSUMED by
  `project-url.md` (chain rev 3 dependency); they should NOT be
  cleaned up at `afterAll` of THIS scenario — `deleting.md`
  (terminal, must_run_last) owns terminal cleanup per chain rev 3
  `dependency_graph.deleting.md.depends_on` which includes
  `projects-copy-clone.md`. Prerequisite: chain orchestrator
  enforces topological ordering at Step 7 per D-STRUCT-01.

All deferred items cite a real prerequisite (D-MERIT-02 satisfied).

## Edge cases

The original lists no explicit edge cases. Implicit edge cases
derivable from the steps:

- **GROK-19750 regression itself.** Sub-flow 4b step 7 explicitly
  tests this. Preserved as scenario step in the migrated body. If
  the original loses viewers after Save-Copy-with-Link, the test
  FAILS — regression coverage achieved.
- **Save Copy with Link vs. Clone semantic distinction.** Sub-flows
  4b and 4c produce semantically different copies: Link copies
  share data with the original; Clone copies are independent. On
  reopen, the link variant should reflect changes to the original's
  data; clone should not. Preserved as scenario step in the
  migrated body's verification ("tables render via the link to the
  original's data" for 4b vs. "tables are independent copies" for
  4c).
- **Personal view customizations apply but data is shared.** Sub-
  flow 4d's variant has customized view-state but shares data with
  the original (per save-mode semantics). On reopen, customizations
  must apply but data tables must be the original's. Preserved as
  scenario step in sub-flow 4d step 5 verification.
- **Re-share permission propagation (step 5).** Each variant is
  shared with recipients in step 5. The recipients must have access
  to the variant's underlying data (which, for Link mode, is the
  original's data). If the original is NOT shared with the same
  recipient, the recipient may hit a permission failure when trying
  to access linked data. This intersects GROK-19403's reproduction
  path. Preserved as scenario step (step 5 recipient open
  verification); cross-cutting full coverage flagged for the
  proposed `projects-grok-19403-spec.ts` cell (not this scenario).
- **Order of variant creation matters for GROK-19750.** Sub-flow 4a
  must complete BEFORE 4b for the GROK-19750 invariant to be
  testable (the original needs at least one viewer before sub-flow
  4b adds another and saves a Link copy). The chained_tests
  ordering enforces this. Preserved as scenario structure
  (sub-flows 4a → 4b → 4c → 4d sequential).

(none additional)

## Unresolved ambiguities

- **"All created projects" in step 1.** The original's step 1 says
  "all created projects" — referring to upstream-produced projects.
  The migrated body interprets this as: at minimum the `demog`
  project from `upload-project.md`, optionally the 18 from
  `uploading.md`. Automator decides which subset to iterate over;
  the GROK-19750 invariant is testable against any single one (the
  `demog` project is the minimum sufficient case). gap_type:
  fixture-scope.
- **"Save Copy" UI control location.** The original assumes a
  "Save Copy" UI control with mode selection (link / clone /
  personal-view-customizations). The exact UI placement (File menu?
  Project context menu? Save dialog dropdown?) is not specified.
  Automator must cross-reference `grok-browser/references/` for
  Save Copy dialog and mode-selection control selectors at spec
  time. gap_type: ui-selector.
- **"Save with personal view customizations" mode trigger.** The
  source-text correction (`mig-2026-04-29-source-text-correction`)
  established that "save personal view customizations" is a
  save-mode flag, but the exact UI trigger (separate menu item?
  toggle in Save dialog? specific button labeled "Save with
  customizations"?) is not specified. Automator must verify the UI
  control name at spec time. gap_type: ui-selector.
- **Recipient open verification depth (step 5).** The original
  says "Check ability to open the shared projects" without
  specifying what depth of verification. The migrated body
  expanded this to "recipient session opens + render + no console
  errors" but does NOT verify deep behavior (linked data updates
  propagate? customizations preserved cross-user?). Automator may
  add deeper assertions at spec time. gap_type: assertion-depth.
- **GROK-19750 invariant: which viewers count as "intact"?** The
  added sub-flow 4b step 7 verifies "the original project's viewers
  ... are still present and intact". Strict interpretation: count
  = original-pre-sub-flow-4b + (sub-flow 4b's additional viewer if
  added BEFORE the Save Copy was triggered). Looser: at minimum the
  original's pre-existing viewers must persist. Per bug-library
  rev 2 reproduction, the bug's failure mode is the original
  losing the viewer added in sub-flow 4b's step 2 — strict
  interpretation is the regression test. Migrated body's wording
  ("viewers ... existing before sub-flow 4b plus the sub-flow 4b
  addition") matches the strict interpretation. gap_type:
  assertion-strictness.
- **GROK-19403 cross-cutting cover via step 5.** The chain rev 3
  `bug_focused_candidates` GROK-19403 entry includes
  `projects-copy-clone.md:Step 5` as a span with rationale "share
  project with un-shared script dependency". This scenario's
  variants are pure dataset projects (no script/query deps), so
  the bug's reproduction path is NOT directly triggered here —
  the cross-cutting cover is by sub_feature intersection
  (`projects.add-relation`, `projects.add-link`), not by repro
  path. The authoritative GROK-19403 regression cover lives in the
  proposed `projects-grok-19403-spec.ts` cell where script-dep
  variants are exercised. gap_type: bug-coverage-routing.
