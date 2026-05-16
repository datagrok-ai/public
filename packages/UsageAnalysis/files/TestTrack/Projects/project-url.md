---
feature: projects
sub_features_covered:
  - projects.url-params.build-share-link
  - projects.url-params.apply
  - projects.shell.open
  - projects.view.browse
target_layer: playwright
coverage_type: regression
pyramid_layer: integration
ui_coverage_responsibility:
  - context-panel-links-url-copy
  - new-tab-open-url
ui_coverage_delegated_to: null
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/project-url.md
migration_date: 2026-05-04
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - click-the-project-single-click-vs-double-click
  - context-panel-links-location
  - url-format-query-parameters
  - new-tab-vs-incognito-tab
  - order-vs-dependency-contradiction
  - source-text-correction-recurrence
scope_reductions: []
related_bugs: []
---

# Project URL

For the `demog` representative project (file-share source class — Variant
C representative source per chain rev 3 `pyramid_layer: integration`,
source-agnostic Context-Panel-Links URL deep-link reopen flow), navigate
to Browse > Dashboards, locate each of the four available variants
(original; copied-with-link; copied-with-clone; saved-with-personal-view-
customizations), copy the deep-link URL from Context Panel > Links, then
open the URL in a new browser tab and verify the corresponding project
loads.

This scenario depends on the composite fixture
`copy-clone-customizations-variants` (the three copy-mode variants
produced by `projects-copy-clone.md`) plus the upstream "original"
project — `demog` from `upload-project.md` is the Variant C
representative source for this scenario; `Test_Case<N>_Sync` /
`_NoSync` from `uploading.md` may also satisfy the "original" role
when the chain runs the matrix path.

## Setup

1. **Fixture prerequisite:** the four project variants must exist on
   the server:
   - **original**: the saved `demog` project from `upload-project.md`
     (Variant C representative source — file-share). The
     `Test_Case<N>_Sync` projects from `uploading.md` are an
     acceptable substitute when `demog` is unavailable in the chain
     run.
   - **copied-with-link**: the `<original>-link` Save-Copy-with-Link
     variant produced by `projects-copy-clone.md`.
   - **copied-with-clone**: the `<original>-clone` Save-Copy-with-Clone
     variant produced by `projects-copy-clone.md`.
   - **saved-with-personal-view-customizations**: the
     `<original>-personal-view-customizations` Save-personal-view-
     customizations variant produced by `projects-copy-clone.md`.
   In `end_to_end_fixtures` strategy, the Automator stage builds the
   composite `copy-clone-customizations-variants` fixture in a
   `beforeAll` block by either reusing fixtures from prior chain
   runs or replaying the relevant save-flow cases via `js-api-replay`
   on top of the `demog-project-with-viewers` baseline fixture.
2. The browser session is authenticated as the project owner (the
   user who produced the four variants upstream). Cross-user share
   verification is OUT of scope here — `share-project.md` covers that.

## Scenarios

### URL deep-link reopen for each variant

The following workflow is parameterised over the 4 project variants.
Implemented per `end_to_end_fixtures` strategy at the Playwright layer
(per `output_plan.project-url.md.strategy = end_to_end_fixtures` in
`scenario-chains/projects.yaml` rev 3). Variant C source-agnostic
discipline: the URL-build / URL-apply / shell-open path is identical
across all four variants; the assertion verifies uniformity, not
source-class-specific behavior.

For each `<variant>` in the 4-variant list:

1. Go to **Browse** > **Dashboards**.
2. Click the project corresponding to `<variant>`:
   - original project
   - copied with the link
   - copied with clone
   - saved with personal view customizations
3. Go to **Context Panel** > **Links** and copy the URL shown there.
4. Open a new tab in the browser and paste/navigate to the copied URL.
5. **Verify on URL load (new tab):** the project corresponding to
   `<variant>` opens — its tables, viewers, and view layout render in
   the new tab; no errors in the browser console (F12).

## Notes

- **Variant C representative source — `demog` (file-share).** Per
  chain rev 3 `pyramid_layer: integration` rationale (Rule 4 Variant
  C), this scenario is source-agnostic — the Context-Panel-Links URL
  copy + new-tab URL-apply path works identically across all source
  classes (files, query, script, spaces, db_table, derived). The
  representative source picked for the test is **file-share
  (`demog`)** because `upload-project.md` produces it as the smallest,
  cheapest baseline fixture (`demog-project-with-viewers`), and
  `projects-copy-clone.md` derives the three copy-mode variants from
  the same `demog` source. Other source classes are covered by
  parallel atlas-driven `proactive_lifecycle_specs` (one per
  source_class × dep_op cell), not by this scenario.
- **UI coverage owned (rev 3 `ui_coverage_responsibility`):**
  - `context-panel-links-url-copy` — Step 3's clipboard-copy of the
    deep-link URL from the Context Panel > Links section.
  - `new-tab-open-url` — Step 4's open-new-tab + paste/navigate flow.
  These two flows are NOT covered by any other scenario in the
  Projects chain (`ui_coverage_delegated_to: null`); save / share /
  open / delete dialogs are owned by other scenarios and are NOT
  exercised here.
- **Original `order: 4`** — runs after `upload-project.md` /
  `uploading.md` (`order: 1`), `share-project.md` (`order: 2`), and
  `opening.md` (`order: 3`). Captured in
  `scenario-chains/projects.yaml` rev 3 `order_from_files`. The
  `order: 4` vs. dependency on `projects-copy-clone.md` (`order: 5`)
  apparent contradiction is recorded in `unresolved_ambiguities`
  (rev 3) — `order` is treated as advisory; the named-variant
  evidence drives the dependency graph.
- **Source-text correction history.** Per `decision-log.yaml ::
  migration_decisions` entry `mig-2026-04-29-source-text-correction`,
  the original step 2 listed a "with layout" variant; this was
  replaced with "saved with personal view customizations" mode
  terminology because the "with layout" variant had no producer in
  the Projects section. The source on disk already reflects the
  corrected wording; this migrated body preserves the corrected
  4-variant list verbatim.
- **Cross-fixture dependency complexity.** This is the only scenario
  in the chain that requires a cross-fixture composite — variants
  produced by both `projects-copy-clone.md` AND an upstream upload.
  The Automator's `beforeAll` block needs to coordinate fixture
  setup across two producers. Surfaced as a candidate for a
  `helpers.playwright.projects.buildVariantsComposite(...)`
  registry entry in the migration report.
- **URL deep-link is generic.** Step 3's "Context Panel > Links"
  surfaces the same URL format used by `projects.url-params.build-
  share-link` for any saved project. The `<variant>` distinction
  is in WHICH project is selected, not in the URL-construction or
  URL-apply path shape (which is identical across variants). The
  4-path matrix exercises that the URL-apply behavior is uniform
  across variant modes (link / clone / personal-customizations /
  original) — directly supporting the `pyramid_layer: integration`
  source-agnostic claim.
- **No right-click Share operation.** Step 3 copies the URL from
  Context Panel > Links section — this is a clipboard-copy of a
  deep-link, NOT the right-click Share dialog. Per-cycle
  Invariant 3 (atlas-aware sub_features_covered for share) does
  NOT apply here; `projects.shell.share-via-context-menu` is
  correctly OMITTED from `sub_features_covered`.
- **Existing `project-url-spec.ts`.** A spec already exists at
  `public/packages/UsageAnalysis/files/TestTrack/Projects/project-url-spec.ts`
  (Wave 1a B70 follow-up — single-variant SCOPE_REDUCTION using
  the `demog` representative source via `page.goto({BASE}/p/{nqName})`).
  Per per-cycle Invariant 2, the existing spec is **READ-ONLY for
  Migrator**; this rev-3 migration aligns the `.md` with the chain's
  rev-3 schema and the spec's existing `demog`-as-representative-
  source choice — a coherent end-state.
