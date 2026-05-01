---
feature: projects
sub_features_covered: [projects.url-params.build-share-link, projects.url-params.apply, projects.shell.open]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/project-url.md
migration_date: 2026-04-30
migration_report: project-url-migration-report.md
related_bugs: []
---

# Project URL

For each of four project variants (original; copied-with-link;
copied-with-clone; saved-with-personal-view-customizations), navigate
to Browse > Dashboards, locate the project, copy its deep-link URL
from Context Panel > Links, then open the URL in a new browser tab
and verify the corresponding project loads.

This scenario depends on a composite fixture `copy-clone-customizations-variants`
produced by `projects-copy-clone.md` (variants: link / clone / personal-view-
customizations) plus an "original" project from either `upload-project.md`
(saved demog) or `uploading.md` (any of `Test_Case<N>_Sync` / `_NoSync`).

## Setup

1. **Fixture prerequisite:** the four project variants must exist on
   the server:
   - **original**: a saved project (e.g. `demog` from `upload-project.md`,
     OR any `Test_Case<N>_Sync` from `uploading.md`).
   - **copied-with-link**: a Save-Copy-with-Link variant produced by
     `projects-copy-clone.md`.
   - **copied-with-clone**: a Save-Copy-with-Clone variant produced by
     `projects-copy-clone.md`.
   - **saved-with-personal-view-customizations**: a Save-personal-view-
     customizations variant produced by `projects-copy-clone.md`.
   In `end_to_end_fixtures` strategy, the Automator stage will build
   the composite `copy-clone-customizations-variants` fixture in a
   `beforeAll` block by either reusing fixtures from prior chain runs
   or replaying the relevant save-flow cases via `js-api-replay`.
2. The browser session must be authenticated as a user who has access
   to all four variants (typically the project owner; cross-user share
   verification is OUT of scope for this scenario — `share-project.md`
   covers that).

## Scenarios

### URL deep-link reopen for each variant

The following workflow is parameterised over the 4 project variants.
Implemented per `end_to_end_fixtures` strategy at the Playwright layer
(per `output_plan.project-url.md.strategy = end_to_end_fixtures` in
`scenario-chains/projects.yaml` rev 2).

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

- **Original `order: 4`** — runs after `upload-project.md` /
  `uploading.md` (`order: 1`), `share-project.md` (`order: 2`), and
  `opening.md` (`order: 3`). Captured in
  `scenario-chains/projects.yaml` rev 2 `order_from_files`.
- **Source-text correction history.** Per
  `decision-log.yaml :: migration_decisions` entry
  `mig-2026-04-29-source-text-correction`, the original step 2 listed
  a "with layout" variant; this was replaced with "saved with personal
  view customizations" mode terminology because the "with layout"
  variant had no producer in the Projects section. The source on disk
  already reflects the corrected wording; this migrated body preserves
  the corrected 4-variant list verbatim.
- **Cross-fixture dependency complexity.** This is the only scenario
  in the chain that requires a cross-fixture composite — variants
  produced by both `projects-copy-clone.md` AND an upstream upload.
  The Automator's `beforeAll` block needs to coordinate fixture
  setup across two producers. Surfaced as a candidate for a
  `helpers.playwright.projects.buildVariantsComposite(...)` registry
  entry in the migration report.
- **URL deep-link is generic.** Step 3's "Context Panel > Links"
  surfaces the same URL format used by `projects.url-params.build-share-link`
  for any saved project. The `<variant>` distinction is in WHICH
  project is selected, not in the URL-construction or URL-apply path
  shape (which is identical across variants). The 4-path matrix
  exercises that the URL-apply behavior is uniform across variant
  modes (link / clone / personal-customizations / original).
- **No right-click Share operation.** Step 3 copies the URL from
  Context Panel > Links section — this is a clipboard-copy of a
  deep-link, NOT the right-click Share dialog. Per-cycle Invariant 3
  (atlas-aware sub_features_covered for share) does NOT apply here;
  `projects.shell.share-via-context-menu` is correctly OMITTED from
  `sub_features_covered`.
- **No existing `project-url-spec.ts`.** Per-cycle Invariant 2
  (existing -spec.ts / -api.ts READ-ONLY) is trivially satisfied —
  no project-url spec file exists at any path under `public/`. The
  Automator stage will create one fresh during a future cycle.
