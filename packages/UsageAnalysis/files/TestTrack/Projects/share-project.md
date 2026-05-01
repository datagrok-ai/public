---
feature: projects
sub_features_covered: [projects.api.search]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/share-project.md
migration_date: 2026-04-30
migration_report: share-project-migration-report.md
ui_companion: share-project-ui.md
related_bugs: []
---

# Share Project

Share the saved `demog` project (produced by `upload-project.md`) with
a registered user via the JS-API-substitutable path
(`grok.dapi.permissions.grant`).

The unregistered-email-invite path (auto-creates an account on the
server), the Sharing-tab listing verification, and the post-share
Context-Panel render verification are NOT in this scenario — they
live in the `share-project-ui.md` companion file as UI-only manual
steps (per D3 bucket-b classification — auto-create-account behavior
and Context Panel render quality require human QA judgment).

This scenario depends on `upload-project.md` having produced a saved
project named `demog` (produces `demog-project-with-viewers` fixture
per scenario-chains rev 2).

## Setup

1. **Fixture prerequisite:** the saved project `demog` (with three
   viewers and Data Sync ON) must exist on the server. This is
   produced by `upload-project.md`. In `end_to_end_fixtures` strategy,
   the Automator stage will set this up via a `beforeAll` block using
   `js-api-replay` of `upload-project.md`'s save flow OR by reusing
   an in-process fixture builder.
2. **Recipient identity for the registered-user share:** `Olena Ahadzhanian`
   (literal account name, as in the original).

## Scenarios

### Share demog with a registered user

1. Go to **Browse** > **Dashboards**.
2. Use the search bar to find the `demog` project.
3. Right-click the project and select **Share**.
4. In the Share dialog, share the project with a registered user —
   `Olena Ahadzhanian`.

> **UI-only steps moved to `share-project-ui.md`** (preserving original
> numbering for cross-reference):
> - Step 4 sub-bullet: "Via email (to an unregistered email recipient)"
> - Step 5: Verify on Context Panel — Sharing tab (auto-created user listed)
> - Step 7: Right-click → Details
> - Step 8: Review every Context Panel tab

## Notes

- **Cross-feature touch:** the actual Share operation (steps 3-4) is
  driven by the platform's permissions / sharing system, not by the
  `projects` feature directly. The migrated frontmatter
  `sub_features_covered` reflects the projects-side touch point
  exercised in this scenario only (`projects.api.search`). Wave 4
  cleanup removed `projects.view.browse` / `projects.view.single`
  overclaims.
- **Hard dependency on `upload-project.md`:** the `demog` project name
  in step 2 is produced ONLY by `upload-project.md`. `uploading.md`
  produces `Test_Case<N>_Sync` / `Test_Case<N>_NoSync` names; those
  are not consumed here.
- **Cleanup responsibility:** the auto-created user account from the
  email invite (UI companion file path) must be cleaned up after the
  manual run — see `share-project-ui.md` cleanup section.
- Original `order: 2` — runs after `upload-project.md` (`order: 1`).
  Captured in `scenario-chains/projects.yaml` rev 2 `dependency_graph`
  (`share-project.md.depends_on: [upload-project.md]`).
