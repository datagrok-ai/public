---
feature: projects
sub_features_covered:
  - projects.api.save
  - projects.api.namespaces
  - projects.move.move-entity
  - projects.move.commit
target_layer: playwright
coverage_type: regression
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
related_bugs: []
---

# Complex — Move project lifecycle

Source-agnostic op (Wave 2B). Decomposition of `complex.md` Step 10:
Move project → file share → Space. Single representative source
(Files + `demog.csv`).

**Important — JS API path is the primary mechanism.** Per decision-log
`b2-2026-05-03-drag-drop-ui-only-reclassification` (2026-05-03), the
two documented UI paths are blocked:
1. **Drag-drop in Browse tree** — drag-handler registration not
   accessible from Playwright (3 mechanisms tried: synthetic
   DragEvent, raw mouse events, CDP Input.dispatchMouseEvent).
2. **Right-click `Move to` menu option** — does NOT exist in current
   UI (verified 2026-05-03 against dev.datagrok.ai).

This scenario therefore uses **`grok.dapi.projects` namespace setters**
as primary path per Plan line 168. JS API is more reliable and
cleaner than drag-drop. UI coverage for any move-related context
menu surface (if/when it lands) is owned by `projects-ui-smoke.md`.

## Setup

1. Authenticate as test user.
2. Project name: `move-test-${Date.now()}`.
3. **Environment dependencies:**
   - File share namespace accessible to test user (e.g.
     `System:AppData/MyShare/`).
   - Space accessible to test user (e.g. inline-created via JS
     API per Spaces prelude pattern).
4. Cleanup: delete the project; delete the inline-created Space.

## Scenarios

### Main flow — move project across namespaces

1. **Open `demog` from File share.** Open
   `System:DemoFiles/demog.csv`. Verify table loaded.
2. **Save project with Data Sync ON.** Save Project, name from
   Setup, Data Sync **ON**, OK. Cancel auto-share. Verify
   project's initial namespace via
   `(await grok.dapi.projects.find(<id>)).fullName` — should
   be in the test user's default namespace.
3. **Move project to a file-share namespace via JS API.**
   ```js
   const fileShareNamespace = '<test-user-namespace>'; // resolve via grok.dapi.namespaces
   await grok.dapi.projects.move(project, fileShareNamespace);
   ```
   - Verify the project's new namespace via
     `(await grok.dapi.projects.find(<id>)).fullName`.
   - Verify the project still opens (re-find by id; load tables;
     no missing-data errors).
4. **Move project to a Space (inline-created).**
   ```js
   const space = await grok.dapi.spaces.createRoot('move-target-${Date.now()}');
   const spaceNamespace = `Spaces:${space.name}`;
   await grok.dapi.projects.move(project, spaceNamespace);
   ```
   - Verify project's new namespace.
   - Verify project still opens; tables still loaded; relations
     preserved.
5. **Cleanup.** Delete the project. Delete the inline Space.

### Expected results

- Move via JS API succeeds in both directions (default namespace
  → file share namespace → Space namespace).
- Project remains openable post-move.
- Tables and relations preserved across moves.

## Notes

- **Origin: `complex.md` Step 10 split (Wave 2B per Plan line
  168).** Step 10 is source-agnostic per Plan; one
  representative source (Files) covers the path for all sources.
- **JS API primary path.** Per Plan line 168 + decision-log
  `b2-2026-05-03-drag-drop-ui-only-reclassification`. The atlas
  `manual_only` block lists `projects.move.move-entity`,
  `projects.move.commit`, `projects.move.move-relations` as
  manual-only (UI-only) — those refer to the UI surface that
  doesn't exist. The JS API path
  (`grok.dapi.projects.move(...)`) IS testable and is what this
  scenario covers.
- **UI surface deferred / atlas-side.** If/when a right-click
  `Move to` UI is added, atlas `manual_only` flagging would be
  re-evaluated and a UI step could be added to
  `projects-ui-smoke.md`. Until then: JS API only.
- **No `related_bugs`.** Move ops have no GROK ticket. Pure
  proactive coverage of the move-namespace path.
- **`projects.move.*` sub_features.** Atlas v11 lists these
  under `manual_only`; this scenario tests the JS API equivalent
  contract (`grok.dapi.projects.move`). Surfacing for atlas:
  the `manual_only` rationale should be refined to "UI-only
  paths are manual; JS API path is testable" for these sub-
  features.
- **Self-cleaning.** Step 5 deletes everything created.
- **Sequencing within Wave 2B.** Second after save-copy per
  Plan execution order — no helper / env blockers beyond
  namespace + Space access.
