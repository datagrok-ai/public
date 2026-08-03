---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: []
realizes: []
realized_as:
  - complex-move-spec.ts
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
related_bugs: []
---

# Complex — Move project lifecycle

Verifies that a saved project can be moved between namespaces — from
its default location to a file-share namespace, then to a Space — and
that it still opens correctly with all its tables and relations
intact after each move. Uses a single file-share source
(`demog.csv`).

Neither UI path for moving a project currently works: drag-and-drop
in the Browse tree can't be driven through Playwright, and the
right-click "Move to" menu option does not exist in the current UI
(verified against dev.datagrok.ai). So this scenario drives the move
through the `grok.dapi.projects` JS API instead. UI coverage for any
move-related context-menu surface, if and when one is added, is owned
by `projects-ui-smoke.md`.

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

- **JS API path is primary; UI move is not currently available.**
  The move-related sub-features are marked UI-only in the feature
  atlas, but that classification refers to the missing UI surface
  (no working drag-drop, no "Move to" menu item) — the JS API move
  contract (`grok.dapi.projects.move(...)`) itself is fully testable
  and is what this scenario covers. If a right-click "Move to" UI is
  added in the future, a corresponding step could be added to
  `projects-ui-smoke.md`.
- **No related bug.** Move operations have no GROK ticket; this is
  proactive coverage of the move-namespace path.
- **Self-cleaning.** Step 5 deletes everything created.
