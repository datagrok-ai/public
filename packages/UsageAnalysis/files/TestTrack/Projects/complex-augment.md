---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p2
realizes: []
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
related_bugs: []
---

# Complex — Augment project (drag-drop add tables)

Verifies that tables can be added to an already-saved project
(augmenting it) using files from `System:DemoFiles` as the source.
Because the drag-and-drop UI for this can't be driven through
Playwright, this scenario exercises the same augmentation through the
JS API (`grok.dapi.projects.addRelation`) instead, and documents the
UI path as a known manual-only gap.

## Setup

1. Authenticate as test user.
2. Project name: `augment-test-${Date.now()}`.
3. Source tables (4 files from `System:DemoFiles`):
   `demog.csv`, `cars.csv`, `iris.csv`, `geo.csv` (or any 4
   reachable .csv files in DemoFiles).
4. Cleanup: delete the project at the end.

## Scenarios

### Main flow — augment via JS API

1. **Open the first table (`demog.csv`) and save baseline
   project.** Open via Files browser. Save Project, name from
   Setup, Data Sync **ON**, OK. Cancel auto-share. Verify
   single-table project saved.
2. **Augment project: add 3 more tables via JS API
   (`addRelation`, Link mode).**
   ```js
   const project = await grok.dapi.projects.find(<id>);
   for (const fileName of ['cars.csv', 'iris.csv', 'geo.csv']) {
     const file = (await grok.dapi.files.list('System:DemoFiles', false, fileName))[0];
     await grok.dapi.projects.addRelation(project, file, /*linkMode=*/true);
   }
   await grok.dapi.projects.save(project);
   ```
   - Verify the project now has 4 relations:
     `(await grok.dapi.projects.find(<id>)).relations.length === 4`.
3. **Re-open the augmented project and verify all 4 tables
   load.** Close all views; reopen project from Browse >
   Dashboards. Verify `grok.shell.tables.length === 4`. Verify
   each table is the expected source file.
4. **(Optional UI drag-drop verification — if/when
   automatable.)** If the drag-drop UI mechanism becomes
   automatable in the future, add an additional step that
   exercises Browse > Files > drag-drop onto an open Dashboard.
   Currently: UI path is documented as a known gap; this
   scenario does NOT exercise the UI drag-drop path.
5. **Cleanup.** Delete the project.

### Expected results

- `addRelation` JS API correctly adds source tables to an
  existing project in Link mode.
- All 4 tables persist across save → close → reopen.
- The project's relations list reflects the 4 sources.

## Notes

- **JS API path is primary; UI drag-drop is a known gap.** The
  drag-drop UI surface in the Browse tree could not be made to work
  with Playwright after trying three different simulated-event
  mechanisms, so this scenario exercises the equivalent
  `addRelation` JS API call instead. If the drag-drop mechanism
  becomes automatable in the future, a UI verification step should
  be added.
- **No related bug.** This is proactive coverage of project
  augmentation; no GROK ticket targets it directly.
- **UI coverage delegated.** No UI surface is exercised here beyond
  Save Project, which is owned by `projects-ui-smoke.md`.
- **Self-cleaning.** Step 5 deletes the project.
