---
feature: projects
sub_features_covered:
  - projects.api.save
  - projects.add-relation
  - projects.tree.add-entity-node
target_layer: playwright
coverage_type: regression
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
migration_report: complex-augment-migration-report.md
related_bugs: []
---

# Complex — Augment project (drag-drop add tables)

Source-agnostic op (Wave 2B). Decomposition of `complex.md` Step 4:
drag-drop 4 tables onto Dashboards (Link mode) to augment an existing
project. Single representative source (Files + multiple `.csv` files
from `System:DemoFiles`).

**Important — JS API fallback is primary.** Per Plan line 169 +
decision-log `b2-2026-05-03-drag-drop-ui-only-reclassification`, the
drag-drop UI surface in Browse tree is NOT automatable from Playwright
(3 mechanisms tried). This scenario uses
`grok.dapi.projects.addRelation` as primary path. UI drag-drop is
documented as a known gap.

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

- **Origin: `complex.md` Step 4 split (Wave 2B per Plan line
  169).** Step 4 is source-agnostic — augment via drag-drop
  works the same regardless of source class.
- **JS API primary path; UI drag-drop deferred.** Per Plan
  line 169 — Plan explicitly says "Plan JS-API fallback
  `grok.dapi.projects.addRelation` from the start" precisely
  because drag-drop is flagged unreliable in atlas
  `tree-drag-drop-actions` AND not automatable per
  `b2-2026-05-03-drag-drop-ui-only-reclassification`.
- **`projects.tree.add-entity-node` sub_feature in atlas
  `manual_only`.** This scenario tests the JS API contract
  equivalent (`addRelation`). Surfacing for atlas: same as
  `complex-move.md` — `manual_only` rationale should be
  refined to "UI drag-drop is manual; JS API addRelation is
  testable".
- **No `related_bugs`.** Drag-drop augment has no GROK
  ticket. Proactive coverage of project augmentation path.
- **UI coverage delegated.** No UI surface exercised in this
  scenario beyond Save Project (owned by
  `projects-ui-smoke.md`).
- **Self-cleaning.** Step 5 deletes the project.
- **Sequencing within Wave 2B.** Third after save-copy + move
  per Plan execution order.
