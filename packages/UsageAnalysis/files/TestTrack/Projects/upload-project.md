---
feature: projects
sub_features_covered: [projects.upload, projects.api.save]
target_layer: playwright
priority: smoke
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/upload-project.md
migration_date: 2026-04-28
migration_report: upload-project-migration-report.md
related_bugs: [GROK-19750, GROK-19212, GROK-19103]
---

# Upload Project

Save a dataset with multiple viewers as a new Datagrok project, with
the **Data Sync** option enabled, then dismiss the post-save Share
dialog without sharing. Golden-path smoke for the Projects section.

## Setup

1. Open the demo dataset `System:DemoFiles/demog.csv`.
2. Add a **Scatter Plot** viewer to the table view.
3. Add a **Bar Chart** viewer to the table view.
4. Add a **Line Chart** viewer to the table view.

## Scenarios

### Save with Data Sync on

Saves the dataset and its three viewers as a new project with
`Data Sync` enabled, then cancels the follow-up Share dialog without
sharing.

1. Trigger **Save Project** (main menu or the standard project-save
   shortcut).
2. In the **Save Project** dialog, leave **Data Sync** toggled ON.
3. Click **OK** and wait for the save to complete (verify the project
   exists on the server before continuing — e.g. via the project's
   appearance in **Browse > Dashboards** or via the underlying
   `POST /projects` succeeding).
4. The **Share demog** dialog opens automatically; click **Cancel** to
   dismiss it (verify the dialog closes without sharing).
5. **Close All** open views and viewers (cleanup).

## Notes

- Dataset path source: original scenario's trailing JSON metadata
  (`datasets: ["System:DemoFiles/demog.csv"]`).
- Original `order: 1` — this scenario produces the saved `demog`
  project consumed by `share-project.md`, `opening.md`,
  `projects-copy-clone.md`, `complex.md`, and others. Surfaced as a
  candidate fixture (`demog-project-with-viewers`) for Step 1 chain
  analysis to consider when the section is migrated as a batch.
- The `uploadProject(projectName, tableInfo, view, df)` helper in
  `helpers-registry.yaml` is `grok_test_layer` (a/b style) only; it
  is NOT Playwright-compatible. At the playwright layer Automator
  must drive the UI dialog directly. Flagged in the migration report
  as a candidate for a Playwright-layer counterpart.
