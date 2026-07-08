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

# Complex — Save Copy mode lifecycle

Verifies the Save Copy round trip: save a project, make an
independent copy of it with Data Sync off (a static snapshot), reopen
the copy, then re-save it with Data Sync on to turn it into a live
link to the source data. Uses a single file-share source
(`demog.csv`), since Save Copy mechanics don't depend on the
underlying data source.

UI coverage delegated to `projects-ui-smoke.md`.

## Setup

1. Authenticate as test user.
2. Project name root: `save-copy-${Date.now()}`. Two derived names:
   `<root>_NoSync` and `<root>_Sync`.
3. Cleanup: delete both projects at the end.
4. **Note on canonical save trigger:** use the ribbon SAVE button
   (NOT Ctrl+S — per `feedback_no_ctrlS_for_layouts`).

## Scenarios

### Main flow — Save Copy round-trip

1. **Open `demog` from File share.** Open
   `System:DemoFiles/demog.csv`. Verify table loaded.
2. **Save baseline project with Data Sync ON.** Save Project,
   name `<root>_Original`, Data Sync **ON**, OK. Cancel auto-share.
   This is the source project for Save Copy operations.
3. **Save Copy without Data Sync.** Trigger **Save Copy** (or
   `File > Save Project As` — the variant that creates a copy).
   In the Save dialog: name `<root>_NoSync`, Data Sync **OFF**
   for all tables, click **OK**. Verify a new project is created
   (`grok.dapi.projects.list().filter('name = "<root>_NoSync"').count()`
   returns 1). Verify the original is untouched.
4. **Close and reopen the NoSync copy.** Close all views
   (`grok.shell.closeAll()`). Open `<root>_NoSync` from
   `Browse > Dashboards`. Verify the table is loaded but Data
   Sync is OFF (i.e. table is a static snapshot from save time;
   modifying source CSV would not propagate).
5. **Re-save the copy with Data Sync ON.** Trigger Save Project
   on the open `<root>_NoSync` copy. In the Save dialog: optionally
   re-name as `<root>_Sync` (or save under same name with sync
   toggle change), Data Sync **ON**, **OK**. Verify the project's
   sync configuration changed (the table is now linked to the
   underlying CSV).
6. **Verify post-resave state.** Reopen the project. Verify Data
   Sync is now ON (modify-and-sync test: change a row in the
   underlying CSV, force Data Sync refresh, verify the table
   reflects the change). Cite `data-sync-refresh-verification`
   UI flow ownership in `complex.md` for the refresh trigger
   surface.
7. **Cleanup.** Delete `<root>_Original`, `<root>_NoSync`,
   `<root>_Sync` (3 projects depending on whether re-save in
   Step 5 used a different name).

### Expected results

- Save Copy creates an independent copy without affecting the
  original.
- Save Copy without Sync produces a static snapshot.
- Re-save with Sync converts the static snapshot to a live
  link.

## Notes

- **No related bug.** Save Copy modes are not directly targeted by
  any GROK ticket. GROK-19750 (Save Copy with Link mode must not
  affect the original project's viewers) is covered separately by
  `projects-copy-clone.md`.
- **UI coverage delegated.** The Save dialog, Data Sync toggle, and
  Browse > Dashboards opens are all owned by `projects-ui-smoke.md`.
  The flow here is mostly UI-driven for the Save dialog (there's no
  JS API substitute for the Data Sync toggle), but assertions are
  made via JS API.
- **Self-cleaning.** Step 7 deletes all projects.
