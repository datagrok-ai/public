---
realizes_atlas: []
realizes: [viewers.form]
priority: p2
target_layer: manual-only
coverage_type: smoke
---

# Form tests (manual)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Form

## Save and load form layout (file)

1. Enter design mode and rearrange fields to a custom layout
2. Click `[name="icon-download"]` (Save form to file) — a JSON file is downloaded
3. Reset the form by clicking `[name="icon-list"]` and re-selecting all columns
4. Click `[name="icon-folder-open"]` (Load form from file) — select the previously saved JSON file
5. Verify the custom layout is restored

## Cleanup

1. Delete the downloaded form JSON file from the local Downloads folder, then Close all —
   the scenario creates no server-side entities.

---
{
  "order": 201,
  "datasets": ["System:DemoFiles/demog.csv"]
}
