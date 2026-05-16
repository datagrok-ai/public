---
feature: projects
ui_only: true
target_layer: manual
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/custom-creation-scripts.md
migration_date: 2026-04-30
moved_to_ui_at: 2026-05-01
moved_to_ui_reason: |
  All 3 in-scope scenario verifications (Step 1 run creation script,
  Step 4 mutate System:DemoFiles/chem, Step 5 verify dynamic-data-
  resolution on reopen) had only environmental/infra blockers per D3
  classification (registered server function, per-test scratch folder,
  qa-pw write perms on shared System: paths). Per 2026-05-01
  decision: rather than tracking these as bucket-a gaps and waiting
  for env/helper backlog to close, mark the entire scenario as
  UI-only manual — the dynamic-data-resolution behavior is most
  efficiently verified by human QA who can set up the chem folder
  mutation, run the script via Tools > Scripting, and inspect the
  result on reopen. The Playwright spec (custom-creation-scripts-spec.ts)
  + verdict file (custom-creation-scripts-spec.verdict.yaml) were
  deleted alongside this rename.
related_bugs: []
---

# Custom Creation Scripts — UI-only manual scenario

Verify that a project bound to a JS creation-script as its data source
correctly re-executes the script on reopen and resolves the latest
data — specifically, that the script picks up file-storage changes
made between save and reopen. The creation script in this scenario
selects the CSV file with the highest numeric suffix from
`System:DemoFiles/chem`; the test mutates that folder between save
and reopen, then verifies the reopened project loads the new highest-
suffixed file.

This scenario is **UI-only manual** — see frontmatter `moved_to_ui_reason`
for rationale. There is no Playwright spec for this scenario. Human QA
runs this scenario manually when validating dynamic-data-resolution
behavior.

This scenario is **self-contained** per
`scenario-chains/projects.yaml` rev 2 (`depends_on: []`,
classification `medium`); no upstream fixture consumed.

The scenario contains a **destructive shared-state mutation** in
step 4 (modifies `System:DemoFiles/chem`, which is a shared file-
storage path). Manual QA must coordinate with concurrent test runs
OR perform the mutation in a controlled time window, AND restore the
folder afterwards (see Cleanup section below).

## Setup

1. **Filesystem prerequisite:** the path `System:DemoFiles/chem`
   must exist on the test server and contain at least one CSV file
   with a numeric suffix in its filename (e.g. `chem_data_001.csv`,
   `compounds_v42.csv`). The creation script in step 1 throws if
   the folder is empty.
2. **Destructive shared-state caveat:** step 4 mutates
   `System:DemoFiles/chem`. Concurrent test runs OR other tests
   reading the same folder will see the mutation. Coordinate timing
   with team OR perform the mutation in a controlled window. Restore
   the folder to a known state in Cleanup.
3. **Browser session authentication:** test-user session with
   read+write access to `System:DemoFiles/chem` and project-creation
   rights (qa-pw on dev currently lacks write access — manual run
   may require elevated credentials).

## Manual scenario

### Creation-script-on-reopen rebuild verification

1. **Run the creation script** (inline below). The script lists
   CSV files in `System:DemoFiles/chem`, sorts by numeric suffix in
   the filename, selects the file with the highest suffix, and
   produces a dataframe from its contents. Run via the platform's
   script-execution UI (e.g. `Tools` > `Scripting` > paste +
   Run, OR `grok.functions.eval(...)` in the browser console for a
   programmatic run).

   ```js
   //name: getLastCreatedFile_by_numeric_suffix
   //language: javascript
   //output: dataframe result

   const csvFiles = await grok.dapi.files.list('System:DemoFiles/chem', true, 'csv');

   if (csvFiles.length === 0)
       throw new Error('No CSV files found in System:DemoFiles/chem');

   // Function to extract a numeric suffix from the file name
   function getNumberSuffix(name) {
       const match = name.match(/(\d+)(?=\.csv$)/);
       return match ? parseInt(match[1], 10) : -1; // -1 if no number is found
   }

   // Sort files by numeric suffix in ascending order
   csvFiles.sort((a, b) => getNumberSuffix(a.fileName) - getNumberSuffix(b.fileName));

   // Select the file with the highest numeric suffix
   const lastFile = csvFiles[csvFiles.length - 1];

   const csv = await grok.dapi.files.readAsText(lastFile.fullPath);
   result = DG.DataFrame.fromCsv(csv);
   ```

   Verify the script produces a dataframe and the result table
   opens in the workspace.

2. **Add some viewers and save the project with Data Sync enabled.**
   On the dataframe from step 1, add 1+ viewers (e.g. a scatter
   plot, bar chart, or any viewer of choice) to the table view.
   Trigger **File** > **Save Project**. In the Save Project dialog,
   ensure the **Data Sync** toggle is **ON** (this binds the script
   as the project's creation script — the script will re-run on
   project open). Name the project (e.g.
   `Custom_Creation_Script_Test`) and click **OK**.
3. **Close All.** Close all open views and viewers. Verify the
   workspace is clear.
4. **Update any CSV file in the `System:DemoFiles/chem` folder.**
   Either:
   - **Add a new CSV file with a HIGHER numeric suffix** than any
     existing file (e.g. if highest existing suffix is 042, add
     `chem_data_999.csv`), OR
   - **Rename an existing file to add or increase its numeric
     suffix** (e.g. rename `compounds.csv` →
     `compounds_999.csv`).
   Either action makes a different file the "highest-suffixed",
   so the creation script's selection should change on reopen.
   *(WARNING: this is a destructive mutation of shared filesystem
   state. Restore in Cleanup.)*
5. **Open the saved project and verify on reopen:**
   - **No errors occur** — the project opens without console
     errors (F12), and no error balloon appears.
   - **The most recently created or modified file in the `chem`
     folder is loaded** — the dataframe in the reopened project
     reflects the contents of the file added or renamed in step 4
     (i.e. the file with the highest numeric suffix at reopen
     time, NOT the file that was highest at save time).

## Cleanup responsibility

After the manual run, restore environment to keep it idempotent:

1. **Delete the saved project** `Custom_Creation_Script_Test` via UI
   (Browse > Dashboards > right-click > Delete) or CLI
   (`grok s projects delete Custom_Creation_Script_Test`).
2. **Restore `System:DemoFiles/chem`** to its pre-mutation state:
   - If a new file was added in step 4, delete it.
   - If a file was renamed in step 4, rename it back.
3. **Log the cleanup** in `custom-creation-scripts-run.md` under a
   new dated entry (preserve the run record + cleanup confirmation).

## Sign-off

PASS — the project's dataframe on reopen reflects the post-mutation
highest-suffixed file. No console errors. Folder restored.
FAIL — list deviation (e.g., dataframe still reflects pre-mutation
file, OR error balloon on reopen, OR script execution failure) +
screenshot. Log in `custom-creation-scripts-run.md` under a new
dated entry.

## Notes

- **Atlas critical_path partial coverage.** This scenario covers the
  creation-script-rebuild-on-reopen behavior; URL-param + share
  aspects of `url-parameterized-share` critical_path are NOT
  exercised here. Atlas may benefit from a dedicated
  `projects.creation-script.*` sub_feature for non-URL-param
  creation-script flows (flagged in migration report Section 1
  for review).
- **No spec / no automation tracking.** Per 2026-05-01 decision,
  this scenario is fully manual. No corresponding `-spec.ts` exists.
  No bucket-a coverage gaps tracked in `coverage-gaps/projects.md`
  for this scenario (excluded — see "Test cases marked UI-only"
  section in coverage-gaps/projects.md).
