---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [upload-save-reopen-golden, derive-then-save-inside-project]
realizes: []
pyramid_layer: source-matrix
ui_coverage_responsibility:
  - save-project-dialog
  - link-tables-dialog
  - join-tables-dialog
  - pivot-table-add-to-workspace
  - aggregate-rows-add-to-workspace
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/uploading.md
migration_date: 2026-05-04
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - selection-to-filter-linking-direction-asymmetry
  - spgi-v2-infinity-csv-source-ambiguity
  - project-naming-collision-risk-in-ci
  - add-to-workspace-semantics
  - original-step-renumbering-across-sync-on-sync-off-blocks
  - order-1-collision-with-upload-project-md
scope_reductions: []
related_bugs:
  - GROK-19103
  - GROK-18345
---

# Uploading — Multi-source project save matrix

Save 8 distinct two-table-source combinations as Datagrok projects,
each with both **Data Sync ON** and **Data Sync OFF** variants — **16
saved projects total**. Verifies that project upload, save, reopen,
and inter-table linking/joining work across all major data sources:
local files, database query results, Spaces, file shares,
pivot-table derivations, and aggregate-row derivations.

DB-table right-click flows (`Get Top 100`, `Get All`) are covered
by their dedicated scenario `../Queries/get-all-get-top-100.md`;
this scenario uses plain DB-table double-click for DB sources.

Save Project dialog and Data Sync toggle UI coverage is owned by
`projects-ui-smoke.md`; this scenario adds the source-specific dialog
flows — Link Tables, Join Tables, Pivot Table > Add, Aggregate Rows
> Add — that don't appear there.

## Setup

A clean Datagrok session is the only shared setup. Each case opens
its own pair of tables, performs its own link/join, then saves and
reopens. Cleanup between cases is the responsibility of the test
harness (close all views before the next case).

**Spaces prelude (Cases 4–6 only):** before running any Spaces-based
case, create a Space named `test-projects-demo` and add `demog.csv`
from `System:DemoFiles` to it via JS API (per decision-log
`sa-2026-05-03-spaces-inline-prelude-pattern`):

```js
const space = await grok.dapi.spaces.createRoot('test-projects-demo');
const client = await grok.dapi.spaces.spaceClient(space.id);
const file = (await grok.dapi.files.list('System:DemoFiles', false, 'demog.csv'))[0];
await client.addEntity(file, /*link=*/true);
```

**Spaces postlude (after Cases 4–6):** delete the Space:

```js
await grok.dapi.spaces.delete(space);
```

## Scenarios

### Matrix: Table source combinations

| Case | Table 1 source                    | Table 2 source                        |
|------|-----------------------------------|---------------------------------------|
| 1    | Browse > Files                    | Browse > Files                        |
| 2    | Query result (Run / Double-click) | Query result (Run / Double-click)     |
| 3    | Query result (Run / Double-click) | Browse > Files                        |
| 4    | Spaces                            | Spaces                                |
| 5    | Spaces                            | Browse > Files                        |
| 6    | Spaces                            | Query result (Run / Double-click)     |
| 8    | Browse > Files                    | Pivot Table > Add to workspace        |
| 9    | DB table (double-click)           | Aggregate Rows > Add to workspace     |

Case numbering preserves 1–6, 8, 9 (no Case 7) so existing spec
mappings (`uploading-spec.ts` covers Cases 1, 3, 4, 9) and run
artefacts stay valid. Case 7 (Get Top 100 / Get All + DB table
double-click + Join Tables) is owned by
`../Queries/get-all-get-top-100.md` and is not duplicated here.

Each case below is run twice — once with **Data Sync ON** (saved as
`Test_Case<N>_Sync`) and once with **Data Sync OFF** (saved as
`Test_Case<N>_NoSync`). Both variants are followed by a
close-and-reopen verification. **Total: 8 cases × 2 sync states = 16
saved projects.**

Per chain `ui_coverage_plan` rev 3 — Save Project dialog and Data
Sync toggle smoke surface is owned by `projects-ui-smoke.md`; this
scenario witnesses the source-specific dialog flows (Link Tables,
Join Tables, Pivot Table > Add to workspace, Aggregate Rows > Add to
workspace).

---

### Case 1: Browse > Files + Browse > Files (Link Tables)

1. Open **Browse** > **Files** > `System:AppData` > `Chem` > `tests`.
2. Double-click `spgi-100.csv` — wait for the table to open.
3. Go back to **Browse** > **Files** > `System:AppData` > `Chem` > `tests`.
4. Double-click `spgi-100.csv` again — wait for the second table tab
   to open (Datagrok auto-disambiguates names: `spgi-100` and
   `spgi-100 (2)`).
5. Go to **Data** > **Link Tables**:
   - Table 1: `spgi-100`, Table 2: `spgi-100 (2)`
   - Link column: `ID`
   - Link type: **Selection to filter**
   - Click **OK**
6. **Verify:** select rows in `spgi-100` — corresponding rows
   in `spgi-100 (2)` are filtered.
7. **Save with Data Sync ON:**
   - **File** > **Save Project** (or Ctrl+S).
   - In the Save Project dialog, ensure the **Data sync** toggle is
     **ON** for both tables.
   - Name the project `Test_Case1_Sync` and click **OK**.
8. Close the project. Reopen it from **Browse** > **Projects**.
9. **Verify on reopen:** both tables are loaded, linking still works
   (selection in one filters the other), no errors in the browser
   console (F12).
10. **Save with Data Sync OFF:**
    - Repeat steps 1–6.
    - In the Save Project dialog, turn the **Data sync** toggle **OFF**
      for both tables.
    - Save as `Test_Case1_NoSync`.
11. Close the project. Reopen it from **Browse** > **Projects**.
12. **Verify on reopen:** both tables are loaded, linking works,
    no console errors.

---

### Case 2: Query result + Query result (Link Tables)

> **Prerequisite:** spec provisions a saved query on `System:Datagrok`
> via `helpers/openers.ts:provisionSystemDatagrokQuery({sql:
> SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE})` before the case starts;
> deletes it via `provisioned.cleanup()` after.

1. Run the provisioned query via
   `helpers/openers.ts:openTableFromDbQuery(page,
   provisioned.queryNqName)` — wait for the result table to open.
2. Run the same query again to get a second result table in the
   workspace.
3. Go to **Data** > **Link Tables**:
   - Table 1: first query result, Table 2: second query result
   - Link column: `id` (always present in `System:Datagrok` metadata
     tables)
   - Link type: **Selection to filter**
   - Click **OK**
4. **Verify:** select rows in Table 1 — corresponding rows in Table 2
   are filtered.
5. **Save with Data Sync ON:**
   - **File** > **Save Project**.
   - Ensure the **Data sync** toggle is **ON** for both tables (the
     CREATION SCRIPT block should be visible).
   - Save as `Test_Case2_Sync`, click **OK**.
6. Close the project. Reopen from **Browse** > **Projects**.
7. **Verify on reopen:** both tables are loaded, linking works,
   no console errors.
8. **Save with Data Sync OFF:**
    - Repeat steps 1–4.
    - Turn the **Data sync** toggle **OFF** for both tables.
    - Save as `Test_Case2_NoSync`.
9. Close and reopen.
10. **Verify on reopen:** tables loaded, linking works, no console
    errors.

---

### Case 3: Query result + Browse > Files (Link Tables)

> **Prerequisite:** spec provisions a saved query on `System:Datagrok`
> (same pattern as Case 2). Deletes it after the case ends.

1. Run the provisioned query via
   `helpers/openers.ts:openTableFromDbQuery(page,
   provisioned.queryNqName)` — wait for the result table to open.
2. Open **Browse** > **Files** > `System:AppData` > `Chem` > `tests`
   and double-click `spgi-100.csv` (or use
   `helpers/openers.ts:openTableFromFile`) — wait for the table to
   open.
3. Go to **Data** > **Link Tables**:
   - Table 1: query result, Table 2: `spgi-100`
   - Link column: any shared key column
   - Link type: **Selection to filter**
   - Click **OK**
4. **Verify:** linking works.
5. **Save with Data Sync ON:**
   - **File** > **Save Project**.
   - **Data sync** toggle **ON** for both tables.
   - Save as `Test_Case3_Sync`.
6. Close and reopen.
7. **Verify on reopen:** tables loaded, linking works, no console
   errors.
8. **Save with Data Sync OFF:**
    - Repeat steps 1–4, **Data sync** toggle **OFF** for both tables.
    - Save as `Test_Case3_NoSync`.
9. Close and reopen.
10. **Verify on reopen:** tables loaded, linking works, no console
    errors.

---

### Case 4: Spaces + Spaces (Link Tables)

> **Prelude:** ensure the `test-projects-demo` Space exists with
> `demog.csv` added (see Setup section above).

1. Open the **test-projects-demo** space: navigate to **Browse** >
   **Spaces** > **test-projects-demo**.
2. Double-click `demog.csv` — wait for the table to open.
3. Go back to the **test-projects-demo** space.
4. Double-click `demog.csv` again — wait for the second table tab
   to open (auto-disambiguated as `demog (2)`).
5. Go to **Data** > **Link Tables**:
   - Table 1: `demog`, Table 2: `demog (2)`
   - Link column: any shared key in `demog.csv`
   - Link type: **Selection to filter**
   - Click **OK**
6. **Verify:** linking works.
7. **Save with Data Sync ON:**
   - **File** > **Save Project**.
   - **Data sync** toggle **ON** for both tables.
   - Save as `Test_Case4_Sync`.
8. Close and reopen.
9. **Verify on reopen:** tables loaded, linking works, no console
   errors.
10. **Save with Data Sync OFF:**
    - Repeat steps 1–6, **Data sync** toggle **OFF** for both tables.
    - Save as `Test_Case4_NoSync`.
11. Close and reopen.
12. **Verify on reopen:** tables loaded, linking works, no console
    errors.

---

### Case 5: Spaces + Browse > Files (Link Tables)

> **Prelude:** ensure the `test-projects-demo` Space exists with
> `demog.csv` added (see Setup section above).

1. Open the **test-projects-demo** space and double-click `demog.csv`.
2. Open **Browse** > **Files** > `System:AppData` > `Chem` > `tests`
   and double-click `spgi-100.csv`.
3. Go to **Data** > **Link Tables**:
   - Table 1: `demog`, Table 2: `spgi-100`
   - Link column: any shared key column
   - Link type: **Selection to filter**
   - Click **OK**
4. **Verify:** linking works.
5. **Save with Data Sync ON:**
   - **Data sync** toggle **ON** for both tables.
   - Save as `Test_Case5_Sync`.
6. Close and reopen.
7. **Verify on reopen:** tables loaded, linking works, no console
   errors.
8. **Save with Data Sync OFF:**
   - Repeat steps 1–4, **Data sync OFF** for both tables.
   - Save as `Test_Case5_NoSync`.
9. Close and reopen.
10. **Verify on reopen:** tables loaded, linking works, no console
    errors.

---

### Case 6: Spaces + Query result (Link Tables)

> **Prelude:** ensure the `test-projects-demo` Space exists with
> `demog.csv` added (see Setup section above).
>
> **Prerequisite:** spec provisions a saved query on `System:Datagrok`
> via `helpers/openers.ts:provisionSystemDatagrokQuery`. Deletes it
> after the case ends.

1. Open the **test-projects-demo** space and double-click `demog.csv`.
2. Run the provisioned query via
   `helpers/openers.ts:openTableFromDbQuery(page,
   provisioned.queryNqName)`.
3. Go to **Data** > **Link Tables**:
   - Table 1: `demog`, Table 2: query result
   - Link column: any shared key column
   - Link type: **Selection to filter**
   - Click **OK**
4. **Verify:** linking works.
5. **Save with Data Sync ON:**
   - **Data sync** toggle **ON** for both tables.
   - Save as `Test_Case6_Sync`.
6. Close and reopen.
7. **Verify on reopen:** tables loaded, linking works, no console
   errors.
8. **Save with Data Sync OFF:**
   - Repeat steps 1–4, **Data sync OFF** for both tables.
   - Save as `Test_Case6_NoSync`.
9. Close and reopen.
10. **Verify on reopen:** tables loaded, linking works, no console
    errors.

---

> **Case 7 (Get Top 100 / Get All + DB table double-click + Join Tables)
> moved to `../Queries/get-all-get-top-100.md`.** This scenario still
> touches two known bugs: GROK-19103 (a derived table must land in the
> active project, not a stray one) via Cases 8 and 9, and GROK-18345
> (Spaces + Data Sync) via Cases 4-6. GROK-19212 (rename losing its
> data-sync link) is no longer covered here — see `complex.md`
> instead.

---

### Case 8: Browse > Files + Pivot Table (Add to workspace)

1. Open **Browse** > **Files** > `System:AppData` > `Chem` > `tests`.
2. Double-click `spgi-100.csv` — wait for the table to open.
3. Go to **Data** > **Pivot Table**.
4. Configure the pivot table:
   - Rows: select a categorical column (e.g. `Country`).
   - Columns: select another categorical column (e.g. `Industry`).
   - Values: select a numeric column with an aggregation
     (e.g. `Count` of `ID`).
5. Click **Add** (Add to workspace).
6. **Verify:** the pivot table opens as a separate table in a new tab.
7. **Save with Data Sync ON:**
   - **File** > **Save Project**.
   - **Data sync** toggle **ON** for the source table.
   - Save as `Test_Case8_Sync`.
8. Close and reopen.
9. **Verify on reopen:** both the source table and the pivot table are
   loaded, no console errors.
10. **Save with Data Sync OFF:**
    - Repeat steps 1–6, **Data sync OFF**.
    - Save as `Test_Case8_NoSync`.
11. Close and reopen.
12. **Verify on reopen:** both tables loaded, no console errors.

---

### Case 9: DB table double-click + Aggregate Rows (Add to workspace)

1. Open `public.groups` on the built-in `System:Datagrok` connection
   via `helpers/openers.ts:openTableFromDbTable(page, {connectionNqName:
   SYSTEM_DATAGROK_NQNAME, schemaName: 'public', tableName: 'groups'})`.
   Mirrors `Browse > Databases > ... > double-click` semantics.
2. Go to **Data** > **Aggregate Rows**.
3. Configure the aggregation:
   - Group by: `personal`.
   - Aggregation: `count` of `id`.
4. Click **Add** (Add to workspace) — or use
   `helpers/openers.ts:addAggregateToWorkspace({via: 'menu'})`.
5. **Verify:** the aggregated table opens as a separate table in a
   new tab.
6. **Save with Data Sync ON:**
   - **File** > **Save Project**.
   - **Data sync** toggle **ON** for the source table.
   - Save as `Test_Case9_Sync`.
7. Close and reopen.
8. **Verify on reopen:** both the source table and the aggregated
   table are loaded, no console errors.
9. **Save with Data Sync OFF:**
    - Repeat steps 1–5, **Data sync OFF**.
    - Save as `Test_Case9_NoSync`.
10. Close and reopen.
11. **Verify on reopen:** both tables loaded, no console errors.

## Notes

- **All 16 paths are the source of truth.** This file spells out all
  8 source-combination cases (each with a Sync-ON and Sync-OFF
  variant, minus Case 7 which moved to
  `../Queries/get-all-get-top-100.md`) — 16 saved projects total. The
  existing `uploading-spec.ts` currently automates only a subset
  (Cases 1, 3, 4, 9, Sync ON only); any such reduction is documented
  in that spec, not here.
- **Project naming.** All 16 saved projects use deterministic names
  (`Test_Case<N>_Sync` / `Test_Case<N>_NoSync`), which would collide
  across concurrent CI runs. The existing spec avoids this with a
  `Date.now()` timestamp suffix — automation should keep doing that.
  Cleanup of these projects is delegated to `projects-ui-smoke.md`.
- **Cases 4-6 create their own Space.** The `test-projects-demo` Space
  is created via JS API at the start of the case bundle and deleted
  at the end, so no Space needs to pre-exist on the server.
- **Cases 2/3/6 provision their own saved query.** Each case creates
  and deletes its own query on the built-in `System:Datagrok`
  connection rather than depending on the Samples package. For the
  two-query-source cases, this means running the same query twice
  rather than two genuinely different queries — an acceptable
  trade-off for not depending on external fixtures.
- **Cases 1/5/8 use `spgi-100.csv`** from `System:AppData/Chem/tests`
  in place of the original `SPGI_v2_infinity.csv` from
  `System:Demo`, for the same self-containment reason. Two-file
  cases become same-file-twice patterns; Datagrok auto-disambiguates
  the tab names (`spgi-100`, `spgi-100 (2)`).
- **Case 9's DB table** (`System:Datagrok / public.groups`) is a
  built-in connection and table — no provisioning needed.
- **The Data Sync toggle only applies to source tables.** Derived
  tables (the Pivot and Aggregate results in Cases 8-9) are always
  persisted as part of the project; they don't have a sync state of
  their own.
- **UI coverage split with `projects-ui-smoke.md`.** The Save Project
  dialog and Data Sync toggle are owned by `projects-ui-smoke.md`; this
  scenario owns the source-specific dialog flows (Link Tables, Join
  Tables, Pivot Table > Add, Aggregate Rows > Add). Get Top 100 /
  Get All flows are owned by `../Queries/get-all-get-top-100.md`.
