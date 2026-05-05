---
feature: projects
sub_features_covered:
  - projects.upload
  - projects.api.save
  - projects.api.files.sync
  - projects.add-relation
target_layer: playwright
coverage_type: regression
pyramid_layer: source-matrix
ui_coverage_responsibility:
  - save-project-dialog
  - link-tables-dialog
  - join-tables-dialog
  - pivot-table-add-to-workspace
  - aggregate-rows-add-to-workspace
ui_coverage_delegated_to: upload-project.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/uploading.md
migration_date: 2026-05-04
migration_report: uploading-migration-report.md
related_bugs:
  - GROK-19103
  - GROK-18345
---

# Uploading

Save 8 distinct two-table-source combinations as Datagrok projects,
each with both **Data Sync ON** and **Data Sync OFF** variants — **16
saved projects total**. Verifies that project upload, save, reopen,
and inter-table linking/joining work across all major data sources:
local files, database query results, Spaces, file shares,
pivot-table derivations, and aggregate-row derivations.

DB-table right-click flows (`Get Top 100`, `Get All`) are covered
by their dedicated scenario `../Queries/get-all-get-top-100.md`;
this scenario uses plain DB-table double-click for DB sources.

`pyramid_layer: source-matrix` per
`scenario-chains/projects.yaml` rev 3 — this scenario enumerates
table-source combinations rather than UI states. Save Project dialog
smoke is owned by `upload-project.md` (chain
`ui_coverage_plan.smoke_scenario`); this scenario adds the
source-specific dialog flows (Link Tables, Join Tables,
Pivot/Aggregate Add to workspace) that do not appear in
upload-project.md's smoke surface.

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
Sync toggle smoke surface is owned by `upload-project.md`; this
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

1. Open **Browse** > **Platform** > **Functions** > **Queries**.
2. Find the query `Samples:PostgresCustomers`.
3. Run it (right-click > **Run** or double-click) — wait for the result
   table to open.
4. Run the same query again to get a second result table in the
   workspace.
5. Go to **Data** > **Link Tables**:
   - Table 1: first query result, Table 2: second query result
   - Link column: `ID`
   - Link type: **Selection to filter**
   - Click **OK**
6. **Verify:** select rows in Table 1 — corresponding rows in Table 2
   are filtered.
7. **Save with Data Sync ON:**
   - **File** > **Save Project**.
   - Ensure the **Data sync** toggle is **ON** for both tables (the
     CREATION SCRIPT block should be visible).
   - Save as `Test_Case2_Sync`, click **OK**.
8. Close the project. Reopen from **Browse** > **Projects**.
9. **Verify on reopen:** both tables are loaded, linking works,
   no console errors.
10. **Save with Data Sync OFF:**
    - Repeat steps 1–6.
    - Turn the **Data sync** toggle **OFF** for both tables.
    - Save as `Test_Case2_NoSync`.
11. Close and reopen.
12. **Verify on reopen:** tables loaded, linking works, no console
    errors.

---

### Case 3: Query result + Browse > Files (Link Tables)

1. Open **Browse** > **Platform** > **Functions** > **Queries**.
2. Run query `Samples:PostgresCustomers` — wait for the result table
   to open.
3. Open **Browse** > **Files** > `System:AppData` > `Chem` > `tests`.
4. Double-click `spgi-100.csv` — wait for the table to open.
5. Go to **Data** > **Link Tables**:
   - Table 1: query result, Table 2: `spgi-100`
   - Link column: `ID`
   - Link type: **Selection to filter**
   - Click **OK**
6. **Verify:** linking works.
7. **Save with Data Sync ON:**
   - **File** > **Save Project**.
   - **Data sync** toggle **ON** for both tables.
   - Save as `Test_Case3_Sync`.
8. Close and reopen.
9. **Verify on reopen:** tables loaded, linking works, no console
   errors.
10. **Save with Data Sync OFF:**
    - Repeat steps 1–6, **Data sync** toggle **OFF** for both tables.
    - Save as `Test_Case3_NoSync`.
11. Close and reopen.
12. **Verify on reopen:** tables loaded, linking works, no console
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

1. Open the **test-projects-demo** space and double-click `demog.csv`.
2. Open **Browse** > **Platform** > **Functions** > **Queries** and
   run query `Samples:PostgresCustomers`.
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
> moved to `../Queries/get-all-get-top-100.md`.** Cross-cutting bug
> spans previously anchored here:
> - `GROK-19103` — derivation lands in active project — still anchored
>   in this scenario via Cases 8 (Pivot Add to workspace) and 9
>   (Aggregate Rows Add to workspace), and externally in
>   `upload-project.md:Step 1`, `projects-copy-clone.md:Step 4`,
>   `complex.md:Step 1`.
> - `GROK-19212` — rename → reference resolution under datasync —
>   no longer anchored in this scenario; sole anchor is now
>   `complex.md:Step 7`.
> - `GROK-18345` — Spaces+sync — anchored in this scenario via
>   Cases 4–6 (Spaces sources).

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

1. Open **Browse** > **Databases** > **Postgres** > **NorthwindTest** >
   **Schema** > **public**.
2. Double-click the `orders` table — wait for the table to open.
3. Go to **Data** > **Aggregate Rows**.
4. Configure the aggregation:
   - Group by: `customerid`.
   - Aggregation: `count` of `orderid`.
5. Click **Add** (Add to workspace).
6. **Verify:** the aggregated table opens as a separate table in a
   new tab.
7. **Save with Data Sync ON:**
   - **File** > **Save Project**.
   - **Data sync** toggle **ON** for the source table.
   - Save as `Test_Case9_Sync`.
8. Close and reopen.
9. **Verify on reopen:** both the source table and the aggregated
   table are loaded, no console errors.
10. **Save with Data Sync OFF:**
    - Repeat steps 1–6, **Data sync OFF**.
    - Save as `Test_Case9_NoSync`.
11. Close and reopen.
12. **Verify on reopen:** both tables loaded, no console errors.

## Notes

- **D-STRUCT-02 — all 16 paths preserved at the .md level.** Per
  chain `ui_coverage_plan` rev 3 and prior decision-log
  `mig-2026-04-30-uploading-migration` (4-of-9 reduction was
  Automator-stage, not Migrator-stage). The migrated .md remains the
  full 16-path source of truth (Case 7 split out to
  `../Queries/get-all-get-top-100.md`); any spec-side reduction
  (e.g. the existing `uploading-spec.ts` covers Cases 1, 3, 4, 9 with
  Sync ON only) is documented in the spec header, not in this .md.
- **Project naming:** all 16 saved projects use deterministic names
  (`Test_Case<N>_Sync` / `Test_Case<N>_NoSync`). Concurrent CI runs
  will collide. The existing `uploading-spec.ts` uses a `Date.now()`
  suffix (`AutoTest-Upload-Case<N>-<timestamp>`) to avoid this; the
  Automator should preserve the timestamp-suffix convention at spec
  time. Cleanup is delegated to `deleting.md` (must_run_last per
  chain rev 3).
- **Cases 4–6 Spaces prelude (decision-log
  `sa-2026-05-03-spaces-inline-prelude-pattern`):** the Space is
  created via JS API at the start of the case bundle and deleted at
  the end, so no pre-existing server-side Space is required. The
  original 2026-03-09 run flagged the prior pre-provisioned-Space
  formulation as blocked ("No Spaces set up on release server"); the
  inline-prelude pattern closes that env dependency.
- **Cases 2/3/6 query substitution (decision-log
  `sa-2026-05-03-postgres-queries-public-data-substitution`):** the
  former `Samples:PostgresAll` reference is substituted with
  `Samples:PostgresCustomers`. Two-table-source cases become
  same-query-twice patterns after substitution; structural test value
  is preserved, semantic distinctness of two query results is
  degraded — flagged for downstream review.
- **Cases 1/5/8 file substitution (decision-log
  `sa-2026-05-03-spgi-file-public-data-substitution`):** the former
  `SPGI_v2_infinity.csv` from `System:Demo` is substituted with
  `spgi-100.csv` from `System:AppData/Chem/tests`. Two-file-source
  cases become same-file-twice patterns after substitution; Datagrok
  auto-disambiguates table-tab names as `spgi-100` and `spgi-100 (2)`.
- **Case 9 DB-table reference** still points at
  `Postgres > Northwind > public > orders` directly via plain
  double-click (the source-class mapping does not cover bare DB-table
  references). May be mapped to `Samples:PostgresOrders` in a future
  pass or accepted as env-dependent skip on dev.
- **Cases 8–9 derivative tables** (pivots, aggregates) are persisted
  as part of the project state alongside the source tables; the
  **Data sync** toggle applies to the source tables only —
  derivative tables are persisted, not synced from disk.
- **UI coverage delegation per chain rev 3:** Save Project dialog +
  Data Sync toggle smoke is owned by `upload-project.md`. This
  scenario owns the source-specific dialog flows (Link Tables, Join
  Tables, Pivot Table > Add, Aggregate Rows > Add) — these are NOT
  covered by `upload-project.md`'s smoke surface. Get Top 100 / Get
  All flows are owned by `../Queries/get-all-get-top-100.md`.
- **Helpers:** the legacy `uploadProject(projectName, tableInfo,
  view, df)` helper at
  `public/packages/UITests/src/gui/gui-utils.ts:100` is registered as
  `grok_test_layer` and is NOT Playwright-compatible. The existing
  `uploading-spec.ts` uses a local `saveProject(page, name)` helper.
  A Playwright-layer counterpart
  (`helpers.playwright.projects.saveAndReopen(page, name, syncOn)`)
  is flagged as a candidate in the migration report (not invented by
  this Migrator).
