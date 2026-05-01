---
feature: projects
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/uploading.md
migration_date: 2026-04-30
migration_report: uploading-migration-report.md
related_bugs: [GROK-19103, GROK-18345]
---

# Uploading

Save 9 distinct two-table-source combinations as Datagrok projects, each
with both **Data Sync ON** and **Data Sync OFF** variants — **18 saved
projects total**. Verifies that project upload, save, reopen, and
inter-table linking/joining work across all major data sources: local
files, database query results, Spaces, file shares, get-top-N / get-all
on DB tables, pivot-table derivations, and aggregate-row derivations.

## Setup

A clean Datagrok session is the only shared setup. Each case opens its
own pair of tables, performs its own link/join, then saves and reopens.
Cleanup between cases is the responsibility of the test harness
(close all views before the next case).

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
| 7    | Get Top 100 / Get All             | DB table (double-click)               |
| 8    | Browse > Files                    | Pivot Table > Add to workspace        |
| 9    | DB table (double-click)           | Aggregate Rows > Add to workspace     |

Each case below is run twice — once with **Data Sync ON** (saved as
`Test_Case<N>_Sync`) and once with **Data Sync OFF** (saved as
`Test_Case<N>_NoSync`). Both variants are followed by a close-and-reopen
verification.

---

### Case 1: Browse > Files + Browse > Files (Link Tables)

1. Open **Browse** > **Files** > `System:Demo`.
2. Double-click `SPGI_v2_infinity.csv` — wait for the table to open.
3. Go back to **Browse** > **Files** > `System:Demo`.
4. Double-click `SPGI_v2.csv` — wait for the table to open.
5. Go to **Data** > **Link Tables**:
   - Table 1: `SPGI_v2_infinity`, Table 2: `SPGI_v2`
   - Link column: `ID`
   - Link type: **Selection to filter**
   - Click **OK**
6. **Verify:** select rows in `SPGI_v2_infinity` — corresponding rows
   in `SPGI_v2` are filtered.
7. **Save with Data Sync ON:**
   - **File** > **Save Project** (or Ctrl+S)
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

1. Open **Browse** > **Databases** > **Postgres** > **NorthwindTest**.
2. Find the query `PostgressAll`.
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

1. Open **Browse** > **Databases** > **Postgres** > **NorthwindTest**.
2. Run query `PostgressAll` — wait for the result table to open.
3. Open **Browse** > **Files** > `System:Demo`.
4. Double-click `SPGI_v2.csv` — wait for the table to open.
5. Go to **Data** > **Link Tables**:
   - Table 1: query result, Table 2: `SPGI_v2`
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

1. Open the **SPGIs** space: navigate to **Browse** > **Spaces** >
   **SPGIs** (or open `https://dev.datagrok.ai/s/SPGIs`).
2. Double-click `SPGI.csv` — wait for the table to open.
3. Go back to the **SPGIs** space.
4. Double-click `SPGI_v2.csv` — wait for the table to open.
5. Go to **Data** > **Link Tables**:
   - Table 1: `SPGI`, Table 2: `SPGI_v2`
   - Link column: `ID`
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

1. Open the **SPGIs** space and double-click `SPGI.csv`.
2. Open **Browse** > **Files** > `System:Demo` and double-click
   `SPGI_v2.csv`.
3. Go to **Data** > **Link Tables**:
   - Table 1: `SPGI`, Table 2: `SPGI_v2`
   - Link column: `ID`
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

1. Open the **SPGIs** space and double-click `SPGI.csv`.
2. Open **Browse** > **Databases** > **Postgres** > **NorthwindTest**
   and run query `PostgressAll`.
3. Go to **Data** > **Link Tables**:
   - Table 1: `SPGI`, Table 2: query result
   - Link column: `ID`
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

### Case 7: Get Top 100 / Get All + DB table double-click (Join Tables)

1. Open **Browse** > **Databases** > **Postgres** > **NorthwindTest** >
   **Schema** > **public**.
2. Right-click on the `orders` table and select **Get Top 100** —
   wait for the table to open.
3. Right-click on the same `orders` table and select **Get All** —
   wait for the second table to open.
4. Add a calculated column to the first table:
   - Select the first table, go to **Data** > **Add New Column**.
   - Name: `key1`, Formula: `${orderid} - 5`, click **OK**.
5. Add a calculated column to the second table:
   - Select the second table, go to **Data** > **Add New Column**.
   - Name: `key2`, Formula: `${orderid} + 5`, click **OK**.
6. Go to **Data** > **Join Tables**:
   - Table 1: first table (`orders` Get Top 100), Table 2: second
     table (`orders` Get All).
   - Key column 1: `key1`, Key column 2: `key2`.
   - Join type: **Inner**.
   - Click **OK**.
7. Repeat step 6 for join types **Outer**, **Left**, **Right**
   (4 joins total).
8. **Verify:** all four joined result tables are visible in the
   workspace.
9. **Save with Data Sync ON:**
   - **File** > **Save Project**.
   - **Data sync** toggle **ON** for applicable tables.
   - Save as `Test_Case7_Sync`.
10. Close and reopen.
11. **Verify on reopen:** all tables (source + 4 joined) are loaded,
    no console errors.
12. **Save with Data Sync OFF:**
    - Repeat steps 1–8, **Data sync OFF**.
    - Save as `Test_Case7_NoSync`.
13. Close and reopen.
14. **Verify on reopen:** all tables loaded, no console errors.

---

### Case 8: Browse > Files + Pivot Table (Add to workspace)

1. Open **Browse** > **Files** > `System:Demo`.
2. Double-click `SPGI_v2.csv` — wait for the table to open.
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

- Original `order: 1` — produces 18 saved projects
  (`Test_Case<N>_Sync` and `Test_Case<N>_NoSync` for N in 1..9).
  These names are not consumed by other migrated scenarios in the
  Projects section as of scenario-chains rev 2; share-project.md and
  later cases consume the `demog` project produced by
  `upload-project.md` instead. Surfaced as candidate fixture set
  `uploading-18-projects` for chain analysis if any future scenario
  needs a wide projects-list fixture.
- Cases 4–6 require a pre-configured **SPGIs** Space on the test
  server; on environments without it, Cases 4–6 will fail at the
  Spaces-navigation step. This is an environment dependency, not a
  test-content dependency. The original 2026-03-09 run flagged this
  by skipping Cases 4–8 ("No Spaces set up on release server").
- Cases 7–9 produce derivative tables (joins, pivots, aggregates)
  alongside the source tables; the `Data sync` toggle applies to the
  source tables only — derivative tables are persisted as part of
  the project state, not synced from disk.
- Helpers: the `uploadProject(projectName, tableInfo, view, df)`
  helper at `public/packages/UITests/src/gui/gui-utils.ts:100` is
  registered as `grok_test_layer` (a/b style) and is NOT
  Playwright-compatible. The existing `uploading-spec.ts` covers
  Cases 1, 3, 4, 9 (Sync ON only) at the playwright layer using a
  local `saveProject(page, name)` helper. A Playwright-layer
  counterpart for the registry is flagged as a candidate in the
  migration report.
