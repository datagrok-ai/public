#### Matrix: Table source combinations

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

---

### Case 1: Browse > Files + Browse > Files (Link Tables)

1. Open **Browse** > **Files** > `System:Demo`
2. Double-click `SPGI_v2_infinity.csv` — wait for the table to open
3. Go back to **Browse** > **Files** > `System:Demo`
4. Double-click `SPGI_v2.csv` — wait for the table to open
5. Go to **Data** > **Link Tables**
   - Table 1: `SPGI_v2_infinity`, Table 2: `SPGI_v2`
   - Link column: `ID`
   - Link type: **Selection to filter**
   - Click **OK**
6. Verify: select rows in `SPGI_v2_infinity` — corresponding rows in `SPGI_v2` should be filtered
7. **Save with Data Sync:**
   - **File** > **Save Project** (or Ctrl+S)
   - In the Save Project dialog, make sure **Data sync** toggle is **ON** for both tables
   - Give the project a name (e.g. `Test_Case1_Sync`) and click **OK**
8. Close the project. Reopen it from **Browse** > **Projects**
9. **Verify:** both tables are loaded, linking works (select rows in one table — the other filters), no errors in the browser console (F12)
10. **Save without Data Sync:**
    - Repeat steps 1–6
    - In the Save Project dialog, turn **Data sync** toggle **OFF** for both tables
    - Save as `Test_Case1_NoSync`
11. Close the project. Reopen it from **Browse** > **Projects**
12. **Verify:** both tables are loaded, linking works, no console errors

---

### Case 2: Query result + Query result (Link Tables)

1. Open **Browse** > **Databases** > **Postgres** > **NorthwindTest**
2. Find the query `PostgressAll`
3. Run it (right-click > **Run** or double-click) — wait for the result table to open
4. Run the same query again to get a second result table in the workspace
5. Go to **Data** > **Link Tables**
   - Table 1: first query result, Table 2: second query result
   - Link column: `ID`
   - Link type: **Selection to filter**
   - Click **OK**
6. Verify: select rows in Table 1 — corresponding rows in Table 2 should be filtered
7. **Save with Data Sync:**
   - **File** > **Save Project**
   - Make sure **Data sync** toggle is **ON** for both tables (CREATION SCRIPT should be visible)
   - Save as `Test_Case2_Sync`, click **OK**
8. Close the project. Reopen from **Browse** > **Projects**
9. **Verify:** both tables are loaded, linking works, no console errors
10. **Save without Data Sync:**
    - Repeat steps 1–6
    - Turn **Data sync** toggle **OFF** for both tables
    - Save as `Test_Case2_NoSync`
11. Close and reopen. **Verify:** tables loaded, linking works, no console errors

---

### Case 3: Query result + Browse > Files (Link Tables)

1. Open **Browse** > **Databases** > **Postgres** > **NorthwindTest**
2. Run query `PostgressAll` — wait for the result table to open
3. Open **Browse** > **Files** > `System:Demo`
4. Double-click `SPGI_v2.csv` — wait for the table to open
5. Go to **Data** > **Link Tables**
   - Table 1: query result, Table 2: `SPGI_v2`
   - Link column: `ID`
   - Link type: **Selection to filter**
   - Click **OK**
6. Verify linking works
7. **Save with Data Sync:**
   - **File** > **Save Project**
   - **Data sync** toggle **ON** for both tables
   - Save as `Test_Case3_Sync`
8. Close and reopen. **Verify:** tables loaded, linking works, no console errors
9. **Save without Data Sync:**
   - Repeat steps 1–6, **Data sync** toggle **OFF** for both tables
   - Save as `Test_Case3_NoSync`
10. Close and reopen. **Verify:** tables loaded, linking works, no console errors

---

### Case 4: Spaces + Spaces (Link Tables)

1. Open the **SPGIs** space: navigate to **Browse** > **Spaces** > **SPGIs** (or go to `https://dev.datagrok.ai/s/SPGIs`)
2. Double-click `SPGI.csv` — wait for the table to open
3. Go back to the **SPGIs** space
4. Double-click `SPGI_v2.csv` — wait for the table to open
5. Go to **Data** > **Link Tables**
   - Table 1: `SPGI`, Table 2: `SPGI_v2`
   - Link column: `ID`
   - Link type: **Selection to filter**
   - Click **OK**
6. Verify linking works
7. **Save with Data Sync:**
   - **File** > **Save Project**
   - **Data sync** toggle **ON** for both tables
   - Save as `Test_Case4_Sync`
8. Close and reopen. **Verify:** tables loaded, linking works, no console errors
9. **Save without Data Sync:**
   - Repeat steps 1–6, **Data sync** toggle **OFF** for both tables
   - Save as `Test_Case4_NoSync`
10. Close and reopen. **Verify:** tables loaded, linking works, no console errors

---

### Case 5: Spaces + Browse > Files (Link Tables)

1. Open the **SPGIs** space and double-click `SPGI.csv`
2. Open **Browse** > **Files** > `System:Demo` and double-click `SPGI_v2.csv`
3. Go to **Data** > **Link Tables**
   - Table 1: `SPGI`, Table 2: `SPGI_v2`
   - Link column: `ID`
   - Link type: **Selection to filter**
   - Click **OK**
4. Verify linking works
5. **Save with Data Sync:**
   - **Data sync** toggle **ON** for both tables
   - Save as `Test_Case5_Sync`
6. Close and reopen. **Verify:** tables loaded, linking works, no console errors
7. **Save without Data Sync:**
   - Repeat steps 1–4, **Data sync OFF**, save as `Test_Case5_NoSync`
8. Close and reopen. **Verify:** tables loaded, linking works, no console errors

---

### Case 6: Spaces + Query result (Link Tables)

1. Open the **SPGIs** space and double-click `SPGI.csv`
2. Open **Browse** > **Databases** > **Postgres** > **NorthwindTest** and run query `PostgressAll`
3. Go to **Data** > **Link Tables**
   - Table 1: `SPGI`, Table 2: query result
   - Link column: `ID`
   - Link type: **Selection to filter**
   - Click **OK**
4. Verify linking works
5. **Save with Data Sync:**
   - **Data sync** toggle **ON** for both tables
   - Save as `Test_Case6_Sync`
6. Close and reopen. **Verify:** tables loaded, linking works, no console errors
7. **Save without Data Sync:**
   - Repeat steps 1–4, **Data sync OFF**, save as `Test_Case6_NoSync`
8. Close and reopen. **Verify:** tables loaded, linking works, no console errors

---

### Case 7: Get Top 100 / Get All + DB table double-click (Join Tables)

1. Open **Browse** > **Databases** > **Postgres** > **NorthwindTest** > **Schema** > **public**
2. Right-click on the `orders` table and select **Get Top 100** — wait for the table to open
3. Right-click on the same `orders` table and select **Get All** — wait for the second table to open
4. Add a calculated column to the first table:
   - Select the first table, go to **Data** > **Add New Column**
   - Name: `key1`, Formula: `${orderid} - 5`, click **OK**
5. Add a calculated column to the second table:
   - Select the second table, go to **Data** > **Add New Column**
   - Name: `key2`, Formula: `${orderid} + 5`, click **OK**
6. Go to **Data** > **Join Tables**
   - Table 1: first table (`orders` Get Top 100), Table 2: second table (`orders` Get All)
   - Key column 1: `key1`, Key column 2: `key2`
   - Join type: **Inner**
   - Click **OK**
7. Repeat step 6 for join types: **Outer**, **Left**, **Right** (4 joins total)
8. **Verify:** all four joined result tables are visible in the workspace
9. **Save with Data Sync:**
   - **File** > **Save Project**
   - **Data sync** toggle **ON** for applicable tables
   - Save as `Test_Case7_Sync`
10. Close and reopen. **Verify:** all tables (source + 4 joined) are loaded, no console errors
11. **Save without Data Sync:**
    - Repeat steps 1–8, **Data sync OFF**, save as `Test_Case7_NoSync`
12. Close and reopen. **Verify:** all tables loaded, no console errors

---

### Case 8: Browse > Files + Pivot Table (Add to workspace)

1. Open **Browse** > **Files** > `System:Demo`
2. Double-click `SPGI_v2.csv` — wait for the table to open
3. Go to **Data** > **Pivot Table**
4. Configure the pivot table:
   - Rows: select a categorical column (e.g. `Country`)
   - Columns: select another categorical column (e.g. `Industry`)
   - Values: select a numeric column with an aggregation (e.g. `Count` of `ID`)
5. Click **Add** (Add to workspace)
6. **Verify:** the pivot table is opened as a separate table in a new tab
7. **Save with Data Sync:**
   - **File** > **Save Project**
   - **Data sync** toggle **ON** for the source table
   - Save as `Test_Case8_Sync`
8. Close and reopen. **Verify:** both the source table and pivot table are loaded, no console errors
9. **Save without Data Sync:**
   - Repeat steps 1–6, **Data sync OFF**, save as `Test_Case8_NoSync`
10. Close and reopen. **Verify:** both tables loaded, no console errors

---

### Case 9: DB table double-click + Aggregate Rows (Add to workspace)

1. Open **Browse** > **Databases** > **Postgres** > **NorthwindTest** > **Schema** > **public**
2. Double-click the `orders` table — wait for the table to open
3. Go to **Data** > **Aggregate Rows**
4. Configure the aggregation:
   - Group by: `customerid`
   - Aggregation: `count` of `orderid`
5. Click **Add** (Add to workspace)
6. **Verify:** the aggregated table is opened as a separate table in a new tab
7. **Save with Data Sync:**
   - **File** > **Save Project**
   - **Data sync** toggle **ON** for the source table
   - Save as `Test_Case9_Sync`
8. Close and reopen. **Verify:** both the source table and aggregated table are loaded, no console errors
9. **Save without Data Sync:**
   - Repeat steps 1–6, **Data sync OFF**, save as `Test_Case9_NoSync`
10. Close and reopen. **Verify:** both tables loaded, no console errors

---
{
  "order": 1
}
