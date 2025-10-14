#### Matrix: Table source combinations

| Case | Table 1 source                 | Table 2 source              |
|------|---------------------------------|-----------------------------|
| 1    | Local storage                  | Local storage               |
| 2    | Local storage                  | Browse > Files              |
| 3    | Browse > Files                 | Browse > Files              |
| 4    | Query result (Run / Double-click) | Query result (Run / Double-click) |
| 5    | Query result (Run / Double-click) | Browse > Files              |
| 6     | Spaces                        | Spaces                     |
| 7     | Spaces                        | Browse > Files             |
| 8     | Spaces                        | Query result (Run / Double-click)              |
| 9    | Get All / Get Top 100          | DB table (double-click)     |
| 10    | Browse > Files| Pivot Table > Add to workspace                            |
| 11    | DB table (double-click) |  Aggregate Rows > Add to workspace                       |

#### 1. Generic steps for cases 1–8 (Link tables)

1. Open **Table 1 and 2** from the defined sources
2. Go to **Data > Link Tables** 
   - Define link columns (e.g. `Sample Name`, `link column 1`, `link column 2`) 
   - Set link type (e.g. **Selection to filter**)
  
3. Save the project:  
   - With data sync (if possible)  
   - Without data sync
4. Reopen and verify the saved projects

### Steps for case 9 (Join tables)

1. Open **Table 1 and 2** from the defined sources
2. Add a calculated column for each table to use it as a Key column to join tables (e.g. `${orderid}-5` and `${orderid}+5`)
2. Use **Data > Join Tables** :  
   - Test all join types: **inner**, **outer**, **left**, **right**  
3. Verify all four resulting joined table are displayed in the workspace
4. Save the project (with and without data sync)
5. Reopen and verify the saved projects

### Steps for case 10: Pivot Table > Add to workspace

1. Open table `SPGI.csv` from **Browse > Files**
2. Add a pivot table and configure it
3. Click **Add**  - Verify pivot table is opened in a separate tab 
4. Save the project (with and without data sync)
5. Reopen and verify the saved projects

### Steps for case 11: Aggregate Rows > Add to workspace

1. Go to Databases > Postgres > Northwind Test > Schema > public and double-click the orders table  
2. Use **Data> Aggregate Rows**
3. Click **Add** - Verify aggregated table is opened in a separate tab.  
4. Save the project (with and without data sync)
5. Reopen and verify the saved projects

---
{
  "order": 1
}