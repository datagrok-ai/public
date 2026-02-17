1. Go to **Browse** > **Databases** > **Postgres** > **NorthwindTest** > **Schemas** > **public**. A list of tables opens
3. Right-click the **customers** table and select **New Visual Query...** from the context menu
5. Set **Group by** to `companyname`
5. Start setting **Order by**  - check that only `companyname` column can be selected
6. Go to **Data** and select `orders.shipcity` and `orders.shipcountry`
7. Set **Where** to `customers.companyname contains C` and select the **Expose as function parameter** checkbox
7. Set **Aggregate** to `sum(orders.freight)`
7. Add **Group by**: `orders.shipcounty`
7. Set **Pivot** to `orders.shipcity`
7. Add **Order by**: `orders.shipcounty`
7. Set different order for the **Order by** fields - check
8. Open the **Debug** tab and press a bug icon to **debug query**. No errors should appear
9. **Toolbar > Actions** > **Run query...** - query result should open in new view
8. Save the query. Close All
8. Run the saved query - check the parameter
8. On the **Toolbar**, change the parameter value and click **Refresh** - check
8. Close all

**Cross-catalog visual query (MS SQL)**
1. Go to **Browse** > **Databases** > **MS SQL** > **NorthwindTest** > **Catalogs** > pick a catalog > pick a schema. A list of tables opens
2. Right-click a table and select **New Visual Query...** from the context menu
3. In **Data**, add a column from a table in a different catalog
4. Verify the cross-catalog join is reflected in the **Debug** tab
5. Run the query — it should return results without errors
---
{
  "order": 11
}
