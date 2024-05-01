1. Go to **Browse** > **Databases**.
2. Expand **PostgresDart** > **northwind** > **Schemas** > **public**. A list of tables opens.
3. Right-click the **products** DB table and select **Visual Q**uery from the context menu.
4. Fill **Name** field with `test_visual_query`.
5. Select:
    - the `productid` column for the **Columns** field
    - the`suplplierid` column for the **Rows** field
    - avg(unitprice) for the **Measures** field
6. On the top toolbar, click the **Play** button.
7. On Toolbox, click the **Run query**... action.
8. On the top toolbar, click the **Save** button.
9. Go to **Data** > **Databases** and run **test_visual_query** (double-click it).
---
{
  "order": 3
}