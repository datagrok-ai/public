1. Go to **Browse** >**Databases** > **Postgres** > **NorthwindTest**
3. Right-click the **Products** query and select **Edit…** from the context menu. A **Query View** opens
4. Click **Transformation** tab in Query View
5. Click **Add new column** from the actions list and add a column which is equal to the **productid** column from the table (${productid}). Press **OK**. Make sure that the action is added to the transformation script
6. On **Actions** tab on **Toolbox**, click **Run query...**
7. On the top toolbar, click the **Save** button to save query
8. Close all. Run the query saved in step 7. Make sure that the function added to the transformation script is executed for the returned dataframe
9. Delete transformations added during this test scenario, save changes. Refresh view, check that deletion was successful
---
{
  "order": 6
}