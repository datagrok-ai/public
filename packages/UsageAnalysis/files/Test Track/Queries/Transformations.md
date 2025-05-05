1. Go to **Browse** > **Databases** > **Postgres** > **NorthwindTest**
2. Right-click the **Products** query and select **Edit…** — the **Query View** opens
3. Open the **Transformation** tab
4. Click **Add new column** and add `${productid}`. Click **OK** — verify that the action is added to the transformation script
5. In **Toolbox** > **Actions**, select **Run query...**
6. Save the query
7. Close all views
8. Run the query — verify that the added transformation function executes
9. Delete the transformations created during this test
10. Save the query
11. Refresh the view — verify that the transformations are deleted

---
{
  "order": 6
}