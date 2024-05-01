1. Go to **Browse** >**Databases**
2. Expand the **Postgres** > **northwindTest**
3. Right-click the **Products** query and select **Edit…** from the context menu. A **Query View** opens.
4. Click **Transformation** tab in Query View.
5. Click **Add new column** from the actions list and add a column which is equal to the **productid** column from the table (${productid})
6. On the menu ribbon select **ML** | **Custer…**. A **Cluster** dialog opens.
7. In Cluster dialog, click OK.
8. Make sure that the actions from steps 6 and 7 are added to the transformation script.
9. On the top toolbar, click the **Play** button.
10. On **Actions** tab on **Toolbox**, click **Run** query…
11. On the top toolbar, click the **Save** button.
12. Run the query saved in step 11
13. Make sure that the functions added to the transformation script are executed for the returned dataframe
---
{
  "order": 5
}