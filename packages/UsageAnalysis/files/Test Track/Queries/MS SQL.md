### 1. Adding
1. Go to **Browse** > **Databases** > **MS SQL**
3. Right-click the `NorthwindTest` connection and select **New query** from the context menu
4. Enter `test_query_ms_sql` to the Name field
5. Enter query: `select * from products`
1. Run the query using both:
  * **Menu Ribbon > Play button** — result appears at the bottom of the current view
  * **Toolbox > Actions > Run query…** — result opens in a new view
8. Save the query

### 2. Editing
1. Refresh view on **Browse**. Right-click the query from the previous step and select **Edit…** from context menu
2. Change name to `new_test_query_ms_sql`
3. Change the query test to: `select * from orders`
1. Run the query using both:
   * **Menu Ribbon > Play button** — result appears at the bottom of the current view
   * **Toolbox > Actions > Run query…** — result opens in a new view
8. Save the query

### 3. Browse
1. Refresh **Browse**
2. Go to **Browse** > **Databases** > **MS SQL** > **NorthwindTest**
3. Type `new_test` in the search field to search for the query from the previous step
4. On the **Context Panel**, check all tabs for the query

### 4. Deleting
1. Go to **Browse** > **Databases** > **MS SQL** > **NorthwindTest** and find query from the previous steps with name `new_test_query_ms_sql`
2. Delete it:
    * Right-click the connection and select **Delete** from the context menu
    * In the confirmation dialog, click DELETE
3. Refresh Browse and check that query has been deleted and is no longer present

### 5. Catalog browsing
1. Go to **Browse** > **Databases** > **MS SQL** > **NorthwindTest** > **Catalogs**
2. Verify that catalogs are listed under the connection
3. Expand a catalog — schemas should appear
4. Expand a schema — tables should appear
5. Click on a catalog node and verify the **Context Panel** shows the catalog preview
6. Check that **Comment** and **LLM comment** meta properties can be assigned to the catalog
---
{
  "order": 5
}