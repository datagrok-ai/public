### Catalog browsing

1. Go to **Browse** > **Databases**
2. Expand a connection that supports catalogs (e.g. **MS SQL** > **NorthwindTest**)
3. Verify that the tree shows **Catalogs** node (with a double-database icon) instead of **Schemas**
4. Expand **Catalogs** — a list of available database catalogs should appear
5. Expand one of the catalogs — schemas should load under it
6. Expand a schema — tables should load as usual

### Catalog preview

7. Click on a catalog node in the tree
8. Check the **Context Panel** — it should show a preview of the catalog
9. Verify the catalog name is displayed correctly

### Catalog meta properties

10. Select a catalog in the tree
11. In the **Context Panel**, check that **Comment** and **LLM comment** meta properties can be assigned to the catalog
12. Add a **Comment** to the catalog, click away, then click the catalog again — the comment should persist
13. Verify **LLM comment** can also be set and is displayed correctly

### Catalog schema view

14. Right-click a catalog node and select **Browse** from the context menu
15. A schema view should open showing all schemas within the catalog
16. Right-click a catalog node and select **Open as table**
17. A table with all database objects in the catalog should open

---
{
  "order": 9
}
