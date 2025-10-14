1. Go to **Databases > Postgres > NorthwindTest**.
2. Right-click the **NorthwindTest** and select **New Query...**:
- name: `Test_Postprocessing`
- body: `select * from products`
3. Run the query – verify that the result appears
4. Switch to the **Post-Process** tab
1. On line 7, add: `grok.shell.info(result.rowCount);`
5. Switch to the **Layout** tab
1. Add a scatterplot and correlation plot
6. Save the query
7. Close all
8. Go to  **Databases > Postgres > NorthwindTest** 
9. Preview and run the **Test_Postprocessing** query - verify that:
     * Layout preview shows both viewers
     * The green balloon with '77' appears
10. Right-click the **Test_Postprocessing** query and select *Edit..**
1. Run the query 
   * From the **Post-Process** tab — verify that a green balloon with **77** appears
   * From the **Layout** tab — verify that a green balloon with **77** appears and both viewers are displayed

---
{
  "order": 13
}