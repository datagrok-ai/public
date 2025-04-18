#### Query Postprocessing 

1. Go to **Databases > Postgres > NorthwindTest**.
2. RMB on NorthwindTest > **New Query...**:
- name: Test_Postprocessing
- body: select * from products
3. Press **Run query** button (F5 or Ctrl+Enter) – result grid should appear.
4. Switch to **Post-Process** tab. On line 7, set: grok.shell.info(result.rowCount);
5. Switch to **Layout** tab. Add Scatter plot and Correlation plot from Toolbox.
6. **Save** the query.
7. **Top menu > Close all**. 
8. In Databases > Postgres > NorthwindTest, locate and click Test_Postprocessing query.
  - **Expected result:** Layout preview shows both viewers. The green balloon with '77' should appear. 
9. Run the Test_Postprocessing query from Browse. 
  - **Expected result:** The green balloon with '77' should appear. The preview should contain two added viewers.
10. Edit the Test_Postprocessing query from Browse. Press **Run query** button. 
  - **Expected result:** Switching to **Post-Process** and **Layout** tabs triggers green balloon with **77**. **Layout** tab opens with both viewers.

---
{
  "order": 13
}