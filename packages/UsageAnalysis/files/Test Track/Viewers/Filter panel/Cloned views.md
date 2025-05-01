
1. Open SPGI
1. Open the **Filter Panel**
1. **View > Layout > Clone View**.
3. On **View 1**, for any numeric filter set:
   - either 'Filter out missing values',
   - or 'Show only missing values'.
4. GO to **View 2** and check missing values settings synchronization.
8. On any view, apply filtering - **rows should be filtered on both views**
9. Go to another view
10. Turn off all filters - **all rows should become visible on both views**
10. Turn on filters again
11. In any view, turn off an individual filter - **rows filtered by this filter should reappear**
12. On **View 2**, sketch any structure
14. Go to **View 1** and check the **Structure** filter
16. Close the **Filter Panel**
17. Go to **View 2** and modify the **Structure** filter value
18. Go to **View 1** and reopen the **Filter Panel** - the structure should match the one in **View 2**
1. Go to **View 2** and remove the **Structure** filter
2. Close / reopen the **Filter Panel** - the removed filter should NOT appear
3. Save the layout 
4. Add some more filtering.
5. Apply the previously saved layout - **both layout and original filtering should be restored correctly**

---
{
"order": 3,
"datasets": ["System:DemoFiles/SPGI.csv"]
}
