1. Open **demog**
2. Open the **Filter Panel**
3. Add a new column: **SEX_bool**: `case ${SEX} when "F" then true else false end`  
4. Go to **Filter Panel > Hamburger menu > Add filter > Combined Boolean**
5. Verify that each boolean column in the dataset appears as a row in the Combined Boolean filter 
6. Verify that the filter displays count indicators for each value.  
7. Click the count number (not the checkbox) in the filter - **Corresponding rows are selected in the grid**
7. Apply the combined filter
9. Apply other filters (categorical, numerical)
10. Remove all filters
11. Close the **Filter Panel**
12. Open the **Filter Panel** - the Combined Boolean filter should be added automatically
12. Save/apply the layout

---
{
"order": 8,
"datasets": ["System:DemoFiles/demog.csv"]
}