1. Open SPGI
1. Open the Filter Panel

#### 1. Filtering and resetting

1. Apply any filters and verify that the data updates accordingly
2. Click the **Reset** icon and confirm that all filters are cleared
3. Apply new filters, disable one of them, and verify that the data updates accordingly
4. Close the **Filter Panel** — **filtering should be reset**
5. Reopen the **Filter Panel** and verify that:
   - The previously disabled filter remains disabled
   - Other filters are applied correctly
6. Remove some filters
7. Close and reopen the **Filter Panel** — **verify that the removed filters are no longer present**

#### 2. Filter Panel UI Behavior

1. Resize the **Filter Panel**
2. Add some viewers
1. Hover over the **categorical** filters and check that related data is highlighted on plots
3. Hover over a **date** filter - **verify that the tooltip shows the correct column format**
4. In the grid, change the date column's format - **verify that it is updated in the filter's tooltip**
4. For a categorical filter, click the **Search categories** icon — the search field should appear
5. Save/apply the layout - **verify that the search field visibility is preserved**

#### 3. Adding and Reordering Filters

1. Reorder filters by dragging them within the **Filter Panel**
2. Alternatively, reorder filters via **Hamburger menu > Reorder Filters**
3. Add or remove filters using **Hamburger menu > Add Columns**
4. Remove all filters using **Hamburger menu > Remove All**
5. Add a new filter to the **Filter Panel** — verify that new filters are added to the top:
    - Drag column headers from the grid
    - Grid: **Column header > Filter > Add Filter**
    - Grid: **Current Value > Use as filter**
    - Filter Panel: **Context action > Add filter > Scaffold Tree Filter**
    - Filter Panel: **Context action > Select Columns**
6. Use **Hamburger menu > Remove All** in the Filter Panel
7. Close the **Filter Panel**

#### 4. Hidden columns

1. Use the **Order or Hide Columns** dialog to hide some columns (e.g., all except 3–5)
2. Open the **Filter Panel** — verify that the hidden columns are no longer shown
3. Use **Order or Hide Columns** to make the previously hidden columns visible again
4. Remove all filters and reopen the **Filter Panel** — once made visible, the columns should reappear correctly

---
{
"order": 1,
"datasets": ["System:DemoFiles/SPGI.csv"]
}
