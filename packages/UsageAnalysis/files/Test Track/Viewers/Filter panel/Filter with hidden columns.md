### Verify common Hidden Columns Behavior in the Filter Panel
Ensure that columns are correctly shown or hidden in the filter panel based on their visibility settings.
* Preconditions: SPGI dataset is opened by pressing 'star' icon in TestTrack (with tooltip "Open test data").
Test Steps:
1. Check Initial Visible Columns: open the filter panel for the dataset and identify the currently visible columns in the filter panel (for e.g., 'Structure' colunm). 
2. Hide one of the columns that is currently visible in the filter panel. (for e.g., 'Structure' colunm).
3. Remove all filters and reopen panel, so that the auto script will execute
* Expected Result: The hidden column should no longer be visible in the filter panel.
Make Hidden Column Visible:
4. Open "Order or Hide Columns". Make the previously hidden column visible again. Reset the filters in the filter panel to refresh the list.
  * Expected Result: The column that was made visible should now appear in the filter panel.
5. Postconditions: The filter panel accurately reflects the visibility status of columns after they are hidden or made visible.
#### Verify special Hidden Columns Behavior in the Filter Panel
* Preconditions: SPGI dataset is loaded into the platform.
Test Steps:
1. Run Chem > R-group analysis
2. Open filter panel. 
Expected Result: Hidden special columns (name starts with "~") are not visible in the filter panel.

---
{
"order": 9,
"datasets": ["System:DemoFiles/SPGI.csv"]
}
