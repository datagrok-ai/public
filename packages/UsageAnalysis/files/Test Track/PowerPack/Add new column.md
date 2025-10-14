### Add new columns

1. Open the Demog Dataset:
- Press the 'star' icon in TestTrack (with tooltip "Open test data"). The Demog dataset opens.
2. Press the "Add new column" icon. A dialog box should open.
3. User Interface Check:
- Hover over the entire dialog box.
- Expected Result: No overlapping text. No unnecessary scrollbars, icons, or lists extending beyond the boundaries. All tooltips should provide clear descriptions.
- Resize the dialog box (increase and decrease the size). The dialog box should adjust appropriately to the changes in size.
4. Add a New Column:
- Set the new column name to "New".
- Set the formula to Round the sum of HEIGHT and WEIGHT columns using autocomplete hints and columns drag-n-drop.
- Press Ok. Expected Result: A new column named "New" should be added to the dataset.
5. Retrieve Recent Activity:
- Press the "Add new column" icon again.
- Locate Recent Activities icon.
- Select the most recent activity (the column you just added).
- Expected Result: The form should autofill with the details of the last activity (i.e., column name "New" and the formula).

---
{
  "order": 2,
  "datasets": ["System:DemoFiles/demog.csv"]  
}