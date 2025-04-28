#### Tile viewer

1. Multiple DataFrames Handling
- Open **SPGI, SPGI-linked1, and SPGI-linked2** datasets.
- Go to Tables > SPGI and add a **Tile viewer**.
- Open the **Context Panel > Data > Table**.
- Switch between the tables (SPGI-linked2, SPGI-linked1, SPGI) and verify the viewer updates correctly.
2. Building from List of Columns:
- Right-click the viewer and select **Edit form** from the list.
- Edit the viewer configuration and **apply the changes**.
- Save the **layout** and verify it is saved correctly.
3. **Edit form with Table Change**:
- Open the **demog** dataset.
- Switch back to SPGI and open the **Edit form** dialog.
- Change the data source to demog, press **Reset**, and select some columns.
- Click **Close and apply** button.
- In Tile viewer properties, manually switch the table to demog.
- **Expected Result**: Tile viewer should now display demog data without errors or red highlights.
4. Main menu (Hamburger Menu) testing
- Open the Hamburger menu and check all options.
- Pay attention: navigate to Properties > Data > Lanes. Hovering over all options should work smoothly without page freezing.
5. Adding and Editing a **Calculated Column**
- Open the SPGI dataset.
- **Create a calculated column**, for example: ${Average Mass} + 5, and name it 'cc'.
- Go to the Hamburger menu > **Edit form...** > press Reset.
- Ensure SPGI is selected. Add the new calculated column 'cc' to the viewer.
- Click **Close and apply** button.
- Click the 'cc' column header. In the Context Pane, find the **Formula editing input**, change the formula, and press **Apply**.
- **Expected Result**: After applying the new formula, both the column values and the Tile viewer should update immediately.
6. Added **Calculated Column with Filters**:
- Apply filters for the 'cc' column, Stereo category, and other columns.
- Expected Result: Viewer updates correctly after formula changes, even when filters are applied.


---
{
  "order": 10
}