### Pivot table

1. Open the **demog** dataset.
2. Click the **Add viewer** icon, find **Pivot table** viewer and press on it. 
* **Expected result**: **Pivot table** viewer opens without errors in console. 
3. Close previously opened viewer. On the **Viewers tab**, click **Pivot table** icon. 
* **Expected result**: **Pivot table** viewer opens without errors in console. 
4. Interact with all elements directly on the viewer.
* **Pivot** and **Group by** fields should be availiable for all columns (ID column).
* **Expected Result**: Changes should be reflected quickly and without errors. 
5. On the **Pivot table** viewer, click the **Gear** icon. The **Property Pane** opens.
6. Modify various properties in **Property Pane**:
* **Expected Result**: Changes to the properties should be reflected in the Pivot table without errors, and the viewer should be updated accordingly.
7. Adding and Saving a **Title** in the Pivot table Viewer:
* Open the Pivot Table viewer on SPGI dataset.
* Add a title to the **Pivot Table** using the **Property Pane**.
* Edit the title directly in the viewer's header.
* Add any coloring settings for the Pivot table (e.g. column coloring).
* Save the layout.
* Switch to another layout and then reapply the saved layout.
* **Expected Result**: The viewer title should be saved and persist in the layout, even if it was edited directly in the viewer's header. After switching and reapplying the layout, the title should remain consistent. Coloring settings are present.
8. Open the **demog** dataset with **Pivot table** viewer on it. 
* Apply colouring for some columns.
* **Property Pane** change row source form 'Filtered' to 'Selected' and than change back to 'Filtered'.
* **Expected Result**: Row source changes shouldn't reset coloring setting. 
9. Editing the viewer properties (changes done in the viewer itself and in its properties panel):
* Columns selected in viewer should be reflected in properties panel
* Columns selected in the panel should be applied to the viewer

---
{
  "order": 12,
  "datasets": ["System:DemoFiles/demog.csv"]
}