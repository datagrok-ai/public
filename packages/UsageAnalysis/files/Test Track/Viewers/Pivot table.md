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
* Save the layout.
* Switch to another layout and then reapply the saved layout.
* **Expected Result**: The viewer title should be saved and persist in the layout, even if it was edited directly in the viewer's header. After switching and reapplying the layout, the title should remain consistent.

Check tickets:

* [#2606](https://github.com/datagrok-ai/public/issues/2606)
* [#2198](https://github.com/datagrok-ai/public/issues/2198)
* [#2497](https://github.com/datagrok-ai/public/issues/2497)
* [#2535](https://github.com/datagrok-ai/public/issues/2535)

---
{
  "order": 12
}