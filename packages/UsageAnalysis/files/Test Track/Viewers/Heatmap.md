1. Open SPGI, SPGI-linked1, SPGI-linked2 by pressing 'star' icon in TestTrack (with tooltip "Open test data").
2. Go to **Tables > SPGI**.
3. On the **Viewers tab**, click **Heatmap**.
4. On the **heatmap** viewer, click the **Gear** icon. The **Property Pane** opens.
5. Go to the **Data** info panel and check all the properties, including the following:
    1. Tables switching (SPGI, SPGI-linked2, SPGI-linked1).
    2. Save to Layout. Check
6. Custom sorting of categorical columns:
    1. Go to the grid.
    2. Right-click the Primary Series Name column’s header and select **Sort > Custom** from the context menu.
    3. Move the empty value to the first place. Click Apply, OK
    4. Right-click the Primary Series Name column’s header and select **Sort > Ascending** (empty value should be in the first place), **Descending**.
7. Layouts saving & visualization:
* Open new **SPGI** dataset. On the **Viewers tab**, click **Heatmap**.
* In the **Property Pane**:
  * increase 'Max Heatmap Columns' to 100.
  * switch off 'Is heatmap' option.
  * move scrollbar position
* Save and re-apply layout.
* Switch on 'Is heatmap' option in the **Property Pane**. 
* **Expected Results**: 
  * After applying the layout, the Heatmap visualization should be identical to the view before the 'Is heatmap' option was switched off. 
  * The scrollbar position should be saved for layout
  * There should be no visible differences in the Heatmap display.

---
{
  "order": 7,
  "datasets": ["System:DemoFiles/SPGI.csv","System:DemoFiles/SPGI-linked1.csv","System:DemoFiles/SPGI-linked2.csv"]
}