#### Heatmap Viewer

1. Load test data
- In TestTrack, click the **Star icon** with the tooltip "Open test data".
- This action should **open the following datasets**:
  - SPGI
  - SPGI-linked1
  - SPGI-linked2
2. Open Heatmap viewer
- Go to Tables > SPGI
- Open the **Viewers tab**, then click **Heatmap**
- Click the **Gear** icon on the viewer to open **the Context Panel**
3. Verify data properties
- In **the Context Panel**, check that table switching works correctly:     
  - Switch between SPGI, SPGI-linked1, and SPGI-linked2
  - The configuration can be saved to a Layout and reopened correctly
4. Test **custom sorting** on **categorical column**
- In the SPGI, right-click the column header **Primary Series Name**
- Select **Sort > Custom**. Move the empty value to the top of the list
- Click Apply, then OK. 
- Right-click the same column again:
  - Choose Sort > Ascending → Confirm the empty value appears first
  - Choose Sort > Descending → Confirm sorting order is reversed correctly
5. Test **Layout saving** and visualization restoration
- Open a new instance of the SPGI dataset
- Open Heatmap from the Viewers tab
- In the Property Pane:
  - Increase Max Heatmap Columns to 100
  - Uncheck Is heatmap
  - Scroll horizontally in the viewer
- Save this view as a Layout
- Switch to another layout, then come back to the saved one
- Enable Is heatmap again in the Property Pane
- **Expected Results:** 
  - The Heatmap should restore to the exact same configuration it had before you disabled Is heatmap
  - The scrollbar position should be preserved 
  - No visual differences should appear after applying the layout


---
{
  "order": 7,
  "datasets": ["System:DemoFiles/SPGI.csv","System:DemoFiles/SPGI-linked1.csv","System:DemoFiles/SPGI-linked2.csv"]
}