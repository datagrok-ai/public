1. Open SPGI
2. Open the **Filter Panel**

#### 1. Apply filters in the Filter Panel

1. For the Structure filter add `c1ccccc1` - (~32 rows filtered)
2. Hover over the Filter Panel and scroll down to find the **Stereo Category** filter
3. On the Stereo Category filter, uncheck `S_UNKN` category — the filter applies (rows without S_UNKN shown)
4. For the **Average Mass** filter set max value to `400`. The filter applies (~4 rows shown).

#### 2. Add Scaffold Tree Filter

1. On the Filter Panel, add a  Scaffold Tree Filter for the Structure column
2. Add the `Cc1ccccc1` scaffold
3. Click the checkbox next to the scaffold — the filter applies and the row count changes

#### 3. Scatter Plot

1. Add a **Scatter Plot** viewer
2. Zoom in to filter the dataset - check the filtering considers the scatter plot and the filters form the Filter Panel
3. Double-click the Scatter Plot — verify the row count returns to the previous filtered value
4. Close the Scatter Plot viewer

#### 4. Bar Chart

1. Add a **Bar Chart** viewer
2. Right-click the bar chart, from the context menu select On Click > Filter
3. Click the S_ACHIR bar in the Bar Chart — verify the data is filtered 
4. Click the white space on the bar chart viewer — verify the row count returns to the previous filtered value
5. Close the Bar Chart viewer

#### 5. Histogram

1. Add a **Histogram** viewer 
2. Use the range slider under the histogram viewer to filter the dataset — verify the row count update to show only the rows within the selected range
3. Double-click the histogram viewer on the white space  — verify the row count returns to the previous filtered value
4. Close the Histogram viewer

#### 6. PC Plot

1. Add a **PC Plot** (Parallel Coordinates Plot) viewer from the Toolbox
2. Drag the range selector on any axis to filter the dataset — verify the row count update
3. Double-click the white space on the PC plot viewer — verify the row count returns to the previous filtered value
4. Close the PC Plot viewer

#### 7. Trellis Plot

1. Add a **Trellis Plot** viewer from the Toolbox
2. Right-click the trellis plot viewer, from the context menu selest On Click > Filter
3. Click on a cell in the Trellis Plot — verify the grid and Filter Panel row count update to show only the rows corresponding to the clicked cell. If the clicked cell is empty, row count expected to be 0
4. Esc to return to the previous filtered state
5. Close the Trellis Plot viewer

#### 8. Pie Chart

1. Add a **Pie Chart** viewer
2. Right-click the bar chart, from the context menu select On Click > Filter
3. Click any segment in the pie chart — verify the data is filtered 
4. Click the white space on the pie chart viewer — verify the row count returns to the previous filtered value
5. Close the pie chart viewer

#### 9. Reset and cleanup

1. Click the **Reset** icon in the Filter Panel header and confirm all filters are cleared and all rows are shown
2. Close All
---
{
"order": 4,
"datasets": ["System:DemoFiles/SPGI.csv"]
}
