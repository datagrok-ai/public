#### 1. Scatterplot: Color/Marker settings and dynamic legends

1. Open SPGI
2. Add a scatterplot
3. Set **Color** and **Marker** to `Series` — legend must be combined
4. Check color picker visibility, change some colors
5. Save and apply the layout
6. Save and open the project — changes should persist
7. On the plot, add new **Color** columns:
   - `if(${Stereo Category}=='S_UNKN', null, ${Average Mass})` → linear legend appears
   - `if(${Stereo Category}=='S_UNKN', null, ${Series})` → categorical legend appears
8. Set **Color** to `ID`, **Marker** to `Core` — verify legend is displayed
9. Close All

#### 2. Scatterplot: Legend update on axis change

1. Open SPGI
2. Add new columns:
   - `col1: if(${Stereo Category}!='S_UNKN', null, ${Average Mass})`
   - `col2: if(${Stereo Category}=='S_UNKN', null, ${Average Mass})`
3. Add a scatterplot
4. Set **X axis** to `col1`, **Color** to `Stereo Category`
5. Change **X axis** to `col2`
6. Verify legend categories update according to data
7. Test zooming/filtering — verify legend stays consistent
8. Close All

#### 3. Scatterplot: In-viewer filtering

1. Open SPGI
2. Add a scatterplot
3. Set **Marker** to `Stereo category`
4. Apply in-viewer filter: `${Stereo Category} in ["R_ONE", "S_UNKN"]`
5. Add another scatterplot with the same filter and marker
6. Save and apply the layout — verify filtered legends
7. Close All

#### 4. Scatterplot: Filtering via Filter Panel

1. Open SPGI
2. Add a scatterplot
3. Set **X/Y axes** to `Chemical space X/Y`
4. Set **Color** to `Primary scaffold name`, **Marker** to `Stereo category`
5. In **Filter Panel**, filter `Primary scaffold name` — verify data and legend update
6. Click 'R_ONE' in the scatterplot legend — verify correct additional filtering, no reappearing filtered-out points
7. Close All

#### 5. Color coding from grid

1. Open SPGI
2. Add a scatterplot, box plot, and PC plot
3. Set **Color** to `Chemical Space X`
4. In the grid, enable linear color coding for `Chemical Space X` — check legends on both viewers
5. Change color schema, invert, apply to text — check legends
6. Save and apply layout — verify color changes
7. Change `Chemical Space X` coding to categorical, modify colors — verify legends updates
8. Save and open the project — verify changes persist

---
{
  "order": 5,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}