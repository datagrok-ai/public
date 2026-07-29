# Scatter plot manual checklist

Human-only checks for the Scatter plot; automated coverage lives in the other
scenarios in this folder — axes and encoding, selection and zoom, zoom-driven
filtering, the legend, regression and formula lines, labels and tooltip, plus
the secondary settings in the main scatter plot scenario. What remains here
needs a human eye or a real second table: switching the viewer between several
open tables, chemical structures used as markers and as an axis, and keeping
the plot's colors in step with the grid's color coding.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add a scatter plot by clicking the Scatter plot icon in the Toolbox > Viewers section

## Multiple tables and switching Data > Table

1. Open a second table alongside demog (for example
   `System:AppData/Chem/tests/spgi-100.csv`) so two tables are loaded
2. Go to the **Context Panel > Data** and switch **Table** to the other table
   -- the viewer rebinds to it and the axes are re-picked for the new columns
3. Switch **Table** back to the original table -- the axes are re-picked again
   and the plot renders the original data
4. Note that the table names shown in the picker are the in-memory table names,
   which do not have to match the file names
5. In the **Context Panel > Data**, set **Filter** to an expression over a
   column of the current table -- verify only the matching points remain on the
   plot
6. Clear the **Filter**

## Structures as markers and as an axis

1. Close all and open a chemical dataset that has a structure column (for
   example `System:AppData/Chem/tests/spgi-100.csv`)
2. Add a scatter plot
3. Go to the **Context Panel > Marker** and set **Markers** to the structure
   column -- verify structure markers render on the plot
4. Set **Markers** to a datetime column -- verify the markers change and the
   legend stays consistent
5. Change the **Markers Map** values for the datetime column -- verify the
   legend and the plot stay consistent with each other
6. Set **Markers** to a categorical column -- verify standard shape markers
   render
7. On the viewer, click the X column selector and set X to the structure column
   -- verify structures render along the X axis
8. Save the layout, close the scatter plot, apply the saved layout -- verify
   the structure axis is restored
9. Delete the saved layout

## Grid color coding synchronized with the plot

1. Go to the grid, right-click the **SEX** column header and choose
   **Color coding > Categorical**
2. On the scatter plot, click the Color selector and set Color to SEX --
   verify the point colors match the grid's categorical scheme
3. Go to the grid, right-click the **AGE** column header and choose
   **Color coding > Linear**, then **Edit** -- change the color scheme and
   invert it
4. On the scatter plot, click the Color selector and set Color to AGE --
   verify the plot uses the same inverted scheme as the grid
5. Go to the grid, right-click the **AGE** column header and choose
   **Color coding > Edit** -- change the scheme again
6. Verify the scatter plot updates to match
7. Remove the color coding from both grid columns and clear the plot's Color
   selector

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
