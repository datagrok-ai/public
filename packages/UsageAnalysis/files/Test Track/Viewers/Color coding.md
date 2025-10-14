1. Open "demog" table
2. For **Age** column:
    1. Enable linear color-coding.
    2. Change color-coding type to conditional.
    3. Change values/colors.
3. For **Sex** column:
    1. Enable categorical color-coding.
    2. Change colors for categories M and F.
4. For **Control** column:
    1. Enable categorical color-coding.
5. For **Started** column:
    1. Enable linear color coding.
    2. Apply new color scheme.
6. Turn on color coding for entire table from **Context menu | Grid Color Coding | All**
7. Change color scheme from Context Menu | Grid Color Coding | Color scheme for table
8. Turn off color coding for table from **Context menu | Grid Color Coding | None**
9. Go back to the color coding state from step 5 from **Context menu | Grid Color Coding | Auto**
10. For **Age**, **Sex**, and  **Started** columns, turn off the color coding from **Context Menu** | **Color coding** | **Off**.
12. For **Age**,**Sex**, and **Started** columns, turn on color coding (it must apply custom colors which were set earlier)
13. Add a random viewer to the table and save layout of view.
14. Reopen demog table and apply the saved layout, close the viewer => you must see the same table state which was after step.
15. Copy **Race** column (Race_copy), apply categorical color coding to it
16. Pick up coloring from Race column through **Context menu** | **Color** **Coding** | **Pick Up Coloring**
17. Apply coloring to **Race_copy** column from **Context menu** | **Color Coding** | **Apply Coloring**
18. Pick up coloring from **Started** column and apply it to **Heigth** column

19. Linear Color-Coding Scheme – Define Colors (Min, Max, Inflection Point)
- Load the SPGI_v2 dataset.
- Locate the 'CAST Idea ID' column.
- Apply Linear color-coding to this column.
- Right-click the column header. Navigate to Color Coding > Edit...
- Click Edit scheme...
- In the Edit Color Scheme dialog:
  - Verify that the 'Use absolute values' checkbox has a tooltip on hover.
  - When the checkbox is checked, confirm that it enables defining absolute values.
  - Modify colors (min, max, inflection point) and confirm that changes are reflected live in the color-coding scheme preview and in the grid.
  - Clicking Cancel reverts any color changes.
  - Add or delete color range entries using + / - buttons and confirm that the scheme updates accordingly.
- Click OK to apply changes and close the dialog.
- Invert the color scheme using the arrows icon next to the scheme drop-down.


---
{
  "order": 1
}