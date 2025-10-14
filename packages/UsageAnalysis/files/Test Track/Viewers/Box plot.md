### Box plot viewer

1. Open SPGI, SPGI-linked1, SPGI-linked2.
2. Go to **Tables > SPGI.**
3. On the **Viewers tab,** click **Box plot.**
4. On the **Box plot** viewer, click the **Gear** icon. The **Property Pane** opens.
5. Go to the **Data** info panel and check all the properties, including the following:
    - Tables switching (SPGI, SPGI-linked2, SPGI-linked1)
    - Set **Filter** to ${Average Mass} > 225
    - Set **Color** to Link Column 1
    - Change Bin Color Aggr Type
    - Save to Layout. Check
6. Axes:
    - Change axes
    - Check the **Invert Y Axis** checkbox
    - Save to Layout. Check
7. Tooltips testing:
    - Right-click the box plot and check all options on the **Tooltip** tab.
    - Save to Layout. Check
    - Go to Property Pane > Tooltip. Check all options.
    - Save to Layout. Check
8. Coloring: (checked in steps 5.c- 5.e)
    - Set **Color** to Link Column 1
    - Change Bin Color Aggr Type
    - Save to Layout. Check 
    - Close all.
9. Visualization zoom:
    - The SPGI_v2 dataset is opened and a Box Plot viewer is added.
    - Adjust some zoom in the Box Plot viewer.
    - Save the current project. 
    - Reopen the saved project.
    - Verify that the zoom level in the Box Plot viewer is preserved.
    - Adjust the zoom level again (if needed).
    - Save the layout and reload this layout.
    - Verify that the zoom level is not preserved when reopening the layout.


---
{
  "order": 6,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}