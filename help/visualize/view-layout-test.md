<!-- TITLE: Tests: View Layout -->
<!-- SUBTITLE: -->

# Tests: View Layout

[View Layout](view-layout.md) contains relative positions of [viewers](../visualize/viewers.md) in a [table view](../overview/table-view.md),
along with the viewers' properties. By separating layouts from the actual data displayed, we now can
save current layout (**View | Layout | Save to Gallery**) and later apply it to a different dataset
(**View | Layout | Open Gallery**). 

To clone current view, either do **View | Layout | Clone**, or click on the plus sign on the view header strip, 
and choose **Clone**.

## Testing scenario

1. Open *"demog"* table

1. Clone table *"demog"* view from menu **View | Layout | Clone** or from menu **+** on tabbed ribbon or use keyboard ```Ctrl + Shift + V```
   * New view created for *"demog"* table 
   * New view is called *"demog copy"*

1. Remove view *"demog copy"* from [Workspace](../overview/workspace.md)

1. Add each one available [Viewer](../visualize/viewers.md) for *"demog"* table 

1. Change style properties for added viewers. (Colors of viewer elements, size, scale, font, etc.)

1. Clone table *"demog"* view from menu **View | Layout | Clone** or from menu **+** on tabbed ribbon or use keyboard ```Ctrl + Shift + V```
   * New view created for *"demog"* table 
   * New view is called *"demog copy"* 
   * On *"demog copy"* view there are added viewers with modified style properties.

1. Remove view *"demog copy"* from [Workspace](../overview/workspace.md)

1. Clear *"demog"* view from menu **View | Layout | Clear** or **View | Reset** from context menu to *demog* tab or use keyboard ```Ctrl + Shift + R```
   * *"demog"* table is in its original form is presented on *"demog"* view
   
1. Change row sorting by *Age* column to decrease
   * Table is sorted by values ​​of *Age* column descending
   
1. Enable Color Coding for *Sex* column
   * Values ​​in *Sex* column are highlighted in different colors

1. Set filter for values ​​only equals ```Asian``` in "Race" column
   * Filter ```"Race" = Asian``` in  applied *"demog"* table
   
1. Swap columns *Age* and *Sex* in table

1. Clone table *"demog"* view from menu **View | Layout | Clone** or from menu **+** on tabbed ribbon or use keyboard ```Ctrl + Shift + V```
   * New view created for *"demog"* table 
   * New view is called *"demog copy"*
   * Row order changes, color coding, filtering and column ordering applied to *"demog copy"* view.

1. Clear *"demog"* view from menu **View | Layout | Clear** or **View | Reset** from context menu to *demog* tab or use keyboard ```Ctrl + Shift + R```
   * *"demog"* table is in its original form is presented on *"demog"* view
  
1. Add [Scatter Plot](../visualize/viewers/scatter-plot.md) on *"demog"* view

1. Change Scatter Plot](../visualize/viewers/scatter-plot) style settings. (Colors of viewer elements, size, scale, font, etc.)

1. Swap columns *Age* and *Sex* in table

1. Change row sorting by *Age* column to decrease
   * Table is sorted by values ​​of *Age* column descending
   
1. Save received view to Layouts Gallery from menu **View | Layout | Save to Gallery** or use ```Ctrl + S``` or use 
   **SAVE** button in "Layouts" tab on [Toolbox](../overview/navigation.md#toolbox)

1. Open Layouts Gallery from **View | Layout | Open Gallery** or use ```Ctrl + L```
   * In Layouts Gallery there is saved early view

1. Restart **Grok** platform and open *"demog"* table in new session

1. Open "Layouts" tab on [Toolbox](../overview/navigation.md#toolbox)
   * Saved layout in step 19 is offered in "Layouts" tab.
   
1. Open Layouts Gallery from **View | Layout | Open Gallery** or use ```Ctrl + L```
   * In Layouts Gallery there is saved view from step 19
   
1. Apply layout from step 19 to *"demog"* view
   * [Scatter Plot](../visualize/viewers/scatter-plot.md) with style settings added to *"demog"* view
   * Row and columns order changes to *"demog"* view.
   
1. Delete layout from step 19 from Layouts Gallery from its context menu   
   

See also:
* [View Layout](view-layout.md)
* [Table view](../overview/table-view.md)
* [Self-learning platform](../learn/self-learning-platform.md)
