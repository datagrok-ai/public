# PowerGrid

PowerGrid contains popular spreadsheet extensions, such as [sparklines](#sparklines)
and support for [images](#images).

## Sparklines

Sparklines are a very small charts, typically without axes or coordinates. The main goal of a sparkline is 
to provide a graphical representation of a data object that allows not only to get an impression of the object 
as a whole but also visually instantly compare objects with each other. To form a representation, any quantitative 
features of objects can be used, which are displayed as the position or size of a graphical element.
Datagrok has a few ready-to-use sparkline chart types 'sparklines' 'bar chart' 'pie chart' and 'radar chart'. 
You can customize which columns / objects' features will be used to form the graph.

![](../../help/develop/how-to/custom-cell-renderers-sparklines-and-settings.gif "Sparklines and settings")

## Images

Datagrok supports two types of images in cells: embedded and linked.

### Embedded images

To programmatically add an image column to the dataframe, use the
following code:  ```dataFrame.columns.addNewBytes('my image')```.

To add an image, double-click on the cell and choose the file.
Supported extensions are `.jpg`, `.png`, and `.jpeg`. 

![](../../help/develop/how-to/binary-cell-renderer.gif "adding image")


See also:

* [Viewers](https://datagrok.ai/help/visualize/viewers)
* [Develop custom cell renderers](https://datagrok.ai/help/develop/how-to/custom-cell-renderers)
* [Packages](https://datagrok.ai/help/develop/develop#packages)
