import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Tutorials
//meta.role: app
//top-menu: Help | Tutorials @Toolbox Help | Tutorials
export function trackOverview() : void {
  PackageFunctions.trackOverview();
}

//output: widget result
export function tutorialWidget() : any {
  return PackageFunctions.tutorialWidget();
}

//meta.role: autostart
export function tutorialsAutostart() : void {
  PackageFunctions.tutorialsAutostart();
}

//meta.role: init
export function tutorialsInit() : void {
  PackageFunctions.tutorialsInit();
}

//name: Demo
//description: Interactive demo of major Datagrok capabilities
//input: string path { meta.url: true; optional: true }
//input: string filter { optional: true }
//output: view result
//meta.browseOnly: true
//meta.role: app
//meta.icon: images/icons/demoapp-icon.png
export function demoApp(path?: string, filter?: string) : any {
  return PackageFunctions.demoApp(path, filter);
}

//input: dynamic treeNode 
export async function demoAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.demoAppTreeBrowser(treeNode);
}

//output: widget result
export function demoAppWidget() : any {
  return PackageFunctions.demoAppWidget();
}

//output: string result
export function getDemoAppHierarchy() : string {
  return PackageFunctions.getDemoAppHierarchy();
}

//name: filesDemo
//description: The File Manager is an interface that allows you to manage connections, browse and preview file content, and perform standard file and folder actions such as opening, downloading, deleting, and renaming.
//meta.demoPath: Data Access | Files
export function _filesDemo() : void {
  PackageFunctions._filesDemo();
}

//name: databasesDemo
//description: Database Manager provides a hierarchical browsing interface for schemas and database objects, such as queries, tables, and table columns (if supported by the providers). You can perform various operations like adding new connections and queries, previewing data, running queries, and managing objects using context actions that are accessible through right-clicking an object.
//meta.demoPath: Data Access | Databases
export async function _databasesDemo() : Promise<void> {
  await PackageFunctions._databasesDemo();
}

//name: scatterPlotDemo
//description: A scatter plot (also called a scatter graph, scatter chart, scattergram, or scatter diagram) is a type of plot or mathematical diagram using Cartesian coordinates to display values for typically two variables for a set of data. If the points are color-coded, you can increase the number of displayed variables to three. The data is displayed as a collection of points, each having the value of one variable determining the position on the horizontal axis and the value of the other variable determining the position on the vertical axis.
//meta.demoPath: Visualization | General | Scatter Plot
export async function _scatterPlotDemo() : Promise<void> {
  await PackageFunctions._scatterPlotDemo();
}

//name: histogramDemo
//description: A histogram is a graphical representation of the distribution of numerical data.
//meta.demoPath: Visualization | General | Histogram
export async function _histogramDemo() : Promise<void> {
  await PackageFunctions._histogramDemo();
}

//name: lineChartDemo
//description: A line chart is a simple, two-dimensional chart with an X and Y axis, each point representing a single value. The data points are joined by a line to depict a trend, usually over time.
//meta.demoPath: Visualization | General | Line Chart
export async function _lineChartDemo() : Promise<void> {
  await PackageFunctions._lineChartDemo();
}

//name: barChartDemo
//description: A bar chart presents grouped data as rectangular bars with lengths proportional to the values that they represent. Unlike histograms which you can apply to display the distribution of numerical data, bar charts are primarily designed for categorical values.
//meta.demoPath: Visualization | General | Bar Chart
export async function _barChartDemo() : Promise<void> {
  await PackageFunctions._barChartDemo();
}

//name: pieChartDemo
//description: A pie chart is useful for reflecting numerical proportions. Conceptually, it is similar to a bar chart in that it represents categorical values. A pie chart shows the relative size of a given category (a slice of the pie) compared to the entire dataset (the whole pie).
//meta.demoPath: Visualization | General | Pie Chart
export async function _pieChartDemo() : Promise<void> {
  await PackageFunctions._pieChartDemo();
}

//name: trellisPlotDemo
//description: Trellis Charts are useful for finding the structure and patterns in complex data. A Trellis Chart is a layout of smaller charts in a grid with consistent scales. Each smaller chart represents rows that belong to a corresponding category. The grid layout looks similar to a garden trellis, hence the name Trellis Chart.
//meta.demoPath: Visualization | Data Separation | Trellis Plot
export async function _trellisPlotDemo() : Promise<void> {
  await PackageFunctions._trellisPlotDemo();
}

//name: matrixPlotDemo
//description: Use a Matrix Plot to assess the relationship among many pairs of columns at the same time.
//meta.demoPath: Visualization | Data Separation | Matrix Plot
export async function _matrixPlotDemo() : Promise<void> {
  await PackageFunctions._matrixPlotDemo();
}

//name: scatterPlot3DDemo
//description: Use a 3D scatter plot to plot data points on three axes to show the relationship between three variables. Each row in the data table is represented by a marker whose position depends on its values in the columns set on the X, Y, and Z axes. Additionally, you can color-code and size-code points, as well as display labels next to markers.
//meta.demoPath: Visualization | General | 3D Scatter Plot
export async function _scatterPlot3DDemo() : Promise<void> {
  await PackageFunctions._scatterPlot3DDemo();
}

//name: densityPlotDemo
//description: Unlike a scatter plot that visualizes each individual data point, the density plot splits 2D area by bins and color-codes it depending on the number of points that fall within this bin. The darker the color, the more points it contains.
//meta.demoPath: Visualization | General | Density Plot
export async function _densityPlotDemo() : Promise<void> {
  await PackageFunctions._densityPlotDemo();
}

//name: pcPlotDemo
//description: Parallel coordinates are a common way of visualizing high-dimensional geometry and analyzing multivariate data. To show a set of points in an n-dimensional space, a backdrop is drawn consisting of n parallel lines, typically vertical and equally spaced. A point in n-dimensional space is represented as a polyline with vertices on the parallel axes; the position of the vertex on the i-th axis corresponds to the i-th coordinate of the point. This visualization is closely related to time series visualization, except that it is applied to data where the axes do not correspond to points in time, and therefore do not have a natural order. Therefore, different axis arrangements may be of interest.
//meta.demoPath: Visualization | Statistical | PC Plot
export async function _pcPlotDemo() : Promise<void> {
  await PackageFunctions._pcPlotDemo();
}

//name: networkDiagramDemo
//description: A network diagram is used to visualize graphs, where values of the specified two columns become nodes, and rows become edges. It is possible to color-code and size-code nodes and columns by choosing the aggregate function that would apply to the values that represent an edge or a Node.js.
//meta.demoPath: Visualization | Data Flow and Hierarchy | Network Diagram
export async function _networkDiagramDemo() : Promise<void> {
  await PackageFunctions._networkDiagramDemo();
}

//name: boxPlotDemo
//description: The box plot (a.k.a. box and whisker diagram) is a standardized way of displaying the distribution of data based on the five-number summary: minimum, first quartile, median, third quartile, and maximum.
//meta.demoPath: Visualization | Statistical | Box Plot
export async function _boxPlotDemo() : Promise<void> {
  await PackageFunctions._boxPlotDemo();
}

//name: treeMapDemo
//description: Treemap displays hierarchical (tree-structured) data as nested rectangles. The branches are rectangles, then tiled with smaller rectangles representing sub-branches. Rectangles have areas proportional to a specified dimension of the data using the specified aggregation function (count by default).
//meta.demoPath: Visualization | Data Flow and Hierarchy | Tree Map
export async function _treeMapDemo() : Promise<void> {
  await PackageFunctions._treeMapDemo();
}

//name: statisticsDemo
//description: Provides specified descriptive statistics for the chosen columns.
//meta.demoPath: Visualization | Statistical | Statistics
export async function _statisticsDemo() : Promise<void> {
  await PackageFunctions._statisticsDemo();
}

//name: correlationPlotDemo
//description: A quick way to assess correlations between all columns at once. Cells are color-coded by the Pearson correlation coefficient or Spearman's rank correlation coefficient. Histograms along the diagonal show the corresponding distribution. Hover over the cell to see the corresponding scatter plot. The grid is sortable. Select columns in the view by selecting corresponding rows.
//meta.demoPath: Visualization | Statistical |Correlation Plot
export async function _correlationPlotDemo() : Promise<void> {
  await PackageFunctions._correlationPlotDemo();
}

//description: Calendar lets you analyze longitudinal data. It needs at least one column of type DateTime.
//meta.demoPath: Visualization | Time and Date | Calendar
export async function calendarDemo() : Promise<void> {
  await PackageFunctions.calendarDemo();
}

//name: gridDemo
//description: A grid table contains a set of data that is structured in rows and columns. It allows the user to scroll in both directions and can handle large numbers of items and columns.
//meta.demoPath: Visualization | Input and Edit | Grid
export async function _gridDemo() : Promise<void> {
  await PackageFunctions._gridDemo();
}

//name: markupDemo
//description: Use this viewer to host any text, arbitrary HTML content, or markdown-formatted text. In most cases, the viewer will auto-detect content type. Use the 'mode' property to explicitly specify it.
//meta.demoPath: Visualization | General | Markup
export async function _markupDemo() : Promise<void> {
  await PackageFunctions._markupDemo();
}

//name: tileViewerDemo
//description: Visualizes rows as a collection of forms that are positioned as tiles.
//meta.demoPath: Visualization | General | Tile Viewer
export async function _tileViewerDemo() : Promise<void> {
  await PackageFunctions._tileViewerDemo();
}

//name: formDemo
//description: Form allows you to customize the appearance of the row by manually positioning the fields, and adding other visual elements, such as pictures or panels. A form can be used either as a stand-alone viewer or as a row template of the Tile Viewer.
//meta.demoPath: Visualization | Input and Edit | Form
export async function _formDemo() : Promise<void> {
  await PackageFunctions._formDemo();
}

//name: pivotTableDemo
//description: A pivot table is a table of grouped values that aggregates the individual items of a more extensive table within one or more discrete categories. This summary might include sums, averages, or other statistics, which the pivot table groups together using a chosen aggregation function applied to the grouped values.
//meta.demoPath: Visualization | Statistical | Pivot Table
export async function _pivotTableDemo() : Promise<void> {
  await PackageFunctions._pivotTableDemo();
}

//name: filtersDemo
//description: Filter is a set of controls for quick filtering, selection, and visual assessment of column values.
//meta.demoPath: Visualization | General | Filters
export async function _filtersDemo() : Promise<void> {
  await PackageFunctions._filtersDemo();
}

//name: tableLinkingDemo
//description: Table linking is based on the key columns. The link type determines which records should be synchronized between the datasets (current record, filter, selection, and mouse-over record).
//meta.demoPath: Data Access | Table Linking
export async function _tableLinkingDemo() : Promise<void> {
  await PackageFunctions._tableLinkingDemo();
}
