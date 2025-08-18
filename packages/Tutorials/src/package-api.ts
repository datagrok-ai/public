import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function trackOverview(): Promise<any> {
    return await grok.functions.call('Tutorials:TrackOverview', {});
  }

  export async function tutorialWidget(): Promise<any> {
    return await grok.functions.call('Tutorials:TutorialWidget', {});
  }

  export async function tutorialsInit(): Promise<any> {
    return await grok.functions.call('Tutorials:TutorialsInit', {});
  }

  //Interactive demo of major Datagrok capabilities
  export async function demoApp(path: string, filter: string): Promise<any> {
    return await grok.functions.call('Tutorials:DemoApp', { path, filter });
  }

  export async function demoAppTreeBrowser(treeNode: any): Promise<any> {
    return await grok.functions.call('Tutorials:DemoAppTreeBrowser', { treeNode });
  }

  export async function demoAppWidget(): Promise<any> {
    return await grok.functions.call('Tutorials:DemoAppWidget', {});
  }

  //The File Manager is an interface that allows you to manage connections, browse and preview file content, and perform standard file and folder actions such as opening, downloading, deleting, and renaming.
  export async function filesDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:FilesDemo', {});
  }

  //Database Manager provides a hierarchical browsing interface for schemas and database objects, such as queries, tables, and table columns (if supported by the providers). You can perform various operations like adding new connections and queries, previewing data, running queries, and managing objects using context actions that are accessible through right-clicking an object.
  export async function databasesDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:DatabasesDemo', {});
  }

  //A scatter plot (also called a scatter graph, scatter chart, scattergram, or scatter diagram) is a type of plot or mathematical diagram using Cartesian coordinates to display values for typically two variables for a set of data. If the points are color-coded, you can increase the number of displayed variables to three. The data is displayed as a collection of points, each having the value of one variable determining the position on the horizontal axis and the value of the other variable determining the position on the vertical axis.
  export async function scatterPlotDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:ScatterPlotDemo', {});
  }

  //A histogram is a graphical representation of the distribution of numerical data.
  export async function histogramDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:HistogramDemo', {});
  }

  //A line chart is a simple, two-dimensional chart with an X and Y axis, each point representing a single value. The data points are joined by a line to depict a trend, usually over time.
  export async function lineChartDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:LineChartDemo', {});
  }

  //A bar chart presents grouped data as rectangular bars with lengths proportional to the values that they represent. Unlike histograms which you can apply to display the distribution of numerical data, bar charts are primarily designed for categorical values.
  export async function barChartDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:BarChartDemo', {});
  }

  //A pie chart is useful for reflecting numerical proportions. Conceptually, it is similar to a bar chart in that it represents categorical values. A pie chart shows the relative size of a given category (a slice of the pie) compared to the entire dataset (the whole pie).
  export async function pieChartDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:PieChartDemo', {});
  }

  //Trellis Charts are useful for finding the structure and patterns in complex data. A Trellis Chart is a layout of smaller charts in a grid with consistent scales. Each smaller chart represents rows that belong to a corresponding category. The grid layout looks similar to a garden trellis, hence the name Trellis Chart.
  export async function trellisPlotDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:TrellisPlotDemo', {});
  }

  //Use a Matrix Plot to assess the relationship among many pairs of columns at the same time.
  export async function matrixPlotDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:MatrixPlotDemo', {});
  }

  //Use a 3D scatter plot to plot data points on three axes to show the relationship between three variables. Each row in the data table is represented by a marker whose position depends on its values in the columns set on the X, Y, and Z axes. Additionally, you can color-code and size-code points, as well as display labels next to markers.
  export async function scatterPlot3DDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:ScatterPlot3DDemo', {});
  }

  //Unlike a scatter plot that visualizes each individual data point, the density plot splits 2D area by bins and color-codes it depending on the number of points that fall within this bin. The darker the color, the more points it contains.
  export async function densityPlotDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:DensityPlotDemo', {});
  }

  //Parallel coordinates are a common way of visualizing high-dimensional geometry and analyzing multivariate data. To show a set of points in an n-dimensional space, a backdrop is drawn consisting of n parallel lines, typically vertical and equally spaced. A point in n-dimensional space is represented as a polyline with vertices on the parallel axes; the position of the vertex on the i-th axis corresponds to the i-th coordinate of the point. This visualization is closely related to time series visualization, except that it is applied to data where the axes do not correspond to points in time, and therefore do not have a natural order. Therefore, different axis arrangements may be of interest.
  export async function pcPlotDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:PcPlotDemo', {});
  }

  //A network diagram is used to visualize graphs, where values of the specified two columns become nodes, and rows become edges. It is possible to color-code and size-code nodes and columns by choosing the aggregate function that would apply to the values that represent an edge or a Node.js.
  export async function networkDiagramDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:NetworkDiagramDemo', {});
  }

  //The box plot (a.k.a. box and whisker diagram) is a standardized way of displaying the distribution of data based on the five-number summary: minimum, first quartile, median, third quartile, and maximum.
  export async function boxPlotDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:BoxPlotDemo', {});
  }

  //Treemap displays hierarchical (tree-structured) data as nested rectangles. The branches are rectangles, then tiled with smaller rectangles representing sub-branches. Rectangles have areas proportional to a specified dimension of the data using the specified aggregation function (count by default).
  export async function treeMapDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:TreeMapDemo', {});
  }

  //Provides specified descriptive statistics for the chosen columns.
  export async function statisticsDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:StatisticsDemo', {});
  }

  //A quick way to assess correlations between all columns at once. Cells are color-coded by the Pearson correlation coefficient or Spearman's rank correlation coefficient. Histograms along the diagonal show the corresponding distribution. Hover over the cell to see the corresponding scatter plot. The grid is sortable. Select columns in the view by selecting corresponding rows.
  export async function correlationPlotDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:CorrelationPlotDemo', {});
  }

  //Calendar lets you analyze longitudinal data. It needs at least one column of type DateTime.
  export async function calendarDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:CalendarDemo', {});
  }

  //A grid table contains a set of data that is structured in rows and columns. It allows the user to scroll in both directions and can handle large numbers of items and columns.
  export async function gridDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:GridDemo', {});
  }

  //Use this viewer to host any text, arbitrary HTML content, or markdown-formatted text. In most cases, the viewer will auto-detect content type. Use the "mode" property to explicitly specify it.
  export async function markupDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:MarkupDemo', {});
  }

  //Visualizes rows as a collection of forms that are positioned as tiles.
  export async function tileViewerDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:TileViewerDemo', {});
  }

  //Form allows you to customize the appearance of the row by manually positioning the fields, and adding other visual elements, such as pictures or panels. A form can be used either as a stand-alone viewer or as a row template of the Tile Viewer.
  export async function formDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:FormDemo', {});
  }

  //A pivot table is a table of grouped values that aggregates the individual items of a more extensive table within one or more discrete categories. This summary might include sums, averages, or other statistics, which the pivot table groups together using a chosen aggregation function applied to the grouped values.
  export async function pivotTableDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:PivotTableDemo', {});
  }

  //Filter is a set of controls for quick filtering, selection, and visual assessment of column values.
  export async function filtersDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:FiltersDemo', {});
  }

  //Table linking is based on the key columns. The link type determines which records should be synchronized between the datasets (current record, filter, selection, and mouse-over record).
  export async function tableLinkingDemo(): Promise<any> {
    return await grok.functions.call('Tutorials:TableLinkingDemo', {});
  }
}
