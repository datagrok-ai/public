import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {TutorialRunner, TutorialSubstitute} from './tutorial-runner';
import {chem} from './tracks/chem';
import {eda} from './tracks/eda';
import {da} from './tracks/data-access';
import {ml} from './tracks/ml';
import {dataTransformation} from './tracks/transform';
import {TutorialWidget} from './widget';
import '../css/tutorial.css';
import {Track} from '@datagrok-libraries/tutorials/src/track';
import {DemoView} from './demo-app/demo-app';
import {viewerDemo} from './demo-app/platform-viewers-demo';
import { DemoAppWidget } from './demo-app/widget';

export const _package = new DG.Package();

const tracks: Track[] = [];
const PATH_START_INDEX: number = 4;

//name: Tutorials
//tags: app
//top-menu: Help | Tutorials @Toolbox Help | Tutorials
export function trackOverview() {
  const tutorialRunners = tracks.map((track) => new TutorialRunner(track));
  const root = ui.div([
    ...tutorialRunners.map((runner) => runner.root),
    ui.panel([], {id: 'tutorial-child-node', style: {paddingTop: '10px'}}),
  ], 'tutorials-root');

  grok.shell.dockManager.dock(root, DG.DOCK_TYPE.RIGHT, null, 'Tutorials', 0.3);
  setPath(window.location.pathname, tutorialRunners);
}

//output: widget tutorial
export function tutorialWidget(): DG.Widget {
  return new TutorialWidget(...tracks.map((track) => new TutorialRunner(track)));
}

//tags: init
export async function tutorialsInit() {
  const properties = await _package.getProperties();
  setProperties(properties as unknown as { [propertyName: string]: boolean });

  const tutorialFuncs = DG.Func.find({tags: ['tutorial']});
  const trackFuncs = DG.Func.find({tags: ['track']});
  const defaultIcon = `${_package.webRoot}package.png`;
  const tutorialIconPaths: { [packageName: string]: { [tutorialName: string]: {
    tutorial: TutorialSubstitute,
    imageUrl: string,
  } } } = {};

  for (const track of tracks) {
    for (const tutorial of track.tutorials) {
      if (tutorial.imageUrl === '')
        tutorial.imageUrl = `${_package.webRoot}images/${tutorial.name.toLowerCase().replace(/ /g, '-')}.png`;
    }
  }

  for (const func of trackFuncs) {
    const t = tracks.find((t) => t.name === func.options['name']);
    if (!t)
      tracks.push(new Track(func.options['name'], [], func.helpUrl ?? ''));
    else
      console.error(`Tutorials: Couldn't add a track. Track ${func.options['name']} already exists.`);
  }

  for (const func of tutorialFuncs) {
    const trackName = func.options['track'] ?? func.package.friendlyName;
    const tutorial = new TutorialSubstitute(func.options['name'], func.description, defaultIcon, func);

    if (func.options['icon']) {
      tutorialIconPaths[func.package.name] = {[func.options['name']]: {
        tutorial,
        imageUrl: func.options['icon'],
      }};
    }

    const track = tracks.find((t) => t.name === trackName);
    if (track) {
      const t = track.tutorials.find((t) => t.name === tutorial.name);
      if (!t)
        track.tutorials.push(tutorial);
      else {
        console.error(`Tutorials: Couldn't add a tutorial to track ${track.name}. ` +
          `Tutorial ${tutorial.name} already exists.`);
      }
    } else
      tracks.push(new Track(trackName, [tutorial], ''));
  }

  grok.events.onPackageLoaded.subscribe((p) => {
    if (p.name in tutorialIconPaths) {
      for (const tutorial in tutorialIconPaths[p.name]) {
        if (tutorialIconPaths[p.name].hasOwnProperty(tutorial)) {
          const data = tutorialIconPaths[p.name][tutorial];
          data.tutorial.imageUrl = `${data.tutorial.func.package.webRoot}${data.imageUrl}`;
        }
      }
    }
  });
}

//name: Demo
//tags: app
//description: Interactive demo of major Datagrok capabilities
export function demoApp() {
  let pathSegments = window.location.pathname.split('/');
  if (!pathSegments[pathSegments.length - 1])
    pathSegments.splice(pathSegments.length - 1, 1);

  const demoView = new DemoView();
  grok.shell.addView(demoView);

  if (pathSegments.length > PATH_START_INDEX) {
    const pathElements = pathSegments.slice(PATH_START_INDEX, pathSegments.length)
      .map((elem) => elem.replaceAll('%20', ' '));
    const path = pathElements.join('/');

    const func = DemoView.findDemoFunc(pathElements.join(' | '));
    func ? demoView.startDemoFunc(func, path) : demoView.nodeView(pathElements[pathElements.length - 1], path);
  }
}

function setProperties(properties: { [propertyName: string]: boolean }): void {
  if (tracks.length > 0)
    return;

  const registry: { [propertyName: string]: Track } = {
    'dataAnalysisTrack': eda,
    'machineLearningTrack': ml,
    'cheminformaticsTrack': chem,
    'dataAccessTrack': da,
    'dataTransformation': dataTransformation,
  };

  for (const property in properties) {
    if (property in registry && properties[property] === true)
      tracks.push(registry[property]);
  }
}

function setPath(path: string, tutorialRunners: TutorialRunner[]): void {
  const pathParts = path.split('/');
  const removeSpaces = (s: string) => s.replaceAll(' ', '');
  const trackShortNames: {[key: string]: Track} = {
    'eda': eda,
    'ml': ml,
    'chem': chem,
    'access': da,
  };

  if (pathParts.length !== 5)
    return;

  const [trackName, tutorialName] = pathParts.slice(3);
  let track: Track | null = null;
  let trackIdx: number | null = null;
  if (trackName in trackShortNames) {
    track = trackShortNames[trackName];
    const _idx = tracks.findIndex((t) => t === track);
    trackIdx = _idx === -1 ? null : _idx;
  } else {
    track = tracks.find((t, idx) => {
      if (removeSpaces(t.name) === trackName) {
        trackIdx = idx;
        return true;
      }
      return false;
    }) ?? null;
  }
  const tutorial = track?.tutorials?.find((t) => removeSpaces(t.name) === tutorialName);
  if (tutorial && trackIdx != null)
    tutorialRunners[trackIdx].run(tutorial);
}

//output: widget tutorial
export function demoAppWidget(): DG.Widget {
  return new DemoAppWidget();
}

//name: scatterPlotDemo
//description: A scatter plot (also called a scatter graph, scatter chart, scattergram, or scatter diagram) is a type of plot or mathematical diagram using Cartesian coordinates to display values for typically two variables for a set of data. If the points are color-coded you can increase the number of displayed variables to three. The data is displayed as a collection of points, each having the value of one variable determining the position on the horizontal axis and the value of the other variable determining the position on the vertical axis..
//meta.demoPath: Viewers | Scatter Plot
//test: _scatterPlotDemo()
export async function _scatterPlotDemo() {
  await viewerDemo(DG.VIEWER.SCATTER_PLOT);
}

//name: histogramDemo
//description: A histogram is a graphical representation of the distribution of numerical data.
//meta.demoPath: Viewers | General | Histogram
//test: _histogramDemo()
export async function _histogramDemo() {
  await viewerDemo(DG.VIEWER.HISTOGRAM);
}

//name: lineChartDemo
//description: Line chart is a simple, two-dimensional chart with an X and Y axis, each point representing a single value. The data points are joined by a line to depict a trend, usually over time.
//meta.demoPath: Viewers | General | Line Chart
//test: _lineChartDemo()
export async function _lineChartDemo() {
  await viewerDemo(DG.VIEWER.LINE_CHART);
}

//name: barChartDemo
//description: A bar chart presents grouped data as rectangular bars with lengths proportional to the values that they represent. Unlike histograms which you can apply to display the distribution of numerical data, bar charts are primarily designed for categorical values.
//meta.demoPath: Viewers | General | Bar Chart
//test: _barChartDemo()
export async function _barChartDemo() {
  await viewerDemo(DG.VIEWER.BAR_CHART);
}

//name: pieChartDemo
//description: Pie chart is useful for reflecting numerical proportions. Conceptually, it is similar to a bar chart in that it represents categorical values. A pie chart shows the relative size of a given category (a slice of the pie) compared to the entire dataset (the whole pie).
//meta.demoPath: Viewers | General | Pie Chart
//test: _pieChartDemo()
export async function _pieChartDemo() {
  await viewerDemo(DG.VIEWER.PIE_CHART);
}

//name: trellisPlotDemo
//description: Trellis Charts are useful for finding the structure and patterns in complex data. A Trellis Chart is a layout of smaller charts in a grid with consistent scales. Each smaller chart represents rows that belong to a corresponding category. The grid layout looks similar to a garden trellis, hence the name Trellis Chart.
//meta.demoPath: Viewers | Data separation | Trellis Plot
//test: _trellisPlotDemo()
export async function _trellisPlotDemo() {
  await viewerDemo(DG.VIEWER.TRELLIS_PLOT);
}

//name: matrixPlotDemo
//description: Use Matrix Plot to assess the relationship among many pairs of columns at the same time.
//meta.demoPath: Viewers | Data separation | Matrix Plot
//test: _matrixPlotDemo()
export async function _matrixPlotDemo() {
  await viewerDemo(DG.VIEWER.MATRIX_PLOT);
}

//name: scatterPlot3DDemo
//description: Use 3D scatter plot to plot data points on three axes to show the relationship between three variables. Each row in the data table is represented by a marker whose position depends on its values in the columns set on the X, Y, and Z axes. Additionally, you can color-code and size-code points, as well as display labels next to markers.
//meta.demoPath: Viewers | General | 3D Scatter Plot
//test: _scatterPlot3DDemo() //skip: GROK-13082
export async function _scatterPlot3DDemo() {
  await viewerDemo(DG.VIEWER.SCATTER_PLOT_3D);
}

//name: densityPlotDemo
//description: Unlike scatter plot that visualizes each individual data point, density plot splits 2D area by bins, and color-codes it depending on the number of points that fall within this bin. The darker the color, the more points it contains.
//meta.demoPath: Viewers | General | Density Plot
//test: _densityPlotDemo()
export async function _densityPlotDemo() {
  await viewerDemo(DG.VIEWER.DENSITY_PLOT);
}

//name: pcPlotDemo
//description: Parallel coordinates is a common way of visualizing high-dimensional geometry and analyzing multivariate data. To show a set of points in an n-dimensional space, a backdrop is drawn consisting of n parallel lines, typically vertical and equally spaced. A point in n-dimensional space is represented as a polyline with vertices on the parallel axes; the position of the vertex on the i-th axis corresponds to the i-th coordinate of the point. This visualization is closely related to time series visualization, except that it is applied to data where the axes do not correspond to points in time, and therefore do not have a natural order. Therefore, different axis arrangements may be of interest.
//meta.demoPath: Viewers | Statistical | PC Plot
//test: _pcPlotDemo()
export async function _pcPlotDemo() {
  await viewerDemo(DG.VIEWER.PC_PLOT);
}

//name: networkDiagramDemo
//description: Network diagram is used to visualize graphs, where values of the specified two columns become nodes, and rows become edges. It is possible to color-code and size-code nodes and columns by choosing the aggregate function that would apply to the values that represent an edge or a Node.js.
//meta.demoPath: Viewers | Data flow and hierarchy | Network Diagram
//test: _networkDiagramDemo() //skip: GROK-13082
export async function _networkDiagramDemo() {
  await viewerDemo(DG.VIEWER.NETWORK_DIAGRAM, {'node1ColumnName': 'Source', 'node2ColumnName': 'Target', useGoogleImage: true});
}

//name: boxPlotDemo
//description: The box plot (a.k.a. box and whisker diagram) is a standardized way of displaying the distribution of data based on the five number summary: minimum, first quartile, median, third quartile, and maximum.
//meta.demoPath: Viewers | Statistical | Box Plot
//test: _boxPlotDemo()
export async function _boxPlotDemo() {
  await viewerDemo(DG.VIEWER.BOX_PLOT);
}

//name: treeMapDemo
//description: Treemap displays hierarchical (tree-structured) data as nested rectangles. The branches are rectangles, then tiled with smaller rectangles representing sub-branches. Rectangles have areas proportional to a specified dimension of the data using the specified aggregation function (count by default).
//meta.demoPath: Viewers | Data flow and hierarchy | Tree map
//test: _treeMapDemo()
export async function _treeMapDemo() {
  await viewerDemo(DG.VIEWER.TREE_MAP);
}

//name: heatMapDemo
//description: Heatmap is a spreadsheet (grid) that contains colors instead of numbers and strings. For numerical data, the higher values are colored red, and the lower ones appear blue. The central value is assigned a light color so that darker colors indicate a larger distance from the center. For categorical data, each possible value is set to one color from a qualitative palette.
//meta.demoPath: Viewers | General | Heat map
//test: _heatMapDemo()
export async function _heatMapDemo() {
  await viewerDemo(DG.VIEWER.HEAT_MAP);
}

//name: statisticsDemo
//description: Provides specified descriptive statistics for the chosen columns.
//meta.demoPath: Viewers | Statistical | Statistics
//test: _statisticsDemo()
export async function _statisticsDemo() {
  await viewerDemo(DG.VIEWER.STATISTICS);
}

//name: correlationPlotDemo
//description: A quick way to assess correlations between all columns at once. Cells are color-coded by the Pearson correlation coefficient or Spearman's rank correlation coefficient . Histograms along the diagonal show the corresponding distribution. Hover over the cell to see the corresponding scatter plot. The grid is sortable. Select columns in the view by selecting corresponding rows.
//meta.demoPath: Viewers | Statistical |Correlation Plot
//test: _correlationPlotDemo()
export async function _correlationPlotDemo() {
  await viewerDemo(DG.VIEWER.CORR_PLOT);
}

//name: calendarDemo
//description: Calendar lets you analyze longitudinal data. It needs at least one column of type DateTime.
//meta.demoPath: Viewers | Time and date | Calendar
//test: calendarDemo()
export async function calendarDemo() {
  await viewerDemo(DG.VIEWER.CALENDAR);
}

//name: gridDemo
//description: A grid table contains a set of data that is structured in rows and columns. It allows the user to scroll in both directions and can handle large numbers of items and columns.
//meta.demoPath: Viewers | Input and edit | Grid
//test: _gridDemo()
export async function _gridDemo() {
  await viewerDemo(DG.VIEWER.GRID);
}

//name: markupDemo
//description: Use this viewer to host any text, arbitrary HTML content, or markdown-formatted text. In most casees, the viewer will auto-detect content type. Use the "mode" property to explicitly specify it.
//meta.demoPath: Viewers | General | Markup
//test: _markupDemo()
export async function _markupDemo() {
  await viewerDemo(DG.VIEWER.MARKUP);
}

//name: tileViewerDemo
//description: Visualizes rows as a collection of forms that are positioned as tiles.
//meta.demoPath: Viewers | General | Tile Viewer
//test: _tileViewerDemo()
export async function _tileViewerDemo() {
  await viewerDemo(DG.VIEWER.TILE_VIEWER);
}

//name: formDemo
//description: Form allows you to customize the appearance of the row by manually positioning the fields, and adding other visual elements, such as pictures or panels. A form can be used either as a stand-alone viewer, or as a row template of the Tile Viewer.
//meta.demoPath: Viewers | Input and edit | Form
//test: _formDemo()
export async function _formDemo() {
  await viewerDemo(DG.VIEWER.FORM);
}

//name: shapeMapDemo
//description: Shows a map that is applicable for the specified dataset. Typically, it would represent a geographical area (countries, states, counties, etc), but it also supports arbitrary shapes (such as a store floor plan, brain regions, or EEG electrodes).
//meta.demoPath: Viewers | Geographical | Shape Map
//test: _shapeMapDemo()
export async function _shapeMapDemo() {
  await viewerDemo(DG.VIEWER.SHAPE_MAP);
}

//name: pivotTableDemo
//description: A pivot table is a table of grouped values that aggregates the individual items of a more extensive table within one or more discrete categories. This summary might include sums, averages, or other statistics, which the pivot table groups together using a chosen aggregation function applied to the grouped values.
//meta.demoPath: Viewers | Pivot Table
//test: _pivotTableDemo()
export async function _pivotTableDemo() {
  await viewerDemo('Pivot table');
}

//name: mapDemo
//description: Map viewer shows geospatial data on a map as either markers, or a heat map.
//meta.demoPath: Viewers | Geographical | Map
//test: _mapDemo() //skip: GROK-13082
export async function _mapDemo() {
  await viewerDemo('Map');
}
