import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {TutorialRunner, TutorialSubstitute} from './tutorial-runner';
import {chem} from './tracks/chem';
import {eda} from './tracks/eda';
import {da} from './tracks/data-access';
import {ml} from './tracks/ml';
import {TutorialWidget} from './widget';
import '../css/tutorial.css';
import {Track} from '@datagrok-libraries/tutorials/src/track';
import {DemoView} from './demo-app/demo-app';
import {viewerDemo} from './demo-app/platform-viewers-demo';


export const _package = new DG.Package();

const tracks: Track[] = [];

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
  const pathSegments = window.location.pathname.split('/');
  const demoView = new DemoView();
  grok.shell.addView(demoView);

  if (pathSegments.length > 4) {
    const category = pathSegments[4];
    if (pathSegments[pathSegments.length - 1] === category) {
      demoView.nodeView(category);
      return;
    }

    const viewerName = pathSegments[5].replaceAll('%20', ' ');
    const f = DemoView.findDemoFunc(`${category} | ${viewerName}`);
    if (f) {
      const viewPath = `${category}/${viewerName}`;
      demoView.startDemoFunc(f, viewPath);
    }
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


//name: scatterPlotDemo
//description: Scatter Plot ....
//meta.demoPath: Viewers | Scatter Plot
export async function _scatterPlotDemo() {
  await viewerDemo(DG.VIEWER.SCATTER_PLOT);
}

//name: histogramDemo
//description: Histogram ....
//meta.demoPath: Viewers | Histogram
export async function _histogramDemo() {
  await viewerDemo(DG.VIEWER.HISTOGRAM);
}

//name: lineChartDemo
//description: Line Chart ....
//meta.demoPath: Viewers | Line Chart
export async function _lineChartDemo() {
  await viewerDemo(DG.VIEWER.LINE_CHART);
}

//name: barChartDemo
//description: Bar Chart ....
//meta.demoPath: Viewers | Bar Chart
export async function _barChartDemo() {
  await viewerDemo(DG.VIEWER.BAR_CHART);
}

//name: pieChartDemo
//description: Pie Chart ....
//meta.demoPath: Viewers | Pie Chart
export async function _pieChartDemo() {
  await viewerDemo(DG.VIEWER.PIE_CHART);
}

//name: trellisPlotDemo
//description: Trellis Plot ....
//meta.demoPath: Viewers | Trellis Plot
export async function _trellisPlotDemo() {
  await viewerDemo(DG.VIEWER.TRELLIS_PLOT);
}

//name: matrixPlotDemo
//description: Matrix Plot ....
//meta.demoPath: Viewers | Matrix Plot
export async function _matrixPlotDemo() {
  await viewerDemo(DG.VIEWER.MATRIX_PLOT);
}

//name: scatterPlot3DDemo
//description: 3D Scatter Plot ....
//meta.demoPath: Viewers | 3D Scatter Plot
export async function _scatterPlot3DDemo() {
  await viewerDemo(DG.VIEWER.SCATTER_PLOT_3D);
}

//name: densityPlotDemo
//description: Density Plot ....
//meta.demoPath: Viewers | Density Plot
export async function _densityPlotDemo() {
  await viewerDemo(DG.VIEWER.DENSITY_PLOT);
}

//name: pcPlotDemo
//description: PC Plot ....
//meta.demoPath: Viewers | PC Plot
export async function _pcPlotDemo() {
  await viewerDemo(DG.VIEWER.PC_PLOT);
}

//name: networkDiagramDemo
//description: Network Diagram ....
//meta.demoPath: Viewers | Network Diagram
export async function _networkDiagramDemo() {
  await viewerDemo(DG.VIEWER.NETWORK_DIAGRAM);
}

//name: boxPlotDemo
//description: Box Plot ....
//meta.demoPath: Viewers | Box Plot
export async function _boxPlotDemo() {
  await viewerDemo(DG.VIEWER.BOX_PLOT);
}

//name: treeMapDemo
//description: Tree map ....
//meta.demoPath: Viewers | Tree map
export async function _treeMapDemo() {
  await viewerDemo(DG.VIEWER.TREE_MAP);
}

//name: heatMapDemo
//description: Heat map ....
//meta.demoPath: Viewers | Heat map
export async function _heatMapDemo() {
  await viewerDemo(DG.VIEWER.HEAT_MAP);
}

//name: statisticsDemo
//description: Statistics ....
//meta.demoPath: Viewers | Statistics
export async function _statisticsDemo() {
  await viewerDemo(DG.VIEWER.STATISTICS);
}

//name: correlationPlotDemo
//description: Correlation Plot ....
//meta.demoPath: Viewers | Correlation Plot
export async function _correlationPlotDemo() {
  await viewerDemo(DG.VIEWER.CORR_PLOT);
}

//name: calendarDemo
//description: Calendar ....
//meta.demoPath: Viewers | Calendar
export async function calendarDemo() {
  await viewerDemo(DG.VIEWER.CALENDAR);
}

//name: gridDemo
//description: Grid ....
//meta.demoPath: Viewers | Grid
export async function _gridDemo() {
  await viewerDemo(DG.VIEWER.GRID);
}

//name: markupDemo
//description: Markup ....
//meta.demoPath: Viewers | Markup
export async function _markupDemo() {
  await viewerDemo(DG.VIEWER.MARKUP);
}

//name: tileViewerDemo
//description: Tile viewer ....
//meta.demoPath: Viewers | Tile Viewer
export async function _tileViewerDemo() {
  await viewerDemo(DG.VIEWER.TILE_VIEWER);
}

//name: formDemo
//description: Form ....
//meta.demoPath: Viewers | Form
export async function _formDemo() {
  await viewerDemo(DG.VIEWER.FORM);
}

//name: shapeMapDemo
//description: Shape Map ....
//meta.demoPath: Viewers | Shape Map
export async function _shapeMapDemo() {
  await viewerDemo(DG.VIEWER.SHAPE_MAP);
}

//name: pivotTableDemo
//description: Pivot Table ....
//meta.demoPath: Viewers | Pivot Table
export async function _pivotTableDemo() {
  await viewerDemo('Pivot table');
}

//name: mapDemo
//description: Map ....
//meta.demoPath: Viewers | Map
export async function _mapDemo() {
  await viewerDemo('Map');
}
