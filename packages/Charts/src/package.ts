import {_MultiPlotViewer} from './package.g';
import {_WordCloudViewer} from './package.g';
import {_TreeViewer} from './package.g';
import {_TimelinesViewer} from './package.g';
import {_SurfacePlot} from './package.g';
import {_SunburstViewer} from './package.g';
import {_SankeyViewer} from './package.g';
import {_RadarViewer} from './package.g';
import {_GroupAnalysisViewer} from './package.g';
import {_GlobeViewer} from './package.g';
import {_ChordViewer} from './package.g';
import * as DG from 'datagrok-api/dg';
import {FlagCellRenderer} from './renderers/flag-cell-renderer';
import {viewerDemo} from './demos/demo';


export const _package = new DG.Package();


//name: flagCellRenderer
//tags: cellRenderer
//meta-cell-renderer-sem-type: flag
//output: grid_cell_renderer result
export function flagCellRenderer() {
  return new FlagCellRenderer();
}

//name: chordViewerDemo
//description: Chord viewer visualizes weighted relationships between several entities
//meta.demoPath: Visualization | General | Chord
//meta.demoWait: 4000
export async function _chordViewerDemo() {
  await viewerDemo('Chord');
}

// //name: globeViewerDemo
// //description: Globe viewer represents data visualization layers on a 3-dimensional globe in a spherical projection
// //meta.demoPath: Visualization | Geographical | Globe
// //meta.demoWait: 4000
// //meta.demoSkip: GROK-14320
// export async function _globeViewerDemo() {
//   await viewerDemo('Globe');
// }

//name: radarViewerDemo
//description: Radar viewer is used on multivariate data to plot groups of values over several common variables
//meta.demoPath: Visualization | General | Radar
//meta.demoWait: 4000
export async function _radarViewerDemo() {
  await viewerDemo('Radar');
}

//name: sankeyViewerDemo
//description: Sankey viewer depicts a flow from one set of values to another
//meta.demoPath: Visualization | General | Sankey
//meta.demoWait: 4000
export async function _sankeyViewerDemo() {
  await viewerDemo('Sankey');
}

//name: sunburstViewerDemo
//description: Sunburst viewer displays hierarchical data
//meta.demoPath: Visualization | General | Sunburst
//meta.demoWait: 4000
export async function _sunburstViewerDemo() {
  await viewerDemo('Sunburst');
}

//name: surfacePlotDemo
//description: Surface plot viewer displays a set of three-dimensional data as a mesh surface
//meta.demoPath: Visualization | General | Surface Plot
//meta.demoWait: 4000
export async function _surfacePlotDemo() {
  await viewerDemo('SurfacePlot');
}

//name: timelinesViewerDemo
//description: Timelines viewer displays the flow of events over time
//meta.demoPath: Visualization | General | Timelines
//meta.demoWait: 4000
export async function _timelinesViewerDemo() {
  await viewerDemo('Timelines', {lineWidth: 4, markerPosition: 'above main line'});
}

//name: treeViewerDemo
//description: Tree viewer visualizes hierarchical data by categories
//meta.demoPath: Visualization | Data Flow and Hierarchy | Tree
//meta.demoWait: 4000
export async function _treeViewerDemo() {
  await viewerDemo('Tree', {left: '40px', right: '75px'});
}

//name: wordCloudViewerDemo
//description: Word Cloud viewer visualizes unstructured text data
//meta.demoPath: Visualization | General | Word Cloud
//meta.demoWait: 4000
export async function _wordCloudViewerDemo() {
  await viewerDemo(DG.VIEWER.WORD_CLOUD, {drawOutOfBound: false});
}

export {_ChordViewer};
export {_GlobeViewer};
export {_GroupAnalysisViewer};
export {_RadarViewer};
export {_SankeyViewer};
export {_SunburstViewer};
export {_SurfacePlot};
export {_TimelinesViewer};
export {_TreeViewer};
export {_WordCloudViewer};
export {_MultiPlotViewer};
