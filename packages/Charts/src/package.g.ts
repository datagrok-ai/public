import {WordCloudViewer} from './viewers/word-cloud/word-cloud-viewer';
import {TreeViewer} from './viewers/tree/tree-viewer';
import {TimelinesViewer} from './viewers/timelines/timelines-viewer';
import {SurfacePlot} from './viewers/surface-plot/surface-plot';
import {SunburstViewer} from './viewers/sunburst/sunburst-viewer';
import {SankeyViewer} from './viewers/sankey/sankey';
import {RadarViewer} from './viewers/radar/radar-viewer';
import {MultiPlotViewer} from './viewers/multiplot/multiplot';
import {GroupAnalysisViewer} from './viewers/group-analysis/group-analysis-viewer';
import {GlobeViewer} from './viewers/globe/globe-viewer';
import {ChordViewer} from './viewers/chord/chord-viewer';
import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: flagCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
export function flagCellRenderer() {
  return PackageFunctions.flagCellRenderer();
}

//name: chordViewerDemo
//description: Chord viewer visualizes weighted relationships between several entities
//meta.demoPath: Visualization | General | Chord
//meta.demoWait: 4000
export async function chordViewerDemo() {
  return PackageFunctions.chordViewerDemo();
}

//name: radarViewerDemo
//description: Radar viewer is used on multivariate data to plot groups of values over several common variables
//meta.demoPath: Visualization | General | Radar
//meta.demoWait: 4000
export async function radarViewerDemo() {
  return PackageFunctions.radarViewerDemo();
}

//name: sankeyViewerDemo
//description: Sankey viewer depicts a flow from one set of values to another
//meta.demoPath: Visualization | General | Sankey
//meta.demoWait: 4000
export async function sankeyViewerDemo() {
  return PackageFunctions.sankeyViewerDemo();
}

//name: sunburstViewerDemo
//description: Sunburst viewer displays hierarchical data
//meta.demoPath: Visualization | General | Sunburst
//meta.demoWait: 4000
export async function sunburstViewerDemo() {
  return PackageFunctions.sunburstViewerDemo();
}

//name: surfacePlotDemo
//description: Surface plot viewer displays a set of three-dimensional data as a mesh surface
//meta.demoPath: Visualization | General | Surface Plot
//meta.demoWait: 4000
export async function surfacePlotDemo() {
  return PackageFunctions.surfacePlotDemo();
}

//name: timelinesViewerDemo
//description: Timelines viewer displays the flow of events over time
//meta.demoPath: Visualization | General | Timelines
//meta.demoWait: 4000
export async function timelinesViewerDemo() {
  return PackageFunctions.timelinesViewerDemo();
}

//name: treeViewerDemo
//description: Tree viewer visualizes hierarchical data by categories
//meta.demoPath: Visualization | Data Flow and Hierarchy | Tree
//meta.demoWait: 4000
export async function treeViewerDemo() {
  return PackageFunctions.treeViewerDemo();
}

//name: wordCloudViewerDemo
//description: Word Cloud viewer visualizes unstructured text data
//meta.demoPath: Visualization | General | Word Cloud
//meta.demoWait: 4000
export async function wordCloudViewerDemo() {
  return PackageFunctions.wordCloudViewerDemo();
}
//name: Chord
//description: Creates a chord viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/chord-viewer.svg
export function _ChordViewer() {
  return new ChordViewer();
}

//name: Globe
//description: Creates a globe viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/globe-viewer.svg
export function _GlobeViewer() {
  return new GlobeViewer();
}

//name: Group Analysis
//description: Creates a group analysis viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/groupanalysis-viewer.svg
export function _GroupAnalysisViewer() {
  return new GroupAnalysisViewer();
}

//name: Multiplot
//description: Creates a multiplot viewer
//tags: viewer
//output: viewer result
export function _MultiPlotViewer() {
  return new MultiPlotViewer();
}

//name: Radar
//description: Creates a radar viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/radar-viewer.svg
export function _RadarViewer() {
  return new RadarViewer();
}

//name: Sankey
//description: Creates a sankey viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/sankey-viewer.svg
export function _SankeyViewer() {
  return new SankeyViewer();
}

//name: Sunburst
//description: Creates a sunburst viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/sunburst-viewer.svg
export function _SunburstViewer() {
  return new SunburstViewer();
}

//name: Surface plot
//description: Creates a surface plot viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/surfaceplot-viewer.svg
export function _SurfacePlot() {
  return new SurfacePlot();
}

//name: Timelines
//description: Creates a timelines viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/timelines-viewer.svg
export function _TimelinesViewer() {
  return new TimelinesViewer();
}

//name: Tree
//description: Creates a tree viewer
//tags: viewer
//output: viewer result
//meta.trellisable: false
//meta.icon: icons/tree-viewer.svg
export function _TreeViewer() {
  return new TreeViewer();
}

//name: Word cloud
//description: Creates a word cloud viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/wordcloud-viewer.svg
//meta.toolbox: true
export function _WordCloudViewer() {
  return new WordCloudViewer();
}

