import {WordCloudViewer} from './viewers/word-cloud/word-cloud-viewer';
import {TreeViewer} from './viewers/tree/tree-viewer';
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

//name: Timelines
//description: Creates a timelines viewer
//output: viewer result
//meta.showInGallery: false
//meta.icon: icons/timelines-viewer.svg
//meta.role: viewer
export function timelinesViewer() : any {
  return PackageFunctions.timelinesViewer();
}

//output: grid_cell_renderer result
//meta.role: cellRenderer
export function flagCellRenderer() {
  return PackageFunctions.flagCellRenderer();
}

//description: Chord viewer visualizes weighted relationships between several entities
//meta.demoPath: Visualization | General | Chord
//meta.demoWait: 4000
export async function chordViewerDemo() : Promise<void> {
  await PackageFunctions.chordViewerDemo();
}

//description: Radar viewer is used on multivariate data to plot groups of values over several common variables
//meta.demoPath: Visualization | General | Radar
//meta.demoWait: 4000
export async function radarViewerDemo() : Promise<void> {
  await PackageFunctions.radarViewerDemo();
}

//description: Sankey viewer depicts a flow from one set of values to another
//meta.demoPath: Visualization | General | Sankey
//meta.demoWait: 4000
export async function sankeyViewerDemo() : Promise<void> {
  await PackageFunctions.sankeyViewerDemo();
}

//description: Sunburst viewer displays hierarchical data
//meta.demoPath: Visualization | General | Sunburst
//meta.demoWait: 4000
export async function sunburstViewerDemo() : Promise<void> {
  await PackageFunctions.sunburstViewerDemo();
}

//description: Surface plot viewer displays a set of three-dimensional data as a mesh surface
//meta.demoPath: Visualization | General | Surface Plot
//meta.demoWait: 4000
export async function surfacePlotDemo() : Promise<void> {
  await PackageFunctions.surfacePlotDemo();
}

//description: Timelines viewer displays the flow of events over time
//meta.demoPath: Visualization | General | Timelines
//meta.demoWait: 4000
export async function timelinesViewerDemo() : Promise<void> {
  await PackageFunctions.timelinesViewerDemo();
}

//description: Tree viewer visualizes hierarchical data by categories
//meta.demoPath: Visualization | Data Flow and Hierarchy | Tree
//meta.demoWait: 4000
export async function treeViewerDemo() : Promise<void> {
  await PackageFunctions.treeViewerDemo();
}

//description: Word Cloud viewer visualizes unstructured text data
//meta.demoPath: Visualization | General | Word Cloud
//meta.demoWait: 4000
export async function wordCloudViewerDemo() : Promise<void> {
  await PackageFunctions.wordCloudViewerDemo();
}
//name: Chord
//description: Creates a chord viewer
//output: viewer result
//meta.role: viewer
//meta.icon: icons/chord-viewer.svg
export function _ChordViewer() {
  return new ChordViewer();
}

//name: Globe
//description: Creates a globe viewer
//output: viewer result
//meta.role: viewer
//meta.icon: icons/globe-viewer.svg
export function _GlobeViewer() {
  return new GlobeViewer();
}

//name: Group Analysis
//description: Creates a group analysis viewer
//output: viewer result
//meta.role: viewer
//meta.icon: icons/groupanalysis-viewer.svg
export function _GroupAnalysisViewer() {
  return new GroupAnalysisViewer();
}

//name: Multiplot
//description: Creates a multiplot viewer
//output: viewer result
//meta.role: viewer
export function _MultiPlotViewer() {
  return new MultiPlotViewer();
}

//name: Radar
//description: Creates a radar viewer
//output: viewer result
//meta.role: viewer
//meta.icon: icons/radar-viewer.svg
export function _RadarViewer() {
  return new RadarViewer();
}

//name: Sankey
//description: Creates a sankey viewer
//output: viewer result
//meta.role: viewer
//meta.icon: icons/sankey-viewer.svg
export function _SankeyViewer() {
  return new SankeyViewer();
}

//name: Sunburst
//description: Creates a sunburst viewer
//output: viewer result
//meta.role: viewer
//meta.icon: icons/sunburst-viewer.svg
export function _SunburstViewer() {
  return new SunburstViewer();
}

//name: Surface plot
//description: Creates a surface plot viewer
//output: viewer result
//meta.role: viewer
//meta.icon: icons/surfaceplot-viewer.svg
export function _SurfacePlot() {
  return new SurfacePlot();
}

//name: Tree
//description: Creates a tree viewer
//output: viewer result
//meta.role: viewer
//meta.trellisable: false
//meta.icon: icons/tree-viewer.svg
export function _TreeViewer() {
  return new TreeViewer();
}

//name: Word cloud
//description: Creates a word cloud viewer
//output: viewer result
//meta.role: viewer
//meta.icon: icons/wordcloud-viewer.svg
//meta.toolbox: true
export function _WordCloudViewer() {
  return new WordCloudViewer();
}

