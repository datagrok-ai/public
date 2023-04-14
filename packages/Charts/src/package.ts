import * as DG from 'datagrok-api/dg';

import {FlagCellRenderer} from './renderers/flag-cell-renderer';

import {TreeViewer} from './viewers/tree/tree-viewer';
import {SunburstViewer} from './viewers/sunburst/sunburst-viewer';
import {RadarViewer} from './viewers/radar/radar-viewer';
import {TimelinesViewer} from './viewers/timelines/timelines-viewer';
import {SankeyViewer} from './viewers/sankey/sankey';
import {ChordViewer} from './viewers/chord/chord-viewer';
import {WordCloudViewer} from './viewers/word-cloud/word-cloud-viewer';
import {GroupAnalysisViewer} from './viewers/group-analysis/group-analysis-viewer';
import {SurfacePlot} from './viewers/surface-plot/surface-plot';
import {GlobeViewer} from './viewers/globe/globe-viewer';

import {viewerDemo} from './demos/demo';


export const _package = new DG.Package();


//name: flagCellRenderer
//tags: cellRenderer
//meta-cell-renderer-sem-type: flag
//output: grid_cell_renderer result
export function flagCellRenderer() {
  return new FlagCellRenderer();
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
//meta.icon: icons/surfaceplot-viewer.svg
//output: viewer result
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
//meta.trellisable: true
//output: viewer result
//meta.icon: icons/tree-viewer.svg
export function _TreeViewer() {
  return new TreeViewer();
}

//name: Word Cloud
//description: Creates a word cloud viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/wordcloud-viewer.svg
export function _WordCloudViewer() {
  return new WordCloudViewer();
}


//name: chordViewerDemo
//description: Chord viewer visualizes weighted relationships between several entities\nPart of the Charts package
//meta.demoPath: Viewers | Chord
export async function _chordViewerDemo() {
  await viewerDemo('Chord');
}

//name: globeViewerDemo
//description: Globe viewer represents data visualization layers on a 3-dimensional globe in a spherical projection\nPart of the Charts package
//meta.demoPath: Viewers | Globe
export async function _globeViewerDemo() {
  await viewerDemo('Globe');
}

//name: groupAnalysisViewerDemo
//description: Group analysis viewer groups data by different options\nPart of the Charts package
//meta.demoPath: Viewers | Group Analysis
export async function _groupAnalysisViewerDemo() {
  await viewerDemo('GroupAnalysis');
}

//name: radarViewerDemo
//description: Radar viewer is used on multivariate data to plot groups of values over several common variables\nPart of the Charts package
//meta.demoPath: Viewers | Radar
export async function _radarViewerDemo() {
  await viewerDemo('Radar');
}

//name: sankeyViewerDemo
//description: Sankey viewer depicts a flow from one set of values to another\nPart of the Charts package
//meta.demoPath: Viewers | Sankey
export async function _sankeyViewerDemo() {
  await viewerDemo('Sankey');
}

//name: sunburstViewerDemo
//description: Sunburst viewer displays hierarchical data\nPart of the Charts package
//meta.demoPath: Viewers | Sunburst
export async function _sunburstViewerDemo() {
  await viewerDemo('Sunburst');
}

//name: surfacePlotDemo
//description: Surface plot viewer displays a set of three-dimensional data as a mesh surface\nPart of the Charts package
//meta.demoPath: Viewers | Surface Plot
export async function _surfacePlotDemo() {
  await viewerDemo('SurfacePlot');
}

//name: timelinesViewerDemo
//description: Timelines viewer displays the flow of events over time\nPart of the Charts package
//meta.demoPath: Viewers | Timelines
export async function _timelinesViewerDemo() {
  await viewerDemo('Timelines', {lineWidth: 4, markerPosition: 'above main line'});
}

//name: treeViewerDemo
//description: Tree viewer visualizes hierarchical data by categories\nPart of the Charts package
//meta.demoPath: Viewers | Tree
export async function _treeViewerDemo() {
  await viewerDemo('Tree', {left: '40px', right: '75px'});
}

//name: wordCloudViewerDemo
//description: Word Cloud viewer visualizes unstructured text data\nPart of the Charts package
//meta.demoPath: Viewers | Word Cloud
export async function _wordCloudViewerDemo() {
  await viewerDemo('WordCloud');
}
