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
//meta.demoPath: Viewers | Chord
export function _chordViewerDemo() {
  viewerDemo('Chord');
}

//name: globeViewerDemo
//meta.demoPath: Viewers | Globe
export function _globeViewerDemo() {
  viewerDemo('Globe');
}

//name: groupAnalysisViewerDemo
//meta.demoPath: Viewers | GroupAnalysis
export function _groupAnalysisViewerDemo() {
  viewerDemo('GroupAnalysis');
}

//name: radarViewerDemo
//meta.demoPath: Viewers | Radar
export function _radarViewerDemo() {
  viewerDemo('Radar');
}

//name: sankeyViewerDemo
//meta.demoPath: Viewers | Sankey
export function _sankeyViewerDemo() {
  viewerDemo('Sankey');
}

//name: sunburstViewerDemo
//meta.demoPath: Viewers | Sunburst
export function _sunburstViewerDemo() {
  viewerDemo('Sunburst');
}

//name: surfacePlotDemo
//meta.demoPath: Viewers | SurfacePlot
export function _surfacePlotDemo() {
  viewerDemo('SurfacePlot');
}

//name: timelinesViewerDemo
//meta.demoPath: Viewers | Timelines
export function _timelinesViewerDemo() {
  viewerDemo('Timelines', {lineWidth: 4, markerPosition: 'above main line'});
}

//name: treeViewerDemo
//meta.demoPath: Viewers | Tree
export function _treeViewerDemo() {
  viewerDemo('Tree');
}

//name: wordCloudViewerDemo
//meta.demoPath: Viewers | WordCloud
export function _wordCloudViewerDemo() {
  viewerDemo('WordCloud');
}
