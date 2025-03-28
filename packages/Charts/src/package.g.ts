import {ChordViewer} from './viewers/chord/chord-viewer';
import {GlobeViewer} from './viewers/globe/globe-viewer';
import {GroupAnalysisViewer} from './viewers/group-analysis/group-analysis-viewer';
import {MultiPlotViewer} from './viewers/multiplot/multiplot';
import {RadarViewer} from './viewers/radar/radar-viewer';
import {SankeyViewer} from './viewers/sankey/sankey';
import {SunburstViewer} from './viewers/sunburst/sunburst-viewer';
import {SurfacePlot} from './viewers/surface-plot/surface-plot';
import {TimelinesViewer} from './viewers/timelines/timelines-viewer';
import {TreeViewer} from './viewers/tree/tree-viewer';
import {WordCloudViewer} from './viewers/word-cloud/word-cloud-viewer';

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

