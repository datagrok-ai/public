import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

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

import {FlagCellRenderer} from './renderers/flag-cell-renderer';

import '../css/styles.css';
import '../css/chord-viewer.css';
import '../css/sankey.css';

export const _package = new DG.Package();


//name: timelinesViewerDemo
export function timelinesViewerDemo() {
  const adverseEvents = DG.DataFrame.fromCsv(
    `USUBJID, AESTDY, AEENDY, SEX, AGE
    s1, 10/02/2018, 10/09/2018, F, 48
    s2, 10/04/2018, 10/07/2018, M, 51
    s3, 10/02/2018, 10/05/2018, F, 39
    s4, 10/07/2018, 10/08/2018, M, 43`);

  const view = grok.shell.addTableView(adverseEvents);
  view.addViewer('TimelinesViewer');
}

//name: RadarViewer
//tags: viewer
//output: viewer result
export function _RadarViewer() {
  return new RadarViewer();
}

//name: radarViewerDemo
export function radarViewerDemo() {
  const view = grok.shell.addTableView(grok.data.demo.demog(100));
  view.addViewer('RadarViewer');
}

//name: TreeViewer
//tags: viewer
//meta.trellisable: true
//output: viewer result
export function _TreeViewer() {
  return new TreeViewer();
}

//name: SunburstViewer
//tags: viewer
//output: viewer result
export function _SunburstViewer() {
  return new SunburstViewer();
}

//name: SankeyViewer
//tags: viewer
//output: viewer result
export function _SankeyViewer() {
  return new SankeyViewer();
}

//name: ChordViewer
//tags: viewer
//output: viewer result
export function _ChordViewer() {
  return new ChordViewer();
}

//name: WordCloudViewer
//tags: viewer
//output: viewer result
export function _WordCloudViewer() {
  return new WordCloudViewer();
}

//name: TimelinesViewer
//tags: viewer
//output: viewer result
export function _TimelinesViewer() {
  return new TimelinesViewer();
}

//name: SurfacePlot
//tags: viewer
//output: viewer result
export function _SurfacePlot() {
  return new SurfacePlot();
}

//name: GroupAnalysisViewer
//tags: viewer
//output: viewer result
export function _GroupAnalysisViewer() {
  return new GroupAnalysisViewer();
}

//name: radarViewerDemo
//meta.demoPath: Viewers | Radar
export function _radarViewerDemo() {
  radarViewerDemo();
}

//name: radarViewerDemo
//meta.demoPath: Viewers | Timelines
export function _timelinesViewerDemo() {
  const tv = grok.shell.addTableView(grok.data.demo.demog());
  tv.addViewer(DG.VIEWER.TIMELINES);
}

//name: GlobeViewer
//description: Creates a globe viewer
//tags: viewer
//output: viewer result
export function _GlobeViewer() {
  return new GlobeViewer();
}

//name: flagCellRenderer
//tags: cellRenderer
//meta-cell-renderer-sem-type: flag
//output: grid_cell_renderer result
export function flagCellRenderer() {
  return new FlagCellRenderer();
}
