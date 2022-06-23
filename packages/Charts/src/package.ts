import {TreeViewer} from './tree-viewer';
import {TreeMapViewer} from './tree-map-viewer';
import {SunburstViewer} from './sunburst-viewer';
import {RadarViewer} from './radar-viewer';
import {TimelinesViewer} from './timelines/timelines-viewer';
import {SankeyViewer} from './sankey-viewer';
import {ChordViewer} from './chord-viewer';
import { WordCloudViewer } from './word-cloud-viewer';
import './timelines/styles.css';

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

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
//meta.trellisable: true
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

//name: TreeMapViewer
//tags: viewer
//meta.trellisable: true
//output: viewer result
export function _TreeMapViewer() {
  return new TreeMapViewer();
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
//meta.trellisable: true
//output: viewer result
export function _ChordViewer() {
  return new ChordViewer();
}

//name: WordCloudViewer
//tags: viewer
//meta.trellisable: true
//output: viewer result
export function _WordCloudViewer() {
  return new WordCloudViewer();
}

//tags: init
export function init() {
  grok.shell.registerViewer('TimelinesViewer', 'Creates TimelinesViewer viewer', () => new TimelinesViewer());
  grok.events.onContextMenu.subscribe((args) => {
    if (args.args.context instanceof TimelinesViewer) {
      args.args.menu.item('Reset View', () => {
        args.args.context.zoomState = [[0, 100], [0, 100], [0, 100], [0, 100]];
        args.args.context.render();
      });
    }
  });
}
