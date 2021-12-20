import {TreeViewer} from './tree-viewer.js';
import {TreeMapViewer} from './tree-map-viewer.js';
import {SunburstViewer} from './sunburst-viewer.js';
import {RadarViewer} from './radar-viewer.js';
import {TimelinesViewer} from './timelines/timelines-viewer.js';
import {SankeyViewer} from './sankey-viewer.js';
import {ChordViewer} from './chord-viewer.js';
import { WordCloudViewer } from './word-cloud-viewer.js';

/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//name: test
//input: string s
export function test(s) {
  grok.shell.info(_package.webRoot);
}


//name: timelinesViewerDemo
export function timelinesViewerDemo() {
  let adverseEvents = DG.DataFrame.fromCsv(
    `USUBJID, AESTDY, AEENDY, SEX, AGE
    s1, 10/02/2018, 10/09/2018, F, 48
    s2, 10/04/2018, 10/07/2018, M, 51
    s3, 10/02/2018, 10/05/2018, F, 39
    s4, 10/07/2018, 10/08/2018, M, 43`);

  let view = grok.shell.addTableView(adverseEvents);
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
  let view = grok.shell.addTableView(grok.data.demo.demog(100));
  view.addViewer('RadarViewer');
}

//name: TreeViewer
//tags: viewer
//output: viewer result
export function _TreeViewer() {
  return new TreeViewer();
}

//name: TreeMapViewer
//tags: viewer
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

//tags: autostart
export function init() {
  grok.shell.registerViewer('TimelinesViewer', 'Creates TimelinesViewer viewer', () => new TimelinesViewer());
  grok.events.onContextMenu.subscribe(args => {
    if (args.args.context instanceof TimelinesViewer) {
      args.args.menu.item('Reset View', () => {
        args.args.context.zoomState = [[0, 100], [0, 100], [0, 100], [0, 100]];
        args.args.context.render();
      });
    }
  });
}
