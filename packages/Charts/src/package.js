import {TreeViewer} from './tree-viewer.js';
import {TreeMapViewer} from './tree-map-viewer.js';
import {SunburstViewer} from './sunburst-viewer.js';
import {RadarViewer} from './radarviewer.js';
import {TimelinesViewer} from './timelinesviewer.js';
import {SankeyViewer} from './sankey-viewer.js';
import {ChordViewer} from './chord-viewer.js';

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

//name: TimelinesViewer
//description: Creates TimelinesViewer viewer
//tags: viewer
//output: viewer result
export function _TimelinesViewer() {
    return new TimelinesViewer();
}

//name: timelinesViewerDemo
export function timelinesViewerDemo() {
    let adverseEvents = DG.DataFrame.fromCsv(
        `subject, start, end
        s1, 10/02/2018, 10/03/2018
        s2, 10/03/2018, 10/04/2018
        s3, 10/02/2018, 10/05/2018`);

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
