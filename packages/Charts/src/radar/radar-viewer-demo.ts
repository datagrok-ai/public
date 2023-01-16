import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export function radarViewerDemo() {
  let tv = grok.shell.addTableView(grok.data.demo.demog());
  tv.addViewer('RadarViewer');
}