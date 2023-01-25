import * as grok from 'datagrok-api/grok';

export function radarViewerDemo() {
  const tv = grok.shell.addTableView(grok.data.demo.demog());
  tv.addViewer('RadarViewer');
}
