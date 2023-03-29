import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';


export class FunctionsView extends UaView {
  // static viewName = 'Functions';

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox, 'function', true);
    this.name = 'Functions';
  }

  async initViewers(): Promise<void> {
    const functionsViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Functions',
      'FunctionsUsage',
      (t: DG.DataFrame) => {
        t.onCurrentRowChanged.subscribe(async () => {
          const rowValues = Array.from(t.currentRow.cells).map((c) => c.value);
          const row = Object.fromEntries(t.columns.names().map((k, i) => [k, rowValues[i]]));
          row.time_start = row.time_start.a;
          row.time_end = row.time_end.a;
          const cp = DG.Accordion.create();
          cp.addPane('Details', () => {
            return ui.tableFromMap({'User': ui.render(`#{x.${row.uid}}`),
              'Package': ui.render(`#{x.${row.pid}}`),
              'From': new Date(row.time_start).toLocaleString(),
              'To': new Date(row.time_end).toLocaleString()});
          }, true);
          grok.shell.o = cp.root;
        });
      }, null, null, this.viewer);
    this.viewers.push(functionsViewer);
    this.root.append(functionsViewer.root);
  }
}
