import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaDataFrameQueryViewer} from '../viewers/ua-data-frame-query-viewer';
// import {UaQueryViewer} from '../viewers/abstract/ua-query-viewer';
import {PropertyPanel} from '../property-panel';


export class PackagesView extends UaView {
  static viewName = 'Packages';

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
  }

  async initViewers(): Promise<void> {
    const packagesViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Packages',
      'PackagesUsage',
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.scatterPlot(t, {
          x: 'time',
          y: 'package',
          size: 'count',
          color: 'user',
        });
        t.onCurrentRowChanged.subscribe((args) => {
          const params = Array.from(t.currentRow.cells).map((c) => c.value);
          params[3] = params[3].a / 1000;
          const pp = new PropertyPanel(
            null,
            null,
            [new UaDataFrameQueryViewer(
              'Events table',
              'PackageUsageAtPoint',
              (t: DG.DataFrame) => DG.Viewer.grid(t).root,
              undefined,
              {date: params[3], user: params[4], package: params[0]},
              this.uaToolbox.filterStream.getValue(),
              false,
            ),
            ], `Info: ${params[1]}, ${params[0]}`, 'Info');
          pp.getRoot().appendChild(ui.button('Functions', () => {
            this.uaToolbox.usersInput.choices.removeActiveItems(-1);
            const user = params[4].charAt(0).toUpperCase() + params[4].slice(1);
            this.uaToolbox.usersInput.choices._addItem({value: user, label: user});
            this.uaToolbox.packagesInput.choices.removeActiveItems(-1);
            this.uaToolbox.packagesInput.choices._addItem({value: params[0], label: params[0]});
            this.uaToolbox.applyFilter();
          }));
          grok.shell.o = pp.getRoot();
        });
        return viewer.root;
      },
    );
    this.viewers.push(packagesViewer);
    this.root.append(packagesViewer.root);
  }
}
