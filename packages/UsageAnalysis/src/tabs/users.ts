/*
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaQueryViewer} from '../viewers/abstract/ua-query-viewer';
import {TopUsersViewer} from '../drilldown_viewers/users/top-users-viewer';

export class UsersView extends UaView {
  static viewName = 'Users';

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
  }

  async initViewers(): Promise<void> {
    const uniqueUsersViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Unique Users',
      'UniqueUsers',
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.lineChart(t, UaQueryViewer.defaultChartOptions);
        viewer.root.style.maxHeight = '150px';
        return viewer;
      },
    );
    this.viewers.push(uniqueUsersViewer);

    const topUsersViewer = new TopUsersViewer(this.uaToolbox.filterStream);
    this.viewers.push(topUsersViewer);

    const usageViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Usage',
      'Usage',
      (t: DG.DataFrame) => DG.Viewer.scatterPlot(t, UaQueryViewer.defaultChartOptions),
    );
    this.viewers.push(usageViewer);

    const topPackagesByUsers = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Packages By Users',
      'TopPackagesByUsers',
      (t: DG.DataFrame) => DG.Viewer.scatterPlot(t, UaQueryViewer.defaultChartOptions),
    );
    this.viewers.push(topPackagesByUsers);

    this.root.append(ui.divV([
      ui.divH([ui.block([uniqueUsersViewer.root])]),
      ui.divH([ui.block([usageViewer.root])]),
      ui.divH([ui.block([topUsersViewer.root])]),
    ]));
  }
}
*/
