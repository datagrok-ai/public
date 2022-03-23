import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../../css/usage_analysis.css';
import {UaToolbox} from "../ua-toolbox";
import {UaView} from "./ua-view";
import {UaFilter} from "../filter2";
import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";
import {UaQueryViewer} from "../viewers/abstract/ua-query-viewer";
import {TopUsersViewer} from "../drilldown_viewers/users/top-users-viewer";

export class UsersView extends UaView {

  constructor(uaToolbox: UaToolbox) {
    super('Users', uaToolbox);
  }

  async initViewers(): Promise<void> {
    let uniqueUsersViewer = new UaFilterableQueryViewer(
        this.uaToolbox.filterStream,
        'Unique Users',
        'UniqueUsers',
        (t: DG.DataFrame) => DG.Viewer.lineChart(t, UaQueryViewer.defaultChartOptions).root
    );
    this.viewers.push(uniqueUsersViewer);

    let topUsersViewer = new TopUsersViewer(this.uaToolbox.filterStream);
    this.viewers.push(topUsersViewer);

    let usageViewer = new UaFilterableQueryViewer(
        this.uaToolbox.filterStream,
        'Usage',
        'Usage',
        (t: DG.DataFrame) => DG.Viewer.scatterPlot(t, UaQueryViewer.defaultChartOptions).root
    );
    this.viewers.push(usageViewer);

    let topPackagesByUsers = new UaFilterableQueryViewer(
        this.uaToolbox.filterStream,
        'Packages By Users',
        'TopPackagesByUsers',
        (t: DG.DataFrame) => DG.Viewer.scatterPlot(t, UaQueryViewer.defaultChartOptions).root
    );
    this.viewers.push(topPackagesByUsers);

    this.root.append(ui.divV([
      ui.divH([ui.block([topUsersViewer.root])]),
      ui.divH([ui.block([usageViewer.root])]),
      ui.divH([ui.block([uniqueUsersViewer.root])])
    ]));

  }
}