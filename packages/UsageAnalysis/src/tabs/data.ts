/*
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UaView} from './ua';
import {UaToolbox} from '../ua-toolbox';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaQueryViewer} from '../viewers/abstract/ua-query-viewer';
import {TopDataSourcesViewer} from '../drilldown_viewers/top-data-sources-viewer';
import {TopQueriesViewer} from '../drilldown_viewers/data/top-queries-viewer';
import {TopConnectionsViewer} from '../drilldown_viewers/data/top-connection-viewer';

export class DataView extends UaView {
  static viewName = 'Data';

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
  }

  async initViewers() : Promise<void> {
    const queriesViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Queries',
      'Queries1',
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.lineChart(t, UaQueryViewer.defaultChartOptions);
        viewer.root.style.maxHeight = '150px';
        return viewer;
      },
    );
    this.viewers.push(queriesViewer);

    const topQueriesViewer = new TopQueriesViewer(this.uaToolbox.filterStream);
    this.viewers.push(topQueriesViewer);

    const topConnectionsViewer = new TopConnectionsViewer(this.uaToolbox.filterStream);
    this.viewers.push(topConnectionsViewer);

    const topDataSourcesViewer = new TopDataSourcesViewer(this.uaToolbox.filterStream);
    this.viewers.push(topDataSourcesViewer);

    this.root.append(ui.divV([
      ui.divH([ui.block([queriesViewer.root])]),
      ui.divH([ui.block50([topQueriesViewer.root]), ui.block50([topConnectionsViewer.root])]),
      ui.divH([ui.block50([topDataSourcesViewer.root])]),
    ]));
  }
}
*/
