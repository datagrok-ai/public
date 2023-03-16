import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaDataFrameQueryViewer} from '../viewers/ua-data-frame-query-viewer';
import {TopPackagesViewer} from '../drilldown_viewers/events/top-packages-viewer';

export class OverviewView extends UaView {
  static viewName = 'Overview';
  topPackagesViewer: TopPackagesViewer | null = null;

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
  }

  async initViewers() : Promise<void> {
    const uniqueUsersViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Unique Users',
      'UniqueUsers',
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.lineChart(t, UaFilterableQueryViewer.splineStyle).root;
        viewer.style.maxHeight = '150px';
        return viewer;
      },
    );
    this.viewers.push(uniqueUsersViewer);

    // let eventsViewer = new UaFilterableQueryViewer(
    //     this.uaToolbox.filterStream,
    //     'Events',
    //     'Events1',
    //     (t: DG.DataFrame) => {
    //       let viewer = DG.Viewer.lineChart(t, UaFilterableQueryViewer.splineStyle).root;
    //       viewer.style.maxHeight = '150px';
    //       return viewer;
    //     }
    // );
    // this.viewers.push(eventsViewer);

    // let errorsViewer = new UaFilterableQueryViewer(
    //     this.uaToolbox.filterStream,
    //     'Errors',
    //     'Errors1',
    //     (t: DG.DataFrame) => {
    //       let viewer = DG.Viewer.lineChart(t, UaFilterableQueryViewer.splineStyle).root;
    //       viewer.style.maxHeight = '150px';
    //       return viewer;
    //     }
    // );
    // this.viewers.push(errorsViewer);

    const uniqueUsersListViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Users',
      'UniqueUsersList',
      (t: DG.DataFrame) => {
        const ids = Array.from(t.getCol('id').values());
        return ui.wait(async () => ui.list(await grok.dapi.getEntities(ids)));
      },
      (host: HTMLElement) => {
        host.style.overflow='auto';
        host.style.height='94.5%';
      },
    );
    this.viewers.push(uniqueUsersListViewer);

    const totalUsersViewer = new UaDataFrameQueryViewer(
      'Total Users',
      'TotalUsersAndGroups',
      (t: DG.DataFrame) => {
        const list = [
          ['Total users', t.get('user_count', 0)],
          ['Total groups', t.get('group_count', 0)],
        ];

        return ui.div([ui.table(list, (item, idx) =>
          [`${item[0]}:`, item[1]],
        )]);
      },
    );
    this.viewers.push(totalUsersViewer);

    this.topPackagesViewer = new TopPackagesViewer('Packages', 'TopPackages', this.uaToolbox.filterStream);
    this.viewers.push(this.topPackagesViewer);

    this.root.className = 'grok-view ui-box';

    this.root.append(ui.splitH([
      ui.splitV([
        ui.panel([ui.h1('Groups')]),
        ui.panel([uniqueUsersListViewer.root]),
      ], {style: {maxWidth: '300px'}}),
      ui.panel([
        ui.block([uniqueUsersViewer.root]),
        ui.block([this.topPackagesViewer.root]),
      ]),
    ]));

    /*
    this.root.append(ui.block25([
      ui.block([totalUsersViewer.root]),
      ui.block([uniqueUsersListViewer.root]),
    ]));

    this.root.append(ui.block75([
      ui.block([uniqueUsersViewer.root]),
      ui.block([this.topPackagesViewer.root]),
    ]));
    */
  }

  handleUrlParams(params: Map<string, string>) : void {
    if (params.has('package'))
      this.topPackagesViewer?.categorySelected(params.get('package')!);
  }
}
