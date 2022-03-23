import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../../css/usage_analysis.css';
import {UaToolbox} from "../ua-toolbox";
import {UaView} from "./ua-view";
import {UaFilter} from "../filter2";
import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";
import {UaDataFrameQueryViewer} from "../viewers/ua-data-frame-query-viewer";

export class OverviewView extends UaView {

  constructor(uaToolbox: UaToolbox) {
    super('Overview', uaToolbox);
  }

  async initViewers() : Promise<void> {
    let uniqueUsersViewer = new UaFilterableQueryViewer(
        this.uaToolbox.filterStream,
        'Unique Users',
        'UniqueUsers',
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.lineChart(t, UaFilterableQueryViewer.splineStyle).root;
          viewer.style.maxHeight = '150px';
          return viewer;
        }
    );
    this.viewers.push(uniqueUsersViewer);

    let eventsViewer = new UaFilterableQueryViewer(
        this.uaToolbox.filterStream,
        'Events',
        'Events1',
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.lineChart(t, UaFilterableQueryViewer.splineStyle).root;
          viewer.style.maxHeight = '150px';
          return viewer;
        }
    );
    this.viewers.push(eventsViewer);

    let errorsViewer = new UaFilterableQueryViewer(
        this.uaToolbox.filterStream,
        'Errors',
        'Errors1',
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.lineChart(t, UaFilterableQueryViewer.splineStyle).root;
          viewer.style.maxHeight = '150px';
          return viewer;
        }
    );
    this.viewers.push(errorsViewer);

    let uniqueUsersListViewer = new UaFilterableQueryViewer(
        this.uaToolbox.filterStream,
        'Unique Users List',
        'UniqueUsersList',
        (t: DG.DataFrame) => {
          let ids = Array.from(t.getCol('id').values());
          return ui.wait(async () =>  ui.list(await grok.dapi.getEntities(ids)));
        },
        (host: HTMLElement) => {
          host.style.overflow="auto";
          host.style.height="94.5%";
        }
    );
    this.viewers.push(uniqueUsersListViewer);

    let totalUsersViewer = new UaDataFrameQueryViewer(
        'Total Users',
        'TotalUsersAndGroups',
        (t: DG.DataFrame) => {
          let list = [
            ['Total users', t.get('user_count', 0)],
            ['Total groups', t.get('group_count', 0)]
          ];

          return ui.div([ui.table(list, (item, idx) =>
              [`${item[0]}:`, item[1]]
          )]);
        }
    );
    this.viewers.push(totalUsersViewer);

    this.root.append(ui.block25([
      ui.block([totalUsersViewer.root]),
      ui.block([uniqueUsersListViewer.root]),
    ]));

    this.root.append(ui.block75([
      ui.block([uniqueUsersViewer.root]),
      ui.block([eventsViewer.root]),
      ui.block([errorsViewer.root])
    ]));
  }

}