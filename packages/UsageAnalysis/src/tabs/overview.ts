import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {Filter, UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaDataFrameQueryViewer} from '../viewers/ua-data-frame-query-viewer';
import {TopPackagesViewer} from '../drilldown_viewers/events/top-packages-viewer';
import {UaQueryViewer} from "../viewers/abstract/ua-query-viewer";

export class OverviewView extends UaView {

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Overview';
  }

  async initViewers() : Promise<void> {

    this.root.className = 'grok-view ui-box';
    const uniqueUsersViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'UniqueUsers',
      'UniqueUsersOverview',
      (t: DG.DataFrame) => {
        return DG.Viewer.lineChart(t, {
          'overviewColumnName': 'date',
          'xColumnName': 'date',
          'showXSelector': false,
          'yColumnNames': ['count'],
          'showYSelectors': false,
          'showAggrSelectors': false,
          'showSplitSelector': false,
          'showMarkers': 'Never',
          'chartTypes': ['Line Chart'],
          'title': 'Unique users',
        });
      }, null, null);

    const packageStatsViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'PackageStats',
      'PackagesUsage',
      (t: DG.DataFrame) => {
        return DG.Viewer.barChart(t, {
          'valueColumnName': 'user',
          'valueAggrType': 'unique',
          'barSortType': 'by value',
          'barSortOrder': 'desc',
          'showValueAxis': false,
          'showValueSelector': false,
          'splitColumnName': 'package',
          'showCategoryValues': false,
          'showCategorySelector': false,
          'stackColumnName': '',
          'showStackSelector': false,
          'title': 'Packages activity',
        });
      }, null, null);

    const userStatsViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'UserStats',
      'UniqueUsersStats',
      (t: DG.DataFrame) => {
        return DG.Viewer.barChart(t, {
          'valueColumnName': 'count',
          'valueAggrType': 'sum',
          'barSortType': 'by value',
          'barSortOrder': 'desc',
          'showValueAxis': false,
          'showValueSelector': false,
          'splitColumnName': 'name',
          'showCategoryValues': false,
          'showCategorySelector': false,
          'showStackSelector': false,
          'title': 'Users activity',
          'legendVisibility': 'Never',
        });
      }, null, null);


    this.viewers.push(uniqueUsersViewer);
    this.viewers.push(packageStatsViewer);
    this.viewers.push(userStatsViewer);
    //this.root.append(uniqueUsersViewer.root);
    this.root.append(ui.splitH([
      ui.splitV([
        ui.box(uniqueUsersViewer.root, {style: {maxHeight: '250px'}}),
        ui.splitH([packageStatsViewer.root, userStatsViewer.root]),
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


}
