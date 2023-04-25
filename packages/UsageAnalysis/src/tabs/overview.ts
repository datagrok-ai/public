// import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {PackagesView} from './packages';

export class OverviewView extends UaView {
  expanded: {[key: string]: boolean} = {f: true, l: true};

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Overview';
  }

  async initViewers() : Promise<void> {
    this.root.className = 'grok-view ui-box';
    const uniqueUsersViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'UniqueUsers',
      queryName: 'UniqueUsersOverview',
      viewerFunction: (t: DG.DataFrame) => {
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
      }});

    const getPackagesViewerName = (user?: string) => {
      if (user === undefined)
        return 'Packages unique users';
      else
        return `Package used by ${user}`;
    };

    const getUsersViewerName = (p?: string) => {
      if (p === undefined)
        return 'Users activity';
      else
        return `Users activity in ${p} package`;
    };

    function resetViewers(skipEvent: boolean, table: DG.DataFrame) {
      if (skipEvent)
        return;
      packageStatsViewer.viewer!.props.title = getPackagesViewerName();
      userStatsViewer.viewer!.props.title = getUsersViewerName();
      table.selection.setAll(false);
      table.filter.setAll(true);
    }

    const packageStatsViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'PackageStats',
      queryName: 'PackagesUsageOverview',
      viewerFunction: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.barChart(t, {
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
          'title': getPackagesViewerName(),
          'rowSource': 'Filtered',
          'filter': '${package} != ""',
          'onClick': 'Filter',
        });
        let skipEvent: boolean = false;
        viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
          t.selection.copyFrom(t.filter);
          PackagesView.showSelectionContextPanel(t, this.uaToolbox, this.expanded, 'Overview');
          userStatsViewer.viewer!.props.title = getUsersViewerName(args.args.categories[0]);
          skipEvent = true;
          console.log('cat 1');
        });
        viewer.root.onclick =(me) => {
          resetViewers(skipEvent, viewer.table);
          skipEvent = false;
          console.log('click 1');
        };
        return viewer;
      }});


    const userStatsViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'UserStats',
      getDataFrame: async () => await packageStatsViewer.dataFrame!,
      viewerFunction: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.barChart(t, {
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
          'title': getUsersViewerName(),
          'legendVisibility': 'Never',
          'rowSource': 'Filtered',
          'onClick': 'Filter',
        });
        let skipEvent: boolean = false;
        viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
          t.selection.copyFrom(t.filter);
          PackagesView.showSelectionContextPanel(t, this.uaToolbox, this.expanded, 'Overview');
          packageStatsViewer.viewer!.props.title = getPackagesViewerName(args.args.categories[0]);
          skipEvent = true;
          console.log('cat 2');
        });
        viewer.root.onclick =(me) => {
          resetViewers(skipEvent, viewer.table);
          console.log('click 2');
          skipEvent = false;
        };
        return viewer;
      }});


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
