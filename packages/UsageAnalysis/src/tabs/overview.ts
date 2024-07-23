import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {PackagesView} from './packages';
import {UaFilter} from '../filter';

export class OverviewView extends UaView {
  expanded: {[key: string]: boolean} = {f: true, l: true};

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Overview';
  }

  async initViewers(path?: string) : Promise<void> {
    this.root.className = 'grok-view ui-box';
    const uniqueUsersViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'UniqueUsers',
      activated: true,
      queryName: 'UniqueUsersOverview',
      createViewer: (t: DG.DataFrame) => {
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

    let df: Promise<DG.DataFrame>;
    let packagesSelection: DG.BitSet;
    let usersSelection: DG.BitSet;

    function query(filter: UaFilter) {
      df = grok.functions.call('UsageAnalysis:PackagesUsageOverview', {...filter});
      df.then((df) => {
        packagesSelection = DG.BitSet.create(df.rowCount, (_) => true);
        usersSelection = DG.BitSet.create(df.rowCount, (_) => true);
      });
    }
    query(this.uaToolbox.getFilter());
    this.uaToolbox.filterStream.subscribe( (filter) => {
      query(filter);
    });
    const packageStatsViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'PackageStats',
      activated: true,
      getDataFrame: () => df,
      createViewer: (t: DG.DataFrame) => {
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
          'rowSource': DG.RowSet.All,
        });
        let skipEvent: boolean = false;
        viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
          skipEvent = true;
          packagesSelection.init((i) => viewer.dataFrame.get('package', i) === args.args.options.categories[0]);
          userStatsViewer.viewer!.props.title = getUsersViewerName(args.args.options.categories[0]);
          userStatsViewer.viewer!.props.rowSource = DG.RowSet.Filtered;
          viewer.dataFrame.filter.copyFrom(usersSelection).and(packagesSelection);
          PackagesView.showSelectionContextPanel(viewer.dataFrame, this.uaToolbox, 'Overview', {showDates: false});
        });
        viewer.root.onclick = (_) => {
          //resetViewers(skipEvent, viewer.table);
          if (skipEvent) {
            skipEvent = false;
            return;
          }
          if (userStatsViewer.viewer!.props.rowSource != DG.RowSet.All) {
            packagesSelection.setAll(true);
            userStatsViewer.viewer!.props.rowSource = DG.RowSet.All;
            userStatsViewer.viewer!.props.title = getUsersViewerName();
          } else if (packageStatsViewer.viewer!.props.rowSource != DG.RowSet.All) {
            usersSelection.setAll(true);
            packageStatsViewer.viewer!.props.rowSource = DG.RowSet.All;
            packageStatsViewer.viewer!.props.title = getPackagesViewerName();
          } else {
            skipEvent = false;
            return viewer;
          }
          viewer.dataFrame.filter.copyFrom(usersSelection).and(packagesSelection);
          viewer.dataFrame.selection.copyFrom(viewer.dataFrame.filter);
          PackagesView.showSelectionContextPanel(viewer.dataFrame, this.uaToolbox, 'Overview', {showDates: false});
          skipEvent = false;
        };
        return viewer;
      }});


    const userStatsViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'UserStats',
      activated: true,
      getDataFrame: () => df,
      createViewer: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.barChart(t, {
          'valueColumnName': 'count',
          'valueAggrType': 'sum',
          'barSortType': 'by value',
          'barSortOrder': 'desc',
          'showValueAxis': false,
          'showValueSelector': false,
          'splitColumnName': 'user',
          'showCategoryValues': false,
          'showCategorySelector': false,
          'showStackSelector': false,
          'title': getUsersViewerName(),
          'legendVisibility': DG.VisibilityMode.Never,
          'onClick': DG.RowGroupAction.Select,
          'rowSource': DG.RowSet.All
        });
        let skipEvent: boolean = false;
        viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
          skipEvent = true;
          usersSelection.init((i) => viewer.dataFrame.get('user', i) == args.args.options.categories[0]);
          packageStatsViewer.viewer!.props.title = getPackagesViewerName(args.args.options.categories[0]);
          packageStatsViewer.viewer!.props.rowSource = DG.RowSet.Filtered;
          viewer.dataFrame.filter.copyFrom(usersSelection).and(packagesSelection);
          viewer.dataFrame.selection.copyFrom(viewer.dataFrame.filter);
          PackagesView.showSelectionContextPanel(viewer.dataFrame, this.uaToolbox, 'Overview', {showDates: false});
        });
        viewer.root.onclick =(_) => {
          //resetViewers(skipEvent, viewer.table);
          if (skipEvent) {
            skipEvent = false;
            return;
          }
          if (packageStatsViewer.viewer!.props.rowSource != DG.RowSet.All) {
            usersSelection.setAll(true);
            packageStatsViewer.viewer!.props.rowSource = DG.RowSet.All;
            packageStatsViewer.viewer!.props.title = getPackagesViewerName();
          } else if (userStatsViewer.viewer!.props.rowSource != DG.RowSet.All) {
            packagesSelection.setAll(true);
            userStatsViewer.viewer!.props.rowSource = DG.RowSet.All;
            userStatsViewer.viewer!.props.title = getUsersViewerName();
          } else {
            skipEvent = false;
            return viewer;
          }
          viewer.dataFrame.filter.copyFrom(usersSelection).and(packagesSelection);
          viewer.dataFrame.selection.copyFrom(viewer.dataFrame.filter);
          PackagesView.showSelectionContextPanel(viewer.dataFrame, this.uaToolbox, 'Overview', {showDates: false});
          skipEvent = false;
        };
        return viewer;
      }});


    this.viewers.push(uniqueUsersViewer);
    this.viewers.push(packageStatsViewer);
    this.viewers.push(userStatsViewer);

    const cardsView = ui.div([
    ], {classes: 'ua-cards'});

    const counters: {[key: string]: string} = {
      'Active users': 'UsageAnalysis:UniqueUsersCount',
      'New users': 'UsageAnalysis:NewUsersCount',
      'Sessions': 'UsageAnalysis:SessionsCount',
      // 'Views': 'UsageAnalysis:ViewsCount',
      'Connections': 'UsageAnalysis:ConnectionsCount',
      'Queries': 'UsageAnalysis:QueriesCount',
    };

    const refresh = (filter: UaFilter): void => {
      cardsView.textContent = '';
      for (const k of Object.keys(counters)) {
        cardsView.append(ui.div([ui.divText(k), ui.wait(async () => {
          const fc: DG.DataFrame = await grok.functions.call(counters[k], filter);
          const valuePrev = fc.get('count1', 0);
          const valueNow = fc.get('count2', 0);
          const d = valueNow - valuePrev;
          return ui.div([ui.divText(`${valueNow}`),
            ui.divText(`${d}`, {classes: d > 0 ? 'ua-card-plus' : d < 0 ? 'ua-card-minus' : ''})]);
        })], 'ua-card'));
      }
    };
    refresh(this.uaToolbox.getFilter());
    this.uaToolbox.filterStream.subscribe( (filter) => {
      refresh(filter);
    });

    this.root.append(ui.splitH([
      ui.splitV([
        ui.box(cardsView, {style: {maxHeight: '100px'}}),
        ui.box(uniqueUsersViewer.root, {style: {maxHeight: '250px'}}),
        ui.splitH([packageStatsViewer.root, userStatsViewer.root]),
      ]),
    ]));
  }
}
