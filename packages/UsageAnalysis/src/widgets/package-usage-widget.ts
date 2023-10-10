import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import ColorHash from 'color-hash';

const colorHash = new ColorHash();

export class PackageUsageWidget extends DG.Widget {
  date: string = 'this week';
  groupsP: Promise<DG.Group[]> = grok.dapi.groups.getGroupsLookup('All users');

  constructor(pack: DG.Package) {
    const scatter = ui.box(null, {style: {maxHeight: '100px', marginTop: '10px'}});
    const line = ui.box(null, {style: {maxHeight: '70px'}});
    super(ui.box(ui.splitV([line, scatter], {classes: 'ua-widget ua-package-widget'})));

    const df = this.groupsP.then((groups) => grok.data.query('UsageAnalysis:PackagesUsage',
      {date: this.date, groups: [groups[0].id], packages: [pack.name]})).then((dataFrame: DG.DataFrame) => {
      dataFrame.tags[DG.TAGS.TOOLTIP] = 'user\ncount\ntime_start\ntime_end';
      if (dataFrame.columns.byName('count') != null)
        dataFrame.columns.byName('count').tags['format'] = '#';
      const userColumn = dataFrame.columns.byName('user');
      if (userColumn != null) {
        const users: {[key: string]: string} = {};
        userColumn.categories.forEach((u: string) => {users[u] = colorHash.hex(u);});
        userColumn.meta.colors.setCategorical(users);
      }
      return dataFrame;
    });

    scatter.appendChild(ui.waitBox(async () => {
      return DG.Viewer.scatterPlot(await df, scatterStyle).root;
    }));

    line.appendChild(ui.waitBox(async () => {
      return DG.Viewer.lineChart(await df, lineStyle).root;
    }));
  }
}

const scatterStyle = {
  innerChartMarginTop: 0,
  innerChartMarginBottom: 0,
  outerChartMarginTop: 0,
  outerChartMarginBottom: 0,
  outerChartMarginLeft: 0,
  outerChartMarginRight: 0,
  autoLayout: false,
  x: 'time_start',
  y: 'package',
  size: 'count',
  color: 'user',
  jitterSize: 5,
  markerMinSize: 10,
  markerMaxSize: 30,
  showColorSelector: false,
  showSizeSelector: false,
  showXSelector: false,
  showYSelector: false,
  showYAxis: false,
  legendVisibility: 'Never',
};

const lineStyle = {
  'aggrType': 'avg',
  'xColumnName': 'time_start',
  'yColumnNames': ['count'],
  'innerChartMarginTop': 0,
  'innerChartMarginBottom': 0,
  'outerChartMarginTop': 0,
  'outerChartMarginBottom': 0,
  'outerChartMarginLeft': 0,
  'outerChartMarginRight': 0,
  'yGlobalScale': false,
  'showTopPanel': false,
  'showMouseOverRowLine': false,
  'showXSelector': false,
  'showYSelectors': false,
  'showAggrSelectors': false,
  'showSplitSelector': false,
  'showYAxis': false,
  'showXAxis': false,
  'showMarkers': 'Never',
  'autoLayout': false,
  'lineWidth': 2,
  // 'chartTypes': ['Area Chart'],
};
