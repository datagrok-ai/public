import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';

const colors = {'passed': '#3CB173', 'failed': '#EB6767', 'skipped': '#FFA24A'};

export class TestsView extends UaView {
  loader = ui.div([ui.loader()], 'grok-wait');
  static filters: HTMLElement = ui.box();

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Tests';
  }

  async initViewers(): Promise<void> {
    // this.root.appendChild(this.loader);
    const df: DG.DataFrame = await grok.data.query('UsageAnalysis:TestsToday');
    df.getCol('id').name = '~id';
    if (TestsView.filters.children.length) TestsView.filters = ui.box();
    const viewer = DG.Viewer.filters(df, filtersStyle);
    const filtersGroup = new DG.FilterGroup(viewer.dart);
    TestsView.filters.appendChild(viewer.root);
    const dfMonth: DG.DataFrame = await grok.data.query('UsageAnalysis:TestsMonth');
    df.getCol('status').colors.setCategorical(colors);
    dfMonth.getCol('status').colors.setCategorical(colors);
    const chart = DG.Viewer.fromType('Line chart', dfMonth, historyStyle);
    // this.root.removeChild(this.loader);
    df.onCurrentRowChanged.subscribe(async () => {
      const row = df.currentRow;
      const acc = DG.Accordion.create();
      const history: DG.DataFrame = await grok.data.query('UsageAnalysis:ScenarioHistory', {id: row.get('~id')});
      acc.addPane('Details', () => {
        const table = ui.tableFromMap({
          Test: row.get('test'),
          Category: row.get('category'),
          Package: row.get('package'),
        });
        table.style.userSelect = 'text';
        return table;
      }, true);
      acc.addPane('History', () => ui.waitBox(async () => {
        history.getCol('status').colors.setCategorical(colors);
        const grid = history.plot.grid();
        // grid.col('status')!.isTextColorCoded = true;
        grid.col('date')!.format = 'MM/dd/yy'; // MM/dd/yyyy HH:mm:ss
        return grid.root;
      }), true);
      acc.addPane('Execution time', () => ui.waitBox(async () => {
        const lct = DG.Viewer.fromType('Line chart', history, execTimeStyle);
        return lct.root;
      }), true);
      grok.shell.o = acc.root;
    });

    const cardsView = ui.div([], {classes: 'ua-cards'});
    const counters = ['passed', 'failed', 'skipped'];
    cardsView.textContent = '';
    const cardsDfP: Promise<DG.DataFrame> = grok.data.query('UsageAnalysis:TestsCount');
    for (let i = 0; i < 3; i++) {
      const c = counters[i];
      const card = ui.div([ui.divText(c), ui.wait(async () => {
        const cardsDf = await cardsDfP;
        const valuePrev = cardsDf.get(c, 0) ?? 0;
        const valueNow = cardsDf.get(c, 1) ?? 0;
        const d = valueNow - valuePrev;
        return ui.div([ui.divText(`${valueNow}`),
          ui.divText(`${d}`, {classes: d > 0 ? 'ua-card-plus' : '', style: {color: getColor(c, d)}})]);
      })], 'ua-card ua-test-card');
      card.addEventListener('click', () => {
        console.log(filtersGroup.table);
        filtersGroup.updateOrAdd({
          type: DG.FILTER_TYPE.CATEGORICAL,
          column: 'status',
          selected: [c],
        });
        console.log(df.filter.falseCount);
      });
      cardsView.append(card);
    }

    this.root.append(ui.splitV([
      ui.splitH([ui.box(cardsView), chart.root], {style: {maxHeight: '150px'}}),
      df.plot.grid().root,
    ], null, true));
  }
}

function getColor(status: string, d: number): string {
  switch (status) {
  case 'passed':
    return d > 0 ? 'var(--green-2)' : d < 0 ? 'var(--red-3)' : '';
  case 'failed':
    return d < 0 ? 'var(--green-2)' : d > 0 ? 'var(--red-3)' : '';
  default:
    return '';
  }
}

const filtersStyle = {
  columnNames: ['type', 'package', 'status', 'ms', 'category'],
};

const historyStyle = {
  'xColumnName': 'date',
  'splitColumnName': 'status',
  'legendVisibility': 'Never',
  // 'aggrType': 'count',
  'innerChartMarginTop': 0,
  'innerChartMarginBottom': 0,
  'outerChartMarginTop': 0,
  'outerChartMarginBottom': 0,
  'outerChartMarginLeft': 0,
  'outerChartMarginRight': 0,
  // 'yGlobalScale': false,
  'showTopPanel': false,
  // 'showMouseOverRowLine': false,
  'showXSelector': false,
  'showYSelectors': false,
  'showAggrSelectors': false,
  'showSplitSelector': false,
  // 'showXAxis': false,
  // 'showYAxis': false,
  // 'showMarkers': 'Never',
  'autoLayout': false,
  // 'lineWidth': 2,
  // 'lineColoringType': 'Custom',
  'chartTypes': ['Area Chart'],
};

const execTimeStyle = {
  x: 'date',
  y: 'ms',
  // showMouseOverRowLine: false,
  showXSelector: false,
  showYSelectors: false,
  showAggrSelectors: false,
  showSplitSelector: false,
  showXAxis: true,
  showYAxis: true,
  showMarkers: 'Never',
  lineWidth: 2,
  autoLayout: false,
};
