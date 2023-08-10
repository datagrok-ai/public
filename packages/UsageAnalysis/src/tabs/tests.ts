import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';

export class TestsView extends UaView {
  loader = ui.div([ui.loader()], 'grok-wait');
  static filters: HTMLElement = ui.box();

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Tests';
  }

  async initViewers(): Promise<void> {
    this.root.appendChild(this.loader);
    const df: DG.DataFrame = await grok.data.query('UsageAnalysis:TestsToday');
    df.getCol('id').name = '~id';
    // todo: fix duplication
    TestsView.filters.appendChild(DG.Viewer.filters(df).root);
    const dfMonth: DG.DataFrame = await grok.data.query('UsageAnalysis:TestsMonth');
    df.getCol('status').colors.setCategorical({passed: '#3CB173', failed: '#EB6767', skipped: '#FFA24A'});
    const slp = DG.Viewer.fromType('Line chart', dfMonth,
      {...historyStyle, ...{yColumnNames: ['passed'], lineColor: 4282167667, description: 'Passed'}});
    const slf = DG.Viewer.fromType('Line chart', dfMonth,
      {...historyStyle, ...{yColumnNames: ['failed'], lineColor: 4293617511, description: 'Failed'}});
    const sls = DG.Viewer.fromType('Line chart', dfMonth,
      {...historyStyle, ...{yColumnNames: ['skipped'], lineColor: 4294943306, description: 'Skipped'}});
    this.root.removeChild(this.loader);
    df.onCurrentRowChanged.subscribe(async () => {
      const acc = DG.Accordion.create();
      const history: DG.DataFrame = await grok.data.query('UsageAnalysis:ScenarioHistory', {id: df.currentRow.get('~id')});
      acc.addPane('History', () => ui.waitBox(async () => {
        history.getCol('status').colors.setCategorical({passed: '#3CB173', failed: '#EB6767', skipped: '#FFA24A'});
        const grid = history.plot.grid();
        // grid.col('status')!.isTextColorCoded = true;
        grid.col('date')!.format = 'MM/dd/yyyy HH:mm:ss';
        return grid.root;
      }), true);
      acc.addPane('Execution time', () => ui.waitBox(async () => {
        const lct = DG.Viewer.fromType('Line chart', history, execTimeStyle);
        return lct.root;
      }), true);
      grok.shell.o = acc.root;
    });

    const cardsView = ui.div([], {classes: 'ua-cards'});
    const counters = ['Passed', 'Failed', 'Skipped'];
    cardsView.textContent = '';
    const cardsDf: DG.DataFrame = await grok.data.query('UsageAnalysis:TestsCount');
    for (let i = 0; i < 3; i++) {
      cardsView.append(ui.div([ui.divText(counters[i]), ui.wait(async () => {
        const valuePrev = cardsDf.get(counters[i], 0) ?? 0;
        const valueNow = cardsDf.get(counters[i], 1) ?? 0;
        const d = valueNow - valuePrev;
        return ui.div([ui.divText(`${valueNow}`),
          ui.divText(`${d}`, {classes: d > 0 ? 'ua-card-plus' : d < 0 ? 'ua-card-minus' : ''})]);
      })], 'ua-card'));
    }

    this.root.append(ui.splitV([
      ui.splitH([ui.box(cardsView), slp.root, slf.root, sls.root], {style: {maxHeight: '110px'}}),
      df.plot.grid().root,
    ]));
  }
}

const historyStyle = {
  'xColumnName': 'date',
  'aggrType': 'count',
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
  'showXAxis': false,
  'showYAxis': false,
  'showMarkers': 'Never',
  'autoLayout': false,
  'lineWidth': 2,
  'lineColoringType': 'Custom',
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
