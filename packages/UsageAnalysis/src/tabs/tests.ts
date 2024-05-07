import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {getDate, colors} from '../utils';

export const filters = ui.box();
filters.id = 'ua-tests-filters';
const counters = ['passed', 'failed', 'skipped'];

export class TestsView extends UaView {
  loader = ui.div([ui.loader()], 'grok-wait');
  cardFilter: string | null = null;
  platformFilter: boolean = false;
  filterGroup: DG.FilterGroup | undefined = undefined;
  grid?: DG.Grid;
  leftDate = ui.label(getDate(new Date()));
  rightDate = ui.label('');
  leftDf?: DG.DataFrame;
  rightDf?: DG.DataFrame;
  current = this.leftDate;
  cardsView: HTMLDivElement = ui.div([], {classes: 'ua-cards'});

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Tests';
    filters.style.display = 'none';
    this.rightDate.style.marginRight = '12px';
  }

  async initViewers(path?: string): Promise<void> {
    // Stacked Bar Chart
    const chart = ui.wait(async () => {
      const dfMonth: DG.DataFrame = await grok.functions.call('UsageAnalysis:TestsMonth');
      dfMonth.getCol('status').colors.setCategorical(colors);
      dfMonth.onSelectionChanged.subscribe(async () => {
        if (!this.grid) return;
        ui.setUpdateIndicator(this.grid.root);
        const date = dfMonth.get('date', dfMonth.selection.getSelectedIndexes()[0]);
        this.updateCards(date);
        this.current.innerText = getDate(date.toDate());
        if (!!this.rightDate.innerText)
          done.disabled = false;
        const df: DG.DataFrame = await grok.functions.call('UsageAnalysis:TestsToday', {'date': date});
        this.updateGrid(df);
        ui.setUpdateIndicator(this.grid.root, false);
      });
      const areaChart = DG.Viewer.fromType('Line chart', dfMonth, historyStyle);
      return areaChart.root;
    });
    chart.style.marginRight = '12px';

    // Table
    const grid = ui.wait(async () => {
      const d = getDate(new Date());
      const df: DG.DataFrame = await grok.functions.call('UsageAnalysis:TestsToday', {'date': d});
      this.updateGrid(df);
      return this.grid!.root;
    });

    // Cards
    this.updateCards(getDate(new Date()));

    // Compare section
    const compare = ui.button('compare to', () => {
      this.leftDf = this.grid?.dataFrame;
      this.current = this.rightDate;
      compare.disabled = true;
    });
    let icon = ui.iconFA('external-link');
    icon.classList.replace('fal', 'fas');
    const done = ui.button(icon, () => {
      const res = grok.data.compareTables(this.leftDf!, this.grid?.dataFrame!,
        ['package', 'test', 'category'], ['package', 'test', 'category'],
        ['status', 'result', 'ms'], ['status', 'result', 'ms']);
      grok.shell.addTableView(res.diffTable);
    }, 'Compare tables');
    done.disabled = true;
    icon = ui.iconFA('undo');
    icon.classList.replace('fal', 'fas');
    const reset = ui.button(icon, () => {
      this.current = this.leftDate;
      this.rightDate.innerText = '';
      compare.disabled = false;
      done.disabled = true;
    }, 'Reset');
    const compareSection = ui.divH([this.leftDate, compare, this.rightDate, done, reset], 'ua-tt-compare');

    // Platfrom filter
    const platform = ui.button('platform', () => {
      this.filterGroup?.updateOrAdd({
        type: DG.FILTER_TYPE.CATEGORICAL,
        column: 'package',
        selected: this.platformFilter ? null : ['ApiTests', 'UiTests', 'Dbtests', 'DevTools'],
      });
      this.platformFilter = !this.platformFilter;
    });

    const leftSide = ui.divV([
      ui.box(this.cardsView, {style: {flexGrow: 0, flexBasis: '35%'}}),
      compareSection,
    ]);
    // leftSide.parentElement!.style.maxWidth = '400px';

    this.root.append(ui.splitV([
      ui.splitH([ui.divH([leftSide, platform]), chart], {style: {height: '150px'}}),
      grid], null, true));
  }

  updateCards(date: any): void {
    this.cardsView.innerHTML = '';
    const cardsDfP: Promise<DG.DataFrame> = grok.functions.call('UsageAnalysis:TestsCount', {'date': date});
    for (let i = 0; i < 3; i++) {
      const c = counters[i];
      const card = ui.div([ui.divText(c), ui.wait(async () => {
        const cardsDf = await cardsDfP;
        const valuePrev = cardsDf.get(c, 0);
        const valueNow = cardsDf.get(c, 1);
        const d = valueNow - valuePrev;
        return ui.div([ui.divText(`${valueNow}`),
          ui.divText(`${d}`, {classes: d > 0 ? 'ua-card-plus' : '', style: {color: getColor(c, d)}})]);
      })], 'ua-card ua-test-card');
      card.addEventListener('click', () => {
        this.cardFilter = this.cardFilter === c ? null : c;
        this.filterGroup?.updateOrAdd({
          type: DG.FILTER_TYPE.CATEGORICAL,
          column: 'status',
          selected: this.cardFilter ? [c] : null,
        });
      });
      this.cardsView.append(card);
    }
  }

  updateGrid(df: DG.DataFrame): void {
    df.getCol('status').colors.setCategorical(colors);
    df.getCol('id').name = '~id';
    if (filters.children.length)
      filters.innerHTML = '';
    const filters_ = DG.Viewer.filters(df, filtersStyle);
    this.filterGroup = new DG.FilterGroup(filters_.dart);
    filters.appendChild(filters_.root);
    if (this.grid === undefined)
      this.grid = df.plot.grid();
    else
      this.grid.dataFrame = df;
    df.onCurrentRowChanged.subscribe(async () => {
      const row = df.currentRow;
      const acc = DG.Accordion.create();
      acc.root.style.width = '100%';
      let history: DG.DataFrame = await grok.functions.call('UsageAnalysis:ScenarioHistory', {id: row.get('~id')});
      history.getCol('id').name = '~id';
      history.getCol('uid').name = '~uid';
      const resName = history.getCol('res_name');
      if (resName.stats.valueCount) {
        const newHistory = history.groupBy(['date', 'status', 'ms', 'result', '~uid']).add('key', '~id').aggregate();
        const categories = resName.categories.filter((c) => c);
        categories.forEach((c) => newHistory.columns.addNewString(c));
        for (const row of newHistory.rows) {
          const id = row.get('~id');
          const df = history.rows.match({'~id': `= ${id}`}).toDataFrame();
          for (const r of df.rows) {
            const name = r.get('res_name');
            if (name)
              newHistory.set(name, row.idx, r.get('res_value'));
          }
        }
        newHistory.columns.byNames(categories).forEach((col) => col.name = col.name.replace('result.', ''));
        history = newHistory;
      } else {
        history.columns.remove('res_name');
        history.columns.remove('res_value');
      }
      const grid = history.plot.grid();
      grid.sort(['date'], [false]);
      acc.addPane('Details', () => {
        const table = ui.tableFromMap({
          Test: row.get('test'),
          Category: row.get('category'),
          Package: row.get('package'),
        });
        table.style.userSelect = 'text';
        return table;
      }, true);
      const hPane = acc.addPane('History', () => ui.waitBox(async () => {
        history.getCol('status').colors.setCategorical(colors);
        // grid.col('status')!.isTextColorCoded = true;
        grid.col('date')!.format = 'MM/dd/yy'; // MM/dd/yyyy HH:mm:ss
        grid.col('date')!.width = 70;
        grid.col('status')!.width = 19;
        grid.onCellTooltip((gc, x, y) => {
          if (gc.isColHeader) return false;
          const row = gc.tableRow;
          ui.tooltip.show(ui.tableFromMap({
            'Time': row?.get('date'),
            'Run by': ui.render(`#{x.${row?.get('~uid')}}`),
            // Package: ui.render(`#{x.${r.uid}}`),
          }), x, y);
          return true;
        });
        return grid.root;
      }), true);
      const info = hPane.root.querySelector('.d4-accordion-pane-header') as HTMLElement;
      info.firstChild?.after(ui.button('Add to workspace', () => grok.shell.addTableView(history)));
      acc.addPane('Execution time', () => ui.waitBox(async () => {
        const lct = DG.Viewer.fromType('Line chart', history, execTimeStyle);
        return lct.root;
      }), true);
      grok.shell.o = acc.root;
    });
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
  columnNames: ['type', 'version', 'package', 'status', 'ms'], // category
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
  'chartTypes': ['Stacked Bar Chart'],
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
