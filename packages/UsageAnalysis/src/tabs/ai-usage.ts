import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaToolbox} from '../ua-toolbox';
import {UaFilter} from '../filter';
import {queries} from '../package-api';
import '../../css/usage_analysis.css';

const CARDS = ['Requests', 'Users', 'Input tokens', 'Output tokens', 'Total tokens', 'Estimated cost'];
const GRID_COLUMNS = ['time', 'user', 'groups', 'feature', 'model', 'status', 'total_tokens', 'cost_usd'];

export class AiUsageView extends UaView {
  private requests?: UaFilterableQueryViewer;

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'AI Usage';
  }

  async exportFiles(): Promise<DG.FileInfo[]> {
    return UaView.csvFiles({'requests': this.requests?.errorDiv == null ? this.requests?.viewer?.dataFrame : null});
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box';
    const cards = ui.div([], {classes: 'ua-cards ua-ai-cards'});
    const showCards = (t: Promise<DG.DataFrame>) => {
      const summary = t.then((t) => AiUsageView.summary(t));
      cards.textContent = '';
      for (const label of CARDS) {
        cards.append(ui.div([ui.divText(label), ui.wait(async () =>
          ui.div((await summary)[label].map((s) => ui.divText(s))))], 'ua-card'));
      }
    };

    let df: Promise<DG.DataFrame>;
    let groups: Promise<DG.DataFrame>;
    this.uaToolbox.filterStream.subscribe((filter: UaFilter) => {
      df = queries.aiUsage(filter.date!, filter.groups).then((t) => {
        t.col('groups')!.setTag(DG.TAGS.MULTI_VALUE_SEPARATOR, ', ');
        t.col('time')!.meta.format = 'yyyy-MM-dd HH:mm';
        t.onFilterChanged.subscribe(() => showCards(Promise.resolve(t)));
        t.onCurrentRowChanged.subscribe(() => AiUsageView.showRequest(t));
        return t;
      });
      groups = df.then((t) => AiUsageView.groupsFrame(t));
      showCards(df);
    });

    const viewer = (name: string, createViewer: (t: DG.DataFrame) => DG.Viewer, getDataFrame = () => df) =>
      new UaFilterableQueryViewer({filterSubscription: this.uaToolbox.filterStream, name, getDataFrame, createViewer});

    const perDay = viewer('Cost per day', (t) => DG.Viewer.lineChart(t, {
      xColumnName: 'day',
      yColumnNames: ['cost_usd'],
      yAggrTypes: ['sum'],
      splitColumnNames: ['model'],
      chartTypes: ['Stacked Bar Chart'],
      legendPosition: 'Top',
      title: 'Cost per day',
    }));
    const byUser = viewer('Cost by user', (t) => AiUsageView.costBarChart(t, 'user', 'Cost by user', {
      onClick: DG.RowGroupAction.Filter,
      rowSource: DG.RowSet.All,
    }));
    const byGroup = viewer('Cost by group', (t) => AiUsageView.costBarChart(t, 'group', 'Cost by group'), () => groups);
    const filters = viewer('Filters', (t) => DG.Viewer.filters(t,
      {columnNames: ['status', 'feature', 'model', 'groups', 'user']}));
    const requests = viewer('Requests', (t) => {
      const grid = DG.Viewer.grid(t, {
        showRowHeader: false,
        showColumnGridlines: false,
        allowRowSelection: false,
        allowBlockSelection: false,
        showCurrentCellOutline: false,
      });
      grid.columns.setOrder(GRID_COLUMNS);
      for (const name of t.columns.names())
        grid.col(name)!.visible = GRID_COLUMNS.includes(name);
      const widths: {[name: string]: number} = {time: 140, user: 110, groups: 140, feature: 90, model: 130, status: 60,
        total_tokens: 90, cost_usd: 70};
      for (const [name, width] of Object.entries(widths))
        grid.col(name)!.width = width;
      return grid;
    });

    filters.root.classList.add('ua-filters');
    requests.root.classList.add('ui-panel');
    this.requests = requests;
    this.viewers.push(perDay, byUser, byGroup, filters, requests);
    const breakdown = ui.tabControl({'By user': byUser.root, 'By group': byGroup.root});
    this.root.append(ui.splitV([
      ui.box(cards, {style: {maxHeight: '96px'}}),
      ui.box(ui.splitH([perDay.root, breakdown.root]), {style: {maxHeight: '290px'}}),
      ui.splitH([ui.box(filters.root, {style: {maxWidth: '230px'}}), requests.root]),
    ]));
  }

  static costBarChart(t: DG.DataFrame, splitColumnName: string, title: string, options: object = {}): DG.Viewer {
    return DG.Viewer.barChart(t, {
      splitColumnName,
      valueColumnName: 'cost_usd',
      valueAggrType: 'sum',
      stackColumnName: '',
      orientation: 'horizontal',
      barSortType: 'by value',
      barSortOrder: 'desc',
      title,
      ...options,
    });
  }

  static summary(t: DG.DataFrame): {[label: string]: string[]} {
    const rows = Array.from(t.filter.getSelectedIndexes());
    const sum = (name: string) => {
      const col = t.col(name)!;
      return rows.reduce((s, i) => col.isNone(i) ? s : s + col.get(i), 0);
    };
    const priced = rows.filter((i) => !t.col('cost_usd')!.isNone(i)).length;
    const n = (x: number) => DG.format(x, 'compact');
    return {
      'Requests': [n(rows.length)],
      'Users': [n(new Set(rows.map((i) => t.get('user', i))).size)],
      'Input tokens': [n(sum('input_tokens'))],
      'Output tokens': [n(sum('output_tokens'))],
      'Total tokens': [n(sum('total_tokens'))],
      'Estimated cost': [`$${n(sum('cost_usd'))}`, priced < rows.length ? `${Math.round(priced / rows.length * 100)}% priced` : ''],
    };
  }

  static groupsFrame(t: DG.DataFrame): DG.DataFrame {
    const rows: number[] = [];
    const names: string[] = [];
    for (let i = 0; i < t.rowCount; i++) {
      for (const name of (t.get('groups', i) as string | null)?.split(', ') ?? []) {
        rows.push(i);
        names.push(name);
      }
    }
    const groups = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('group', names),
      DG.Column.fromStrings('user', rows.map((i) => t.get('user', i))),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'total_tokens', rows.map((i) => Number(t.get('total_tokens', i) ?? 0))),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'cost_usd', rows.map((i) => t.get('cost_usd', i))),
    ]);
    t.onFilterChanged.subscribe(() => groups.filter.init((j) => t.filter.get(rows[j])));
    return groups;
  }

  static showRequest(t: DG.DataFrame): void {
    const i = t.currentRowIdx;
    if (i < 0)
      return;
    const s = (name: string, row: number = i) => t.col(name)!.getString(row);
    const acc = ui.accordion();
    acc.addPane('Details', () => ui.tableFromMap({'User': s('user'), 'Groups': s('groups'), 'Feature': s('feature'),
      'Model': s('model'), 'Status': s('status'), 'Cost, $': s('cost_usd'), 'Duration, ms': s('duration_ms'),
      'Chat': s('chat_id')}), true);
    acc.addPane('Tokens', () => ui.tableFromMap({'Input': s('input_tokens'), 'Output': s('output_tokens'),
      'Cache read': s('cache_read_tokens'), 'Cache write': s('cache_write_tokens'), 'Total': s('total_tokens')}), true);
    const chat = t.get('chat_id', i);
    if (chat != null) {
      acc.addPane('Chat', () => {
        const others = [...Array(t.rowCount).keys()].filter((j) => j !== i && t.get('chat_id', j) === chat);
        return others.length === 0 ? ui.divText('No other requests') : ui.table(others,
          (j) => [ui.link(s('time', j), () => t.currentRowIdx = j), s('total_tokens', j), s('cost_usd', j)],
          ['Time', 'Tokens', 'Cost, $']);
      }, true);
    }
    const content = queries.aiTurnContent(s('turn_id'));
    const text = (name: string) => ui.wait(async () => {
      const value = (await content).get(name, 0);
      return value ? ui.divText(value, 'ua-ai-text') : ui.divText('Not logged', 'ua-ai-text-missing');
    });
    acc.addPane('Prompt', () => text('prompt'), true);
    acc.addPane('Answer', () => text('response'), true);
    acc.addPane('Tool calls', () => text('tool_calls'), true);
    grok.shell.o = acc.root;
  }
}
