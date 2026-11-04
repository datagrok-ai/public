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
const COLUMNS = ['time', 'period', 'user', 'groups', 'feature', 'model', 'status', 'chat_id', 'turn_id', 'input_tokens',
  'output_tokens', 'cache_read_tokens', 'cache_write_tokens', 'total_tokens', 'cost_usd', 'duration_ms'] as const;
type AiColumn = typeof COLUMNS[number];
const GRID_COLUMNS: AiColumn[] = ['time', 'user', 'feature', 'model', 'status', 'total_tokens', 'cost_usd'];
// The narrowest widths that fit the pane without scrolling sideways; a wider grid scales them up.
const GRID_WIDTHS: Partial<Record<AiColumn, number>> = {time: 108, user: 86, feature: 82, model: 116, status: 46,
  total_tokens: 60, cost_usd: 52};
const col = (t: DG.DataFrame, name: AiColumn): DG.Column => t.getCol(name);

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
    this.root.style.overflow = 'hidden';
    const cards = ui.div([], {classes: 'ua-cards ua-ai-cards'});
    const showCards = (t: Promise<DG.DataFrame>) => {
      const summary = t.then((t) => AiUsageView.summary(t));
      cards.textContent = '';
      for (const label of CARDS) {
        cards.append(ui.div([ui.divText(label), ui.wait(async () =>
          ui.div((await summary)[label].map((s) => ui.divText(s))))], 'ua-card'));
      }
      summary.catch((e) => {
        cards.textContent = '';
        cards.append(ui.divText(`AI usage could not be loaded: ${e instanceof Error ? e.message : e}`, 'd4-viewer-error'));
      });
    };

    let df: Promise<DG.DataFrame>;
    let groups: Promise<DG.DataFrame>;
    this.uaToolbox.filterStream.subscribe((filter: UaFilter) => {
      df = queries.aiUsage(filter.date ?? 'this week', filter.groups).then((t) => {
        for (const name of COLUMNS)
          col(t, name);
        col(t, 'groups').setTag(DG.TAGS.MULTI_VALUE_SEPARATOR, ', ');
        col(t, 'time').meta.format = 'yyyy-MM-dd HH:mm';
        col(t, 'cost_usd').meta.format = '0.0000';
        t.onFilterChanged.subscribe(() => showCards(Promise.resolve(t)));
        t.onCurrentRowChanged.subscribe(() => AiUsageView.showRequest(t));
        return t;
      });
      groups = df.then((t) => AiUsageView.groupsFrame(t));
      showCards(df);
    });

    const viewer = (name: string, createViewer: (t: DG.DataFrame) => DG.Viewer, getDataFrame = () => df) =>
      new UaFilterableQueryViewer({filterSubscription: this.uaToolbox.filterStream, name, getDataFrame, createViewer});

    const perDay = viewer('Cost over time', (t) => DG.Viewer.lineChart(t, {
      xColumnName: 'period',
      yColumnNames: ['cost_usd'],
      yAggrTypes: ['sum'],
      splitColumnNames: ['model'],
      chartTypes: ['Stacked Bar Chart'],
      legendPosition: 'Top',
      title: 'Cost over time',
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
      for (const name of t.columns.names()) {
        const column = grid.col(name);
        if (column)
          column.visible = (GRID_COLUMNS as string[]).includes(name);
      }
      const minimum = Object.values(GRID_WIDTHS).reduce((s, w) => s + w, 0);
      ui.onSizeChanged(grid.root).subscribe(() => {
        const scale = Math.max(1, (grid.root.clientWidth - 20) / minimum);
        for (const [name, width] of Object.entries(GRID_WIDTHS)) {
          const column = grid.col(name);
          if (column)
            column.width = Math.floor(width * scale);
        }
      });
      return grid;
    });

    filters.root.classList.add('ua-filters');
    requests.root.classList.add('ui-panel');
    this.requests = requests;
    this.viewers.push(perDay, byUser, byGroup, filters, requests);
    const breakdown = ui.tabControl({'By user': byUser.root, 'By group': byGroup.root});
    this.root.append(ui.splitV([
      ui.box(cards, {style: {maxHeight: '96px'}}),
      ui.box(ui.splitH([perDay.root, breakdown.root]), {style: {maxHeight: '38%'}}),
      ui.splitH([ui.box(filters.root, {style: {maxWidth: '190px'}}), requests.root]),
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
    const sum = (name: AiColumn) => {
      const values = col(t, name);
      return rows.reduce((s, i) => values.isNone(i) ? s : s + values.get(i), 0);
    };
    const cost = col(t, 'cost_usd');
    const priced = rows.filter((i) => !cost.isNone(i)).length;
    const n = (x: number) => DG.format(x, 'compact');
    return {
      'Requests': [n(rows.length)],
      'Users': [n(new Set(rows.map((i) => col(t, 'user').get(i))).size)],
      'Input tokens': [n(sum('input_tokens'))],
      'Output tokens': [n(sum('output_tokens'))],
      'Total tokens': [n(sum('total_tokens'))],
      'Estimated cost': [`$${n(sum('cost_usd'))}`, priced < rows.length ? `${Math.round(priced / rows.length * 100)}% priced` : ''],
    };
  }

  static groupsFrame(t: DG.DataFrame): DG.DataFrame {
    const rows: number[] = [];
    const names: string[] = [];
    const column = col(t, 'groups');
    const separator = column.getTag(DG.TAGS.MULTI_VALUE_SEPARATOR);
    for (let i = 0; i < t.rowCount; i++)
      for (const name of column.isNone(i) ? [] : (column.get(i) as string).split(separator)) {
        rows.push(i);
        names.push(name);
      }
    const groups = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('group', names),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'cost_usd', rows.map((i) => col(t, 'cost_usd').get(i))),
    ]);
    t.onFilterChanged.subscribe(() => groups.filter.init((j) => t.filter.get(rows[j])));
    return groups;
  }

  static showRequest(t: DG.DataFrame): void {
    const i = t.currentRowIdx;
    if (i < 0)
      return;
    const s = (name: AiColumn, row: number = i) => col(t, name).getString(row);
    const acc = ui.accordion();
    acc.addPane('Details', () => ui.tableFromMap({'User': s('user'), 'Groups': s('groups'), 'Feature': s('feature'),
      'Model': s('model'), 'Status': s('status'), 'Cost, $': s('cost_usd'), 'Duration, ms': s('duration_ms'),
      'Chat': s('chat_id')}), true);
    acc.addPane('Tokens', () => ui.tableFromMap({'Input': s('input_tokens'), 'Output': s('output_tokens'),
      'Cache read': s('cache_read_tokens'), 'Cache write': s('cache_write_tokens'), 'Total': s('total_tokens')}), true);
    const chat = col(t, 'chat_id').get(i);
    if (chat != null) {
      acc.addPane('Chat', () => {
        const others = [...Array(t.rowCount).keys()].filter((j) => j !== i && col(t, 'chat_id').get(j) === chat);
        return others.length === 0 ? ui.divText('No other requests') : ui.table(others,
          (j) => [ui.link(s('time', j), () => t.currentRowIdx = j), s('total_tokens', j), s('cost_usd', j)],
          ['Time', 'Tokens', 'Cost, $']);
      }, true);
    }
    const turnId: string | null = col(t, 'turn_id').get(i);
    if (turnId != null) {
      const content = queries.aiTurnContent(turnId);
      const text = (name: 'prompt' | 'response' | 'tool_calls') => ui.wait(async () => {
        try {
          const value: string | null = (await content).col(name)?.get(0) ?? null;
          return value ? ui.divText(value, 'ua-ai-text') : ui.divText('Not logged', 'ua-ai-text-missing');
        }
        catch (e) {
          return ui.info(`${e instanceof Error ? e.message : e}`, 'Could not be loaded');
        }
      });
      acc.addPane('Prompt', () => text('prompt'), true);
      acc.addPane('Answer', () => text('response'), true);
      acc.addPane('Tool calls', () => text('tool_calls'), true);
    }
    grok.shell.setCurrentObject(acc.root, true);
  }
}
