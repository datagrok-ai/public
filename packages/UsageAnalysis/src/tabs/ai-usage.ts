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

type TurnText = {prompt: string | null, response: string | null, toolCalls: string | null};

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
    const body = ui.box();
    const failure = ui.div();
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
      df.then(() => {
        failure.textContent = '';
        ui.setDisplay(body, true);
      }).catch((e) => {
        ui.setDisplay(body, false);
        failure.textContent = '';
        failure.append(ui.info([ui.divText('Check the date and groups in the toolbox, then apply them again.'),
          ui.divText(`${e instanceof Error ? e.message : e}`, 'ua-ai-text-missing')], 'AI usage could not be loaded'));
      });
    });

    const viewer = (name: string, createViewer: (t: DG.DataFrame) => DG.Viewer, getDataFrame = () => df) =>
      new UaFilterableQueryViewer({filterSubscription: this.uaToolbox.filterStream, name, getDataFrame, createViewer});

    const perDay = viewer('Usage over time', (t) => DG.Viewer.lineChart(t, {
      xColumnName: 'period',
      yColumnNames: [AiUsageView.measure(t)],
      yAggrTypes: ['sum'],
      splitColumnNames: ['model'],
      chartTypes: ['Stacked Bar Chart'],
      legendPosition: 'Top',
      title: 'Usage over time',
    }));
    const byUser = viewer('Usage by user', (t) => AiUsageView.barChart(t, 'user', 'Usage by user', AiUsageView.measure(t), {
      onClick: DG.RowGroupAction.Filter,
      rowSource: DG.RowSet.All,
    }));
    const byGroup = viewer('Usage by group', (t) => AiUsageView.barChart(t, 'group', 'Usage by group', AiUsageView.measure(t)),
      () => groups);
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
    body.append(ui.splitV([
      ui.box(cards, {style: {maxHeight: '96px'}}),
      ui.box(ui.splitH([perDay.root, breakdown.root]), {style: {maxHeight: '38%'}}),
      ui.splitH([ui.box(filters.root, {style: {maxWidth: '190px'}}), requests.root]),
    ]));
    this.root.append(failure, body);
  }

  /** Cost when most requests are priced; total tokens otherwise, as for models without a known price. */
  static measure(t: DG.DataFrame): 'cost_usd' | 'total_tokens' {
    const cost = t.col('cost_usd');
    return cost != null && cost.stats.valueCount > t.rowCount / 2 ? 'cost_usd' : 'total_tokens';
  }

  static barChart(t: DG.DataFrame, splitColumnName: string, title: string, valueColumnName: string,
    options: object = {}): DG.Viewer {
    return DG.Viewer.barChart(t, {
      splitColumnName,
      valueColumnName,
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
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'total_tokens', rows.map((i) => col(t, 'total_tokens').get(i))),
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
    const chatId: string | null = col(t, 'chat_id').get(i);
    if (chatId != null) {
      const turns = queries.aiChatContent(chatId).then((c) => {
        const texts = new Map<string, TurnText>();
        for (let j = 0; j < c.rowCount; j++) {
          texts.set(c.getCol('turn_id').get(j), {prompt: c.getCol('prompt').get(j), response: c.getCol('response').get(j),
            toolCalls: c.getCol('tool_calls').get(j)});
        }
        return texts;
      });
      const turnId: string | null = col(t, 'turn_id').get(i);
      const pane = (render: (texts: Map<string, TurnText>) => HTMLElement) => ui.wait(async () => {
        try {
          return render(await turns);
        }
        catch (e) {
          return ui.info(`${e instanceof Error ? e.message : e}`, 'Could not be loaded');
        }
      });
      const text = (value: string | null | undefined) =>
        value ? ui.divText(value, 'ua-ai-text') : ui.divText('Not logged', 'ua-ai-text-missing');
      const current = (texts: Map<string, TurnText>) => turnId == null ? undefined : texts.get(turnId);
      acc.addPane('Prompt', () => pane((texts) => text(current(texts)?.prompt)), true);
      acc.addPane('Answer', () => pane((texts) => text(current(texts)?.response)), true);
      acc.addPane('Tool calls', () => pane((texts) => {
        const calls = current(texts)?.toolCalls?.split('\n') ?? [];
        return calls.length === 0 ? text(null) : ui.table(calls, (call) => {
          const at = call.indexOf(': ');
          return at < 0 ? [call, ''] : [call.slice(0, at), call.slice(at + 2)];
        }, ['Tool', 'Details']);
      }), true);
      const requests = [...Array(t.rowCount).keys()].filter((j) => col(t, 'chat_id').get(j) === chatId)
        .sort((a, b) => col(t, 'time').get(a) - col(t, 'time').get(b));
      if (requests.length > 1) {
        acc.addPane(`Chat (${requests.length} requests)`, () => pane((texts) => ui.table(requests, (j) => {
          const asked = texts.get(col(t, 'turn_id').get(j))?.prompt ?? '';
          const time = s('time', j).slice(5);
          return [j === i ? ui.divText(time) : ui.link(time, () => t.currentRowIdx = j),
            asked.length > 60 ? `${asked.slice(0, 60)}…` : asked, s('total_tokens', j), s('cost_usd', j)];
        }, ['Time', 'Asked', 'Tokens', 'Cost, $'])), true);
      }
    }
    grok.shell.setCurrentObject(acc.root, true);
  }
}
