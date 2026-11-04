import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaToolbox} from '../ua-toolbox';
import {UaFilter} from '../filter';
import {queries} from '../package-api';
import '../../css/usage_analysis.css';
import {Subscription} from 'rxjs';

const CARDS = ['Requests', 'Users', 'Input tokens', 'Output tokens', 'Total tokens', 'Estimated cost'];
const COLUMNS = ['time', 'period', 'user', 'uid', 'groups', 'feature', 'model', 'status', 'chat_id', 'turn_id',
  'input_tokens', 'output_tokens', 'cache_read_tokens', 'cache_write_tokens', 'total_tokens', 'cost_usd',
  'duration_ms'] as const;
type AiColumn = typeof COLUMNS[number];
const GRID_COLUMNS: AiColumn[] = ['time', 'user', 'feature', 'model', 'status', 'total_tokens', 'cost_usd'];
// The narrowest widths that fit the pane without scrolling sideways; a wider grid scales them up.
const GRID_WIDTHS: Partial<Record<AiColumn, number>> = {time: 108, user: 86, feature: 82, model: 116, status: 46,
  total_tokens: 60, cost_usd: 52};
const col = (t: DG.DataFrame, name: AiColumn): DG.Column => t.getCol(name);
const BAR_CHART = {valueAggrType: 'sum', stackColumnName: '', orientation: 'horizontal', barSortType: 'by value',
  barSortOrder: 'desc'};

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
    let measure: 'cost_usd' | 'total_tokens' = 'cost_usd';
    let subs: Subscription[] = [];
    this.uaToolbox.filterStream.subscribe((filter: UaFilter) => {
      for (const sub of subs)
        sub.unsubscribe();
      subs = [];
      const load = df = queries.aiUsage(filter.date ?? 'this week', filter.groups).then((t) => {
        for (const name of COLUMNS)
          col(t, name);
        col(t, 'groups').setTag(DG.TAGS.MULTI_VALUE_SEPARATOR, ', ');
        col(t, 'time').meta.format = 'yyyy-MM-dd HH:mm';
        col(t, 'cost_usd').meta.format = '0.0000';
        if (df === load) {
          // An unpriced request would read as $0 in a cost chart, so cost is the default only when all are priced.
          const priced = t.rowCount > 0 && AiUsageView.unpriced(t, [...t.filter.getSelectedIndexes()]) === 0;
          measure = priced ? 'cost_usd' : 'total_tokens';
          subs.push(t.onFilterChanged.subscribe(() => showCards(Promise.resolve(t))),
            t.onCurrentRowChanged.subscribe(() => AiUsageView.showRequest(t)));
        }
        return t;
      });
      groups = df.then((t) => df === load ? AiUsageView.groupsFrame(t, subs) : t);
      failure.textContent = '';
      ui.setDisplay(body, true);
      showCards(df);
      const fail = (title: string, hint: string, detail = '') => {
        if (df !== load)
          return;
        ui.setDisplay(body, false);
        failure.append(ui.info([ui.divText(hint), ui.divText(detail, 'ua-ai-text-missing')], title));
      };
      df.then((t) => t.rowCount > 0 ? null :
        fail('No AI requests in this period', 'Try a longer period or other groups in the toolbox.'))
        .catch((e) => fail('AI usage could not be loaded',
          'Check the date and groups in the toolbox, then apply them again.', `${e instanceof Error ? e.message : e}`));
    });

    const viewer = (name: string, createViewer: (t: DG.DataFrame) => DG.Viewer, getDataFrame = () => df) =>
      new UaFilterableQueryViewer({filterSubscription: this.uaToolbox.filterStream, name, getDataFrame, createViewer});

    const perDay = viewer('Usage over time', (t) => DG.Viewer.lineChart(t, {
      xColumnName: 'period',
      yColumnNames: [measure],
      yAggrTypes: ['sum'],
      splitColumnNames: ['model'],
      chartTypes: ['Stacked Bar Chart'],
      legendPosition: 'Top',
      title: 'Usage over time',
    }));
    const byUser = viewer('Usage by user', (t) => DG.Viewer.barChart(t, {...BAR_CHART, splitColumnName: 'user',
      valueColumnName: measure, title: 'Usage by user', onClick: DG.RowGroupAction.Filter, rowSource: DG.RowSet.All}));
    const byGroup = viewer('Usage by group', (t) => DG.Viewer.barChart(t, {...BAR_CHART, splitColumnName: 'group',
      valueColumnName: measure, title: 'Usage by group, a user counts in each of their groups'}), () => groups);
    const filters = viewer('Filters', (t) => DG.Viewer.filters(t,
      {columnNames: ['user', 'groups', 'model', 'feature', 'status']}));
    let resize: Subscription | undefined;
    const requests = viewer('Requests', (t) => {
      const grid = DG.Viewer.grid(t, {
        showRowHeader: false,
        showColumnGridlines: false,
        allowRowSelection: false,
        allowBlockSelection: false,
        showCurrentCellOutline: false,
      });
      grid.columns.setOrder(GRID_COLUMNS);
      grid.columns.setVisible(GRID_COLUMNS);
      const minimum = Object.values(GRID_WIDTHS).reduce((s, w) => s + w, 0);
      resize?.unsubscribe();
      resize = ui.onSizeChanged(grid.root).subscribe(() => {
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

  /** Requests that used tokens but have no cost; requests without any token counts have nothing to price. */
  static unpriced(t: DG.DataFrame, rows: number[]): number {
    return rows.filter((i) => !col(t, 'total_tokens').isNone(i) && col(t, 'cost_usd').isNone(i)).length;
  }

  static summary(t: DG.DataFrame): {[label: string]: string[]} {
    const rows = Array.from(t.filter.getSelectedIndexes());
    const sum = (name: AiColumn) => {
      const values = col(t, name);
      return rows.reduce((s, i) => values.isNone(i) ? s : s + values.get(i), 0);
    };
    const used = rows.filter((i) => !col(t, 'total_tokens').isNone(i)).length;
    const priced = used - AiUsageView.unpriced(t, rows);
    const n = (x: number) => DG.format(x, 'compact');
    return {
      'Requests': [n(rows.length)],
      'Users': [n(new Set(rows.map((i) => col(t, 'uid').get(i))).size)],
      'Input tokens': [n(sum('input_tokens'))],
      'Output tokens': [n(sum('output_tokens'))],
      'Total tokens': [n(sum('total_tokens')), `${n(sum('cache_read_tokens') + sum('cache_write_tokens'))} cache`],
      'Estimated cost': priced === 0 ? ['—', 'no prices'] :
        [`$${n(sum('cost_usd'))}`, priced < used ? `${Math.floor(priced * 100 / used)}% priced` : ''],
    };
  }

  static groupsFrame(t: DG.DataFrame, subs: Subscription[]): DG.DataFrame {
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
    subs.push(t.onFilterChanged.subscribe(() => groups.filter.init((j) => t.filter.get(rows[j]))));
    return groups;
  }

  static showRequest(t: DG.DataFrame): void {
    const i = t.currentRowIdx;
    if (i < 0)
      return;
    const s = (name: AiColumn, row: number = i) => col(t, name).getString(row);
    const acc = ui.accordion();
    acc.addPane('Details', () => ui.tableFromMap({'User': s('user'), 'Groups': s('groups'), 'Feature': s('feature'),
      'Model': s('model'), 'Status': s('status'), 'Cost, $': s('cost_usd'),
      'Duration, s': col(t, 'duration_ms').isNone(i) ? '' : (col(t, 'duration_ms').get(i) / 1000).toFixed(1),
      'Chat': s('chat_id')}), true);
    acc.addPane('Tokens', () => ui.tableFromMap({'Input': s('input_tokens'), 'Output': s('output_tokens'),
      'Cache read': s('cache_read_tokens'), 'Cache write': s('cache_write_tokens'), 'Total': s('total_tokens')}), true);
    const chatId: string | null = col(t, 'chat_id').get(i);
    if (chatId != null) {
      const chat = queries.aiChatContent(chatId);
      // A turn's text; an empty string when the chat has no logged content for it.
      const said = (c: DG.DataFrame, row: number, name: 'prompt' | 'response' | 'tool_calls'): string => {
        const j = c.getCol('turn_id').toList().indexOf(col(t, 'turn_id').get(row));
        return j < 0 ? '' : c.getCol(name).get(j) ?? '';
      };
      const pane = (render: (c: DG.DataFrame) => HTMLElement) => ui.wait(async () => {
        try {
          return render(await chat);
        }
        catch (e) {
          return ui.divText(`Could not be loaded: ${e instanceof Error ? e.message : e}`);
        }
      });
      const text = (value: string) =>
        value ? ui.divText(value, 'ua-ai-text') : ui.divText('Not logged', 'ua-ai-text-missing');
      acc.addPane('Prompt', () => pane((c) => text(said(c, i, 'prompt'))), true);
      acc.addPane('Answer', () => pane((c) => text(said(c, i, 'response'))), true);
      acc.addPane('Tool calls', () => pane((c) => {
        const calls = said(c, i, 'tool_calls').split('\n').filter((call) => call !== '');
        // Grokky logs tool calls only when there are some and logging is on, so none logged can mean either.
        if (calls.length === 0)
          return ui.divText('None logged', 'ua-ai-text-missing');
        // ui.table maps each item to its row, so a call repeated in one turn needs its index as the item.
        return ui.table(calls.map((_, k) => k), (k) => {
          const call = calls[k];
          const at = call.indexOf(': ');
          return at < 0 ? [call, ''] : [call.slice(0, at), call.slice(at + 2)];
        }, ['Tool', 'Details']);
      }), true);
      const requests = [...Array(t.rowCount).keys()].filter((j) => col(t, 'chat_id').get(j) === chatId)
        .sort((a, b) => col(t, 'time').get(a) - col(t, 'time').get(b));
      if (requests.length > 1) {
        acc.addPane(`Chat (${requests.length} requests)`, () => pane((c) => ui.table(requests, (j) => {
          const asked = said(c, j, 'prompt');
          const time = s('time', j).slice(5);
          return [j === i ? ui.divText(time) : ui.link(time, () => t.currentRowIdx = j),
            asked.length > 60 ? `${asked.slice(0, 60)}…` : asked, s('total_tokens', j), s('cost_usd', j)];
        }, ['Time', 'Asked', 'Tokens', 'Cost, $'])), true);
      }
    }
    grok.shell.setCurrentObject(acc.root, false);
  }
}
