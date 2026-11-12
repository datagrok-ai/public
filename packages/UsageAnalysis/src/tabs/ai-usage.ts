import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Subscription} from 'rxjs';

import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaToolbox} from '../ua-toolbox';
import {UaFilter} from '../filter';
import {queries} from '../package-api';
import '../../css/usage_analysis.css';

const COLUMNS = ['time', 'period', 'user', 'uid', 'groups', 'feature', 'model', 'status', 'chat_id', 'turn_id',
  'duration_ms', 'input_tokens', 'output_tokens', 'cache_read_tokens', 'cache_write_tokens', 'tokens', 'cost'] as const;
type AiColumn = typeof COLUMNS[number];
type Measure = 'cost' | 'tokens';
type TurnText = (row: number, name: 'prompt' | 'response' | 'tool_calls') => string;
const GRID_WIDTHS: Partial<Record<AiColumn, number>> = {time: 108, user: 86, feature: 82, model: 116, status: 46,
  tokens: 60, cost: 52};
const BAR_CHART = {valueAggrType: 'sum', stackColumnName: '', orientation: 'horizontal', barSortType: 'by value',
  barSortOrder: 'desc'};
const GROUP_SEPARATOR = '\n';
const col = (t: DG.DataFrame, name: AiColumn): DG.Column => t.getCol(name);
const message = (e: unknown) => `${e instanceof Error ? e.message : e}`;

export class AiUsageView extends UaView {
  private requests?: UaFilterableQueryViewer;
  private df?: Promise<DG.DataFrame>;
  private groups?: Promise<DG.DataFrame>;
  private measure: Measure = 'cost';
  private frameSubs: Subscription[] = [];
  private contents = new Map<string, Promise<TurnText>>();
  private panel?: HTMLElement;
  private cards = ui.div([], {classes: 'ua-cards'});
  private body = ui.box();
  private failure = ui.div();

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
    this.uaToolbox.filterStream.subscribe((filter: UaFilter) => this.load(filter));

    const viewer = (name: string, createViewer: (t: DG.DataFrame) => DG.Viewer, getDataFrame = () => this.df!) =>
      new UaFilterableQueryViewer({filterSubscription: this.uaToolbox.filterStream, name, getDataFrame, createViewer});
    const perDay = viewer('Usage over time', (t) => DG.Viewer.lineChart(t, {xColumnName: 'period',
      yColumnNames: [this.measure], yAggrTypes: ['sum'], splitColumnNames: ['model'], chartTypes: ['Stacked Bar Chart'],
      legendPosition: 'Top', title: 'Usage over time'}));
    const byUser = viewer('Usage by user', (t) => DG.Viewer.barChart(t, {...BAR_CHART, splitColumnName: 'user',
      valueColumnName: this.measure, title: 'Usage by user', onClick: DG.RowGroupAction.Filter,
      rowSource: DG.RowSet.All}));
    const byGroup = viewer('Usage by group', (t) => DG.Viewer.barChart(t, {...BAR_CHART, splitColumnName: 'group',
      valueColumnName: this.measure, title: 'Usage by group, a user counts in each of their groups'}),
      () => this.groups!);
    const filters = viewer('Filters', (t) => DG.Viewer.filters(t,
      {columnNames: ['user', 'groups', 'model', 'feature', 'status']}));
    const requests = viewer('Requests', (t) => {
      const grid = DG.Viewer.grid(t, {showRowHeader: false, showColumnGridlines: false, allowRowSelection: false,
        allowBlockSelection: false, showCurrentCellOutline: false});
      grid.columns.setVisible(Object.keys(GRID_WIDTHS));
      for (const [name, width] of Object.entries(GRID_WIDTHS))
        grid.col(name)!.width = width;
      return grid;
    });

    filters.root.classList.add('ua-filters');
    requests.root.classList.add('ui-panel');
    this.requests = requests;
    this.viewers.push(perDay, byUser, byGroup, filters, requests);
    this.body.append(ui.splitV([
      ui.box(this.cards, {style: {maxHeight: '96px'}}),
      ui.box(ui.splitH([perDay.root, ui.tabControl({'By user': byUser.root, 'By group': byGroup.root}).root]),
        {style: {maxHeight: '38%'}}),
      ui.splitH([ui.box(filters.root, {style: {maxWidth: '190px'}}), requests.root]),
    ]));
    this.root.append(this.failure, this.body);
  }

  private load(filter: UaFilter): void {
    for (const sub of this.frameSubs)
      sub.unsubscribe();
    this.frameSubs = [];
    this.contents.clear();
    ui.empty(this.cards);
    ui.empty(this.failure);
    ui.setDisplay(this.body, true);
    if (this.panel != null && grok.shell.o === this.panel)
      grok.shell.setCurrentObject(null, false);
    const query = this.df = queries.aiUsage(filter.date ?? 'this week', filter.groups).then((t) => {
      if (this.df === query)
        this.prepare(t);
      return t;
    });
    this.groups = query.then((t) => this.df === query ? this.groupsFrame(t) : t);
    const fail = (title: string, hint: string, detail = '') => {
      if (this.df !== query)
        return;
      ui.setDisplay(this.body, false);
      this.failure.append(ui.info([ui.divText(hint), ui.divText(detail, 'ua-ai-text-missing')], title));
    };
    query.then((t) => t.rowCount > 0 ? null :
      fail('No AI requests in this period', 'Try a longer period or other groups in the toolbox.'))
      .catch((e) => fail('AI usage could not be loaded',
        'Check the date and groups in the toolbox, then apply them again.', message(e)));
  }

  private prepare(t: DG.DataFrame): void {
    for (const name of COLUMNS)
      col(t, name);
    col(t, 'groups').setTag(DG.TAGS.MULTI_VALUE_SEPARATOR, GROUP_SEPARATOR);
    col(t, 'time').meta.format = 'yyyy-MM-dd HH:mm';
    col(t, 'cost').meta.format = '0.0000';
    const cost = col(t, 'cost').stats;
    this.measure = cost.sum > 0 && cost.valueCount >= col(t, 'tokens').stats.valueCount ? 'cost' : 'tokens';
    this.showCards(t);
    this.frameSubs.push(t.onFilterChanged.subscribe(() => this.showCards(t)),
      t.onCurrentRowChanged.subscribe(() => this.showRequest(t)));
  }

  private showCards(t: DG.DataFrame): void {
    const stats = (name: AiColumn) => DG.Stats.fromColumn(col(t, name), t.filter);
    const n = (x: number) => DG.format(x, 'compact');
    const tokens = stats('tokens');
    const cost = stats('cost');
    const cards: {[label: string]: string[]} = {
      'Requests': [n(t.filter.trueCount)],
      'Users': [n(stats('uid').uniqueCount)],
      'Input tokens': [n(stats('input_tokens').sum)],
      'Output tokens': [n(stats('output_tokens').sum)],
      'Total tokens': [n(tokens.sum), `${n(stats('cache_read_tokens').sum + stats('cache_write_tokens').sum)} cache`],
      'Estimated cost': cost.valueCount === 0 || cost.sum === 0 ? ['—', 'no prices'] : [`$${n(cost.sum)}`,
        cost.valueCount < tokens.valueCount ? `${Math.floor(cost.valueCount * 100 / tokens.valueCount)}% priced` : ''],
    };
    ui.empty(this.cards);
    for (const [label, lines] of Object.entries(cards)) {
      this.cards.append(ui.div([ui.divText(label),
        ui.wait(async () => ui.div(lines.map((s) => ui.divText(s))))], 'ua-card'));
    }
  }

  private groupsFrame(t: DG.DataFrame): DG.DataFrame {
    const lists = col(t, 'groups').categories.map((list) => list === '' ? [] : list.split(GROUP_SEPARATOR));
    const listOf = col(t, 'groups').getRawData();
    const rows: number[] = [];
    const names: string[] = [];
    for (let i = 0; i < t.rowCount; i++) {
      for (const name of lists[listOf[i]]) {
        rows.push(i);
        names.push(name);
      }
    }
    const values = (name: 'cost' | 'tokens') => {
      const raw = col(t, name).getRawData();
      const nulls = raw instanceof Int32Array ? DG.INT_NULL : DG.FLOAT_NULL;
      const out = Float64Array.from(rows, (i) => raw[i] === nulls ? DG.FLOAT_NULL : raw[i]);
      return DG.Column.fromFloat64Array(name, out);
    };
    const groups = DG.DataFrame.fromColumns([DG.Column.fromStrings('group', names), values('cost'), values('tokens')]);
    this.frameSubs.push(t.onFilterChanged.subscribe(() => {
      const bits = t.filter.getBuffer();
      groups.filter.init((j) => (bits[rows[j] >>> 5] & (1 << (rows[j] & 31))) !== 0);
    }));
    return groups;
  }

  private showRequest(t: DG.DataFrame): void {
    const i = t.currentRowIdx;
    if (i < 0)
      return;
    const s = (name: AiColumn) => col(t, name).getString(i);
    const acc = ui.accordion();
    acc.addPane('Details', () => ui.tableFromMap({'User': s('user'),
      'Groups': s('groups').split(GROUP_SEPARATOR).join(', '), 'Feature': s('feature'), 'Model': s('model'),
      'Status': s('status'), 'Cost, $': s('cost'),
      'Duration, s': col(t, 'duration_ms').isNone(i) ? '' : (col(t, 'duration_ms').get(i) / 1000).toFixed(1),
      'Chat': s('chat_id')}), true);
    acc.addPane('Tokens', () => ui.tableFromMap({'Input': s('input_tokens'), 'Output': s('output_tokens'),
      'Cache read': s('cache_read_tokens'), 'Cache write': s('cache_write_tokens'), 'Total': s('tokens')}), true);
    if (s('chat_id') !== '')
      this.addChatPanes(acc, t, i, s('chat_id'));
    this.panel = acc.root;
    grok.shell.setCurrentObject(acc.root, false);
  }

  private addChatPanes(acc: DG.Accordion, t: DG.DataFrame, i: number, chatId: string): void {
    const content = this.contents.get(chatId) ?? queries.aiChatContent(chatId).then((c): TurnText => {
      const turns = new Map<string, number>(c.getCol('turn_id').toList().map((id, j) => [id, j]));
      return (row, name) => {
        const turnId: string | null = col(t, 'turn_id').get(row);
        const j = turnId == null ? undefined : turns.get(turnId);
        return j === undefined ? '' : c.getCol(name).get(j) ?? '';
      };
    });
    this.contents.set(chatId, content);
    const pane = (name: string, render: (turnText: TurnText) => HTMLElement) =>
      acc.addPane(name, () => ui.wait(async () => {
        try {
          return render(await content);
        }
        catch (e) {
          return ui.divText(`Could not be loaded: ${message(e)}`);
        }
      }), true);
    const text = (value: string) =>
      value ? ui.divText(value, 'ua-ai-text') : ui.divText('Not logged', 'ua-ai-text-missing');
    pane('Prompt', (turnText) => text(turnText(i, 'prompt')));
    pane('Answer', (turnText) => text(turnText(i, 'response')));
    pane('Tool calls', (turnText) => {
      const calls = turnText(i, 'tool_calls').split('\n').filter((call) => call !== '');
      if (calls.length === 0)
        return ui.divText('None logged', 'ua-ai-text-missing');
      // ui.table maps each item to its row, so a call repeated in one turn needs its index as the item.
      return ui.table(calls.map((_, k) => k), (k) => {
        const at = calls[k].indexOf(': ');
        return at < 0 ? [calls[k], ''] : [calls[k].slice(0, at), calls[k].slice(at + 2)];
      }, ['Tool', 'Details']);
    });
    const chatIndex = col(t, 'chat_id').categories.indexOf(chatId);
    const chatOf = col(t, 'chat_id').getRawData();
    const chat: number[] = [];
    for (let j = t.rowCount - 1; j >= 0; j--) {
      if (chatOf[j] === chatIndex && t.filter.get(j))
        chat.push(j);
    }
    if (chat.length > 1) {
      pane(`Chat (${chat.length} requests)`, (turnText) => ui.table(chat, (j) => {
        const asked = turnText(j, 'prompt');
        const time = col(t, 'time').getString(j).slice(5);
        return [j === i ? ui.divText(time) : ui.link(time, () => t.currentRowIdx = j),
          asked.length > 60 ? `${asked.slice(0, 60)}…` : asked, col(t, 'tokens').getString(j),
          col(t, 'cost').getString(j)];
      }, ['Time', 'Asked', 'Tokens', 'Cost, $']));
    }
  }
}
