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
type Said = (row: number, name: 'prompt' | 'response' | 'tool_calls') => string;
// Without fixed widths the grid is wider than its pane and scrolls sideways.
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
  private contents = new Map<string, Promise<Said>>();
  private panel?: HTMLElement;
  private cards = ui.div([], {classes: 'ua-cards ua-ai-cards'});
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
    this.cards.textContent = '';
    this.failure.textContent = '';
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
    // An unpriced request would read as $0 in a cost chart, so cost is the default only when all are priced.
    this.measure = col(t, 'cost').stats.valueCount === col(t, 'tokens').stats.valueCount ? 'cost' : 'tokens';
    this.showCards(t);
    this.frameSubs.push(t.onFilterChanged.subscribe(() => this.showCards(t)),
      t.onCurrentRowChanged.subscribe(() => this.showRequest(t)));
  }

  private showCards(t: DG.DataFrame): void {
    const stats = (name: AiColumn) => DG.Stats.fromColumn(col(t, name), t.filter);
    const n = (x: number) => DG.format(x, 'compact');
    // Requests without any token counts have nothing to price.
    const used = stats('tokens').valueCount;
    const priced = stats('cost').valueCount;
    const cards: {[label: string]: string[]} = {
      'Requests': [n(t.filter.trueCount)],
      'Users': [n(stats('uid').uniqueCount)],
      'Input tokens': [n(stats('input_tokens').sum)],
      'Output tokens': [n(stats('output_tokens').sum)],
      'Total tokens': [n(stats('tokens').sum),
        `${n(stats('cache_read_tokens').sum + stats('cache_write_tokens').sum)} cache`],
      'Estimated cost': priced === 0 ? ['—', 'no prices'] :
        [`$${n(stats('cost').sum)}`, priced < used ? `${Math.floor(priced * 100 / used)}% priced` : ''],
    };
    this.cards.textContent = '';
    for (const [label, lines] of Object.entries(cards)) {
      this.cards.append(ui.div([ui.divText(label),
        ui.wait(async () => ui.div(lines.map((s) => ui.divText(s))))], 'ua-card'));
    }
  }

  private groupsFrame(t: DG.DataFrame): DG.DataFrame {
    const rows: number[] = [];
    const names: string[] = [];
    const column = col(t, 'groups');
    for (let i = 0; i < t.rowCount; i++) {
      for (const name of column.isNone(i) ? [] : (column.get(i) as string).split(GROUP_SEPARATOR)) {
        rows.push(i);
        names.push(name);
      }
    }
    const groups = DG.DataFrame.fromColumns([DG.Column.fromStrings('group', names),
      ...(['cost', 'tokens'] as const).map((name) => {
        const values = col(t, name);
        return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, rows.map((i) => values.get(i)));
      })]);
    this.frameSubs.push(t.onFilterChanged.subscribe(() => groups.filter.init((j) => t.filter.get(rows[j]))));
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
    // Expanded panes render at once, so the content is kept per chat to query it once per load.
    const content = this.contents.get(chatId) ?? queries.aiChatContent(chatId).then((c): Said => {
      const turns = new Map<string, number>(c.getCol('turn_id').toList().map((id, j) => [id, j]));
      return (row, name) => {
        const turnId: string | null = col(t, 'turn_id').get(row);
        const j = turnId == null ? undefined : turns.get(turnId);
        return j === undefined ? '' : c.getCol(name).get(j) ?? '';
      };
    });
    this.contents.set(chatId, content);
    const pane = (name: string, render: (said: Said) => HTMLElement) =>
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
    pane('Prompt', (said) => text(said(i, 'prompt')));
    pane('Answer', (said) => text(said(i, 'response')));
    pane('Tool calls', (said) => {
      const calls = said(i, 'tool_calls').split('\n').filter((call) => call !== '');
      // Grokky logs tool calls only when there are some and logging is on, so none logged can mean either.
      if (calls.length === 0)
        return ui.divText('None logged', 'ua-ai-text-missing');
      // ui.table maps each item to its row, so a call repeated in one turn needs its index as the item.
      return ui.table(calls.map((_, k) => k), (k) => {
        const at = calls[k].indexOf(': ');
        return at < 0 ? [calls[k], ''] : [calls[k].slice(0, at), calls[k].slice(at + 2)];
      }, ['Tool', 'Details']);
    });
    // Rows come newest first; the chat reads oldest first, and only the requests the filters show.
    const chat = col(t, 'chat_id').toList().flatMap((id, j) => id === chatId && t.filter.get(j) ? [j] : [])
      .reverse();
    if (chat.length > 1) {
      pane(`Chat (${chat.length} requests)`, (said) => ui.table(chat, (j) => {
        const asked = said(j, 'prompt');
        const time = col(t, 'time').getString(j).slice(5);
        return [j === i ? ui.divText(time) : ui.link(time, () => t.currentRowIdx = j),
          asked.length > 60 ? `${asked.slice(0, 60)}…` : asked, col(t, 'tokens').getString(j),
          col(t, 'cost').getString(j)];
      }, ['Time', 'Asked', 'Tokens', 'Cost, $']));
    }
  }
}
