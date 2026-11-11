import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {boolColumn, textColumn} from '../columns';
import {METRIC_DESCRIPTIONS, METRIC_LABELS, MetricId} from '../metrics/metrics';
import {candidateOutcomes, CandidateRecord, cautiousCaption, cautiousScore, isRankedBefore, Leader, SearchCandidate,
  SearchOutcome, SearchPlan} from '../training/search';
import {numberColumn, readOnlyGrid} from './data-grid';
import {ResultsPane} from './results-pane';
import {PANEL_TEXTS, resultsTitle} from './results-viewers';

/** A search to show; with its [outcome], the methods without a completed trial are listed after the leaders. */
interface ShownSearch { plan: SearchPlan; leaders: Leader[]; outcome?: SearchOutcome }
/** What a click on a Search row or bar sets: the method and its settings. */
export type Restored = Pick<Leader, 'method' | 'hyperparameters'>;

/** The chart's initial share of the panel's width; the grid gets the rest. */
const CHART_WEIGHTS = [2, 1];
const NO_SETTINGS = 'No settings to tune';

/** **Search**: the leaders' cautious scores (`<Metric> − SE`) in a bar chart left of the grid of the settings chosen
 * per method, in the chart's order, in a platform splitter; both over one frame, so the current row and the selection
 * are shared; the chosen one is marked; then the methods without a completed trial, their reason in the Method tooltip;
 * a click on a row or a bar hands its method and settings (the defaults for a method without a leader) to
 * [onRestore]. Hidden without a search of two methods or more. */
export class SearchPanel {
  readonly root: HTMLElement;
  private readonly host = ui.div([], 'forge-panel-fill');
  private currentGrid: DG.Grid | null = null;
  private currentChart: DG.BarChartViewer | null = null;
  private subs: rxjs.Subscription[] = [];
  /** The metric of the last search shown. */
  private metric: MetricId | undefined;

  constructor(private readonly onRestore: (leader: Restored) => void) {
    this.root = ui.divV([resultsTitle('search', '', undefined, () => this.metric === undefined ?
      PANEL_TEXTS.search[1] : `Each method's best settings, ranked by ${cautiousCaption(this.metric)}.`), this.host],
    'forge-results-panel forge-search');
    ui.setDisplay(this.root, false);
  }

  get grid(): DG.Grid | null {
    return this.currentGrid;
  }

  /** The bar chart of the leaders' cautious scores; null without a leader. */
  get chart(): DG.BarChartViewer | null {
    return this.currentChart;
  }

  /** Shows the leaders of [search], the chosen one current, then its methods without a completed trial; hides the
   * panel with fewer than two rows (a Tune, a Find of one method), which Performance and Runs describe. */
  update(search: ShownSearch | undefined): void {
    for (const sub of this.subs)
      sub.unsubscribe();
    this.subs = [];
    this.currentGrid = null;
    this.currentChart = null;
    ui.empty(this.host);
    const direction = search?.plan.direction ?? 1;
    const cautiousOf = (l: Leader) => cautiousScore(l, direction);
    const leaders = [...search?.leaders ?? []].sort((a, b) => isRankedBefore(a, b, direction) ? -1 :
      isRankedBefore(b, a, direction) ? 1 : 0);
    const records = search?.outcome === undefined ? [] : candidateOutcomes(search.plan, search.outcome);
    const others = records.filter((r) => !leaders.some((l) => l.method === r.method));
    const count = leaders.length + others.length;
    ui.setDisplay(this.root, count > 1);
    if (search === undefined || count < 2)
      return;
    const {metric, candidates} = search.plan;
    // The leaders, then the methods without a completed trial with their records.
    const rows: {method: string; leader?: Leader; record?: CandidateRecord}[] = [
      ...leaders.map((leader) => ({method: leader.method, leader})),
      ...others.map((record) => ({method: record.method, record}))];
    const candidateOf = (method: string) => candidates.find((c) => c.engine.name === method);
    const eliminated = (method: string): string => {
      const round = records.find((r) => r.method === method)?.eliminatedRound;
      return round === undefined ? '' : `Eliminated in round ${round + 1}`;
    };
    const isHoldOut = search.plan.request.validation.scheme === 'holdout';
    const score = `${METRIC_LABELS[metric]} (${isHoldOut ? 'hold-out' : 'CV'})`;
    const label = METRIC_LABELS[metric];
    const cautious = cautiousCaption(metric);
    this.metric = metric;
    const highest = direction === 1 ? 'highest' : 'lowest';
    const frame = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, '#', rows.map((_, i) => i + 1)),
      textColumn('Method', rows.map((r) => r.method)),
      numberColumn(score, rows.map((r) => r.leader?.score ?? null)),
      numberColumn('SE', rows.map((r) => r.leader?.se ?? null)),
      numberColumn(cautious, rows.map((r) => r.leader === undefined ? null : cautiousOf(r.leader))),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Trials', rows.map((r) => (r.leader ?? r.record)!.trials)),
      boolColumn('Chosen', rows.map((r) => r.leader?.chosen ?? false)),
    ]);
    const headers = {
      '#': `Rank by ${cautious}.`,
      'Method': 'The method and, in the cell\'s tooltip, its chosen settings and the space searched.',
      [score]: `${METRIC_DESCRIPTIONS[metric]} ` + (isHoldOut ? 'On the hold out of the search rows, the same for ' +
        'every trial.' : 'By cross-validation on the search rows, the same folds for every trial.'),
      'SE': 'The standard error of the metric over the folds: a difference within it is noise.',
      [cautious]: `${label} ${direction === 1 ? 'minus' : 'plus'} its standard error ` +
        `${isHoldOut ? '(0 on a hold out: one score)' : 'over the folds'}: a cautious estimate of ${label}. ` +
        `The method with the ${highest} value is chosen.`,
      'Trials': 'Every trial of the method, pruned and failed ones included.',
      'Chosen': `The method with the ${highest} ${cautious}.`,
    };
    // An equal split's short count: the method's space ran out, else its share of the limit did.
    const shareNote = (method: string): string => {
      const k = candidates.findIndex((c) => c.engine.name === method);
      const record = records.find((r) => r.method === method);
      return record === undefined || record.outcome === 'failed' ? '' :
        search.outcome?.resume?.exhausted[k] ? 'Its space was exhausted' : 'Its share of the limit ran out';
    };
    // Successive halving drops a method without a space after the first round.
    const trialsNote = (method: string) => search.plan.strategy === 'equal' ? shareNote(method) :
      (candidateOf(method)?.space ?? []).length === 0 ? NO_SETTINGS : eliminated(method);
    const cells: {[column: string]: (row: number) => string} = {
      'Method': (row) => {
        const {method, leader, record} = rows[row];
        const round = eliminated(method);
        return [...leader === undefined ? [SearchPanel.reasonOf(record!)] :
          [...Object.entries(leader.hyperparameters).map(([name, value]) => `${name}: ${value}`),
            ...(leader.note === undefined ? [] : [leader.note])],
        ...(round === '' ? [] : [round]), ...SearchPanel.spaceLines(candidateOf(method))].join('\n');
      },
      [score]: (row) => rows[row].leader === undefined ? '' :
        `Folds: ${rows[row].leader!.foldScores.map((s) => s.toFixed(3)).join(', ')}`,
      'Trials': (row) => trialsNote(rows[row].method),
      'Chosen': (row) => rows[row].leader?.chosen ? `Chosen: the ${highest} ${cautious}.` : '',
    };
    const restore = (method: string) => {
      const leader = leaders.find((l) => l.method === method);
      const hyperparameters = leader?.hyperparameters ?? candidateOf(method)?.start;
      if (hyperparameters !== undefined)
        this.onRestore({method, hyperparameters: {...hyperparameters}});
    };
    const grid = readOnlyGrid(frame, {headers, cells});
    grid.props.showCurrentRowIndicator = true;
    frame.currentRowIdx = Math.max(0, leaders.findIndex((l) => l.chosen));
    this.subs.push(grid.onCellClick.subscribe((cell) => {
      const row = cell.isTableCell ? cell.tableRowIndex : null;
      if (row !== null && row >= 0)
        restore(rows[row].method);
    }));
    this.currentGrid = grid;
    if (leaders.length === 0) {
      this.host.append(grid.root);
      return;
    }
    const chart = DG.Viewer.barChart(frame, {splitColumnName: 'Method', valueColumnName: cautious,
      valueAggrType: DG.AGG.SUM, filter: `IsNotEmpty(\${${cautious}})`, barSortType: 'by value',
      barSortOrder: direction === 1 ? 'desc' : 'asc', showValueSelector: false,
      showCategorySelector: false, showStackSelector: false, legendVisibility: 'Never', showTitle: false,
      description: `${label} ${direction === 1 ? 'minus' : 'plus'} its standard error: a cautious estimate`,
      descriptionVisibilityMode: 'Always', showSelectedRows: false, showFilteredRows: false, onClick: 'Select'});
    // A bar click selects its method's row.
    this.subs.push(frame.onSelectionChanged.subscribe(() => {
      const [row] = frame.selection.getSelectedIndexes();
      if (frame.selection.trueCount === 1 && rows[row]?.leader !== undefined)
        restore(rows[row].method);
    }));
    this.currentChart = chart;
    // Each viewer fills its pane, a box, which the splitter takes as it is (it wraps anything else in a box); the
    // viewer's root as the pane would keep the grid's own height.
    this.host.append(ResultsPane.weightedRow([chart.root, grid.root], CHART_WEIGHTS, 'forge-panel-fill').root);
  }

  /** Why a method has no leader: its first error, or the budget that abandoned or never started its trials. */
  private static reasonOf(record: CandidateRecord): string {
    return record.outcome === 'failed' ? `No completed trial: ${record.note}` : record.outcome === 'abandoned' ?
      'Abandoned: the budget ran out' : 'Not started: the budget ran out';
  }

  /** `Searched: <dimension>, ...` and, when some are, `Pinned: <setting>, ...` of [candidate]; nothing without a
   * space. */
  private static spaceLines(candidate: SearchCandidate | undefined): string[] {
    const space = candidate?.space;
    if (space === undefined)
      return [];
    const searched = space.map((d) => d.type === 'categorical' ? `${d.name} (${d.choices!.join(', ')})` :
      `${d.name} ${d.low}–${d.high}`);
    const pinned = candidate?.pinned ?? [];
    return [`Searched: ${searched.length > 0 ? searched.join(', ') : 'nothing'}`,
      ...(pinned.length > 0 ? [`Pinned: ${pinned.join(', ')}`] : [])];
  }
}
