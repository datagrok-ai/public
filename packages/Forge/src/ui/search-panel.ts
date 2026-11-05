import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {boolColumn, textColumn} from '../columns';
import {METRIC_DESCRIPTIONS, METRIC_LABELS} from '../metrics/metrics';
import {Leader, SearchCandidate, SearchPlan} from '../training/search';
import {numberColumn, readOnlyGrid} from './data-grid';
import {resultsTitle} from './results-viewers';

interface ShownSearch { plan: SearchPlan; leaders: Leader[] }

/** **Search**: the settings chosen per method, best first, in a grid; the chosen one is marked; a click on a row
 * hands its method and settings to [onRestore]. Hidden without a search. */
export class SearchPanel {
  readonly root: HTMLElement;
  private readonly host = ui.div([]);
  private currentGrid: DG.Grid | null = null;
  private sub: rxjs.Subscription | null = null;

  constructor(private readonly onRestore: (leader: Leader) => void) {
    this.root = ui.divV([resultsTitle('search'), this.host], 'forge-search');
    ui.setDisplay(this.root, false);
  }

  get grid(): DG.Grid | null {
    return this.currentGrid;
  }

  /** Shows the leaders of [search], the chosen one current; hides the panel without one. */
  update(search: ShownSearch | undefined): void {
    this.sub?.unsubscribe();
    this.sub = null;
    this.currentGrid = null;
    ui.empty(this.host);
    const leaders = search?.leaders ?? [];
    ui.setDisplay(this.root, leaders.length > 0);
    if (search === undefined || leaders.length === 0)
      return;
    const {metric, candidates} = search.plan;
    const isHoldOut = search.plan.request.validation.scheme === 'holdout';
    const score = `${METRIC_LABELS[metric]} (${isHoldOut ? 'hold-out' : 'CV'})`;
    const frame = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, '#', leaders.map((_, i) => i + 1)),
      textColumn('Method', leaders.map((l) => l.method)),
      numberColumn(score, leaders.map((l) => l.score)),
      numberColumn('SE', leaders.map((l) => l.se)),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Trials', leaders.map((l) => l.trials)),
      boolColumn('Chosen', leaders.map((l) => l.chosen)),
    ]);
    const headers = {
      '#': 'Rank by the metric.',
      'Method': 'The method and, in the cell\'s tooltip, its chosen settings and the space searched.',
      [score]: `${METRIC_DESCRIPTIONS[metric]} ` + (isHoldOut ? 'On the hold out of the search rows, the same for ' +
        'every trial.' : 'By cross-validation on the search rows, the same folds for every trial.'),
      'SE': 'The standard error of the metric over the folds: a difference within it is noise.',
      'Trials': 'Settings tried for the method, pruned and failed ones included.',
      'Chosen': 'The model trained: the simplest settings within one standard error of the best.',
    };
    const cells: {[column: string]: (row: number) => string} = {
      'Method': (row) => [...Object.entries(leaders[row].hyperparameters).map(([name, value]) => `${name}: ${value}`),
        ...(leaders[row].note === undefined ? [] : [leaders[row].note!]),
        ...SearchPanel.spaceLines(candidates.find((c) => c.engine.name === leaders[row].method))]
        .join('\n'),
      [score]: (row) => `Folds: ${leaders[row].foldScores.map((s) => s.toFixed(3)).join(', ')}`,
      'Chosen': (row) => leaders[row].chosen ? 'Chosen by the one-standard-error rule.' : '',
    };
    const grid = readOnlyGrid(frame, {headers, cells});
    grid.props.showCurrentRowIndicator = true;
    frame.currentRowIdx = Math.max(0, leaders.findIndex((l) => l.chosen));
    this.sub = grid.onCellClick.subscribe((cell) => {
      const row = cell.isTableCell ? cell.tableRowIndex : null;
      if (row !== null && row >= 0)
        this.onRestore(leaders[row]);
    });
    this.currentGrid = grid;
    this.host.append(grid.root);
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
