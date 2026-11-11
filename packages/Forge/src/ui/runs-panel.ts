import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {textColumn} from '../columns';
import {Engine, Hyperparameters, visibleValuesOf} from '../engines/engine';
import {TrainingRunStatus} from '../generated/db';
import {METRIC_DESCRIPTIONS, METRIC_IDS, METRIC_LABELS, MetricId, MetricValues} from '../metrics/metrics';
import {MissingValuesSettings} from '../preparation/missing-values';
import {PreparationSteps} from '../preparation/pipeline';
import {ValidationSettings} from '../training/k-fold';
import {SEARCH_CAPTIONS, SearchRecord} from '../training/search';
import {MetricsRecord} from '../training/train-model';
import {numberColumn, readOnlyGrid} from './data-grid';
import {resultsTitle} from './results-viewers';
import {TuneStatus} from './tune-settings';

const CANCELLED_TOOLTIP = 'Stopped from the task bar.';
const HEADERS: {[column: string]: string} = {
  '#': 'The order of the trainings of this view.',
  'Method': 'The method and, in the cell\'s tooltip, its settings.',
  'Target': 'The column the model learned to predict.',
  'Features': 'How many columns the model learned from; hover for their names.',
  'Preparation': 'The preparation steps of the training.',
  'Rows': 'Rows the model learned from.',
  'Status': 'completed, failed or cancelled.',
  'Split seed': 'The seed that shuffled the rows into the folds or the hold out; empty when they were not shuffled.',
};

/** One recorded training of a Train view, with what it was trained on. */
export interface SessionRun {
  /** 1 for the first recorded training of the view's table. */
  number: number;
  status: TrainingRunStatus;
  engine: Engine;
  targetName: string;
  featureNames: string[];
  hyperparameters: Hyperparameters;
  missingValues: MissingValuesSettings;
  steps: PreparationSteps;
  /** The target has exactly two classes. */
  twoClasses: boolean;
  /** The features had missing values. */
  featureGaps: boolean;
  /** The record's training rows; 0 when nothing was recorded. */
  rowCount: number;
  validation: ValidationSettings;
  /** The split seed. */
  seed: number;
  metrics?: MetricsRecord;
  error?: string;
  /** The search that chose the method and settings. */
  search?: SearchRecord;
  /** Default, Customized or Auto-tuned, as the settings read when the training started. */
  settingsStatus: TuneStatus;
}

/** What a run trained on: the inputs a restore sets again. */
export type RunInputs = Pick<SessionRun, 'engine' | 'targetName' | 'featureNames' | 'hyperparameters' |
  'missingValues' | 'steps' | 'twoClasses' | 'validation' | 'seed'>;

/** Whether [a] and [b] train the same model: every input but Predict probability, the threshold and its metric, which
 * only re-cut a result (a method whose probabilities are a setting differs by that setting). */
export function isSameInputs(a: RunInputs, b: RunInputs): boolean {
  const text = (x: RunInputs) => JSON.stringify([x.engine.name, x.targetName, x.featureNames, x.hyperparameters,
    x.missingValues, {...x.steps, predictProbability: undefined, threshold: undefined, thresholdMetric: undefined},
    x.validation, x.seed]);
  return text(a) === text(b);
}

/** **Runs**: the trainings of a Train view in a grid, newest first, the newest current; a click on a row hands its
 * run to [onRestore]. Hidden while there is no run. */
export class RunsPanel {
  readonly root: HTMLElement;
  private readonly host = ui.div([], 'forge-panel-fill');
  private currentGrid: DG.Grid | null = null;
  private sub: rxjs.Subscription | null = null;

  constructor(private readonly onRestore: (run: SessionRun) => void) {
    this.root = ui.divV([resultsTitle('runs'), this.host], 'forge-results-panel forge-runs');
    ui.setDisplay(this.root, false);
  }

  get grid(): DG.Grid | null {
    return this.currentGrid;
  }

  /** Shows [list], newest first, the first one current. */
  update(list: SessionRun[]): void {
    this.sub?.unsubscribe();
    this.sub = null;
    this.currentGrid = null;
    // Kills the replaced grid too (the platform's widget disposal).
    ui.empty(this.host);
    ui.setDisplay(this.root, list.length > 0);
    if (list.length === 0)
      return;
    const runs = [...list];
    const metricIds = METRIC_IDS.filter((id) => runs.some((run) => run.metrics?.validation[id] !== undefined));
    const frame = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, '#', runs.map((run) => run.number)),
      textColumn('Method', runs.map((run) => run.engine.name)),
      textColumn('Target', runs.map((run) => run.targetName)),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Features', runs.map((run) => run.featureNames.length)),
      textColumn('Preparation', runs.map((run) => RunsPanel.preparationItems(run).join(', '))),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Rows', runs.map((run) => run.rowCount)),
      textColumn('Status', runs.map((run) => run.status)),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Split seed',
        runs.map((run) => run.validation.shuffle ? run.seed : null)),
      ...metricIds.map((id) => numberColumn(RunsPanel.metricColumn(id),
        runs.map((run) => run.metrics?.validation[id] ?? null))),
    ]);
    const headers = {...HEADERS};
    const cells: {[column: string]: (row: number) => string} = {
      'Method': (row) => [...RunsPanel.searchLines(runs[row].search),
        ...Object.entries(visibleValuesOf(runs[row].engine, runs[row].hyperparameters))
          .map(([name, value]) => `${name}: ${value}`)].join('\n'),
      'Features': (row) => runs[row].featureNames.join('\n'),
      'Status': (row) => runs[row].status === 'failed' ? runs[row].error ?? '' :
        runs[row].status === 'cancelled' ? CANCELLED_TOOLTIP : '',
    };
    for (const id of metricIds) {
      headers[RunsPanel.metricColumn(id)] = METRIC_DESCRIPTIONS[id];
      cells[RunsPanel.metricColumn(id)] = (row) => RunsPanel.metricTooltip(runs[row].metrics, id);
    }
    const grid = readOnlyGrid(frame, {headers, cells});
    grid.props.showCurrentRowIndicator = true;
    frame.currentRowIdx = 0;
    // A click, the current row too; keys and code that move the current row restore nothing.
    this.sub = grid.onCellClick.subscribe((cell) => {
      const row = cell.isTableCell ? cell.tableRowIndex : null;
      if (row !== null && row >= 0)
        this.onRestore(runs[row]);
    });
    this.currentGrid = grid;
    this.host.append(grid.root);
  }

  /** `keep missing values` (only when the features had gaps) or `impute (<neighbors> neighbors, <distance>)`,
   * `one-hot`, `skip unique categories`, `probability (threshold <t>, best <metric>)` (without the metric for a
   * threshold set by hand, `probability` for more than two classes): the steps the training used, as the Preparation
   * cell lists them; empty for none. */
  static preparationItems(run: SessionRun): string[] {
    const {steps} = run;
    const best = steps.thresholdMetric === undefined || run.status !== 'completed' ? '' :
      `, best ${METRIC_LABELS[steps.thresholdMetric]}`;
    return [
      ...RunsPanel.missingValuesText(run),
      ...(steps.oneHot ? ['one-hot'] : []),
      ...(steps.skipUniqueCategories ? ['skip unique categories'] : []),
      ...(steps.predictProbability ? [run.twoClasses ? `probability (threshold ${steps.threshold}${best})` :
        'probability'] : []),
    ];
  }

  private static missingValuesText({missingValues: {mode, impute}, featureGaps}: SessionRun): string[] {
    if (mode === 'keep')
      return featureGaps ? ['keep missing values'] : [];
    if (mode !== 'impute')
      return [];
    return [impute === undefined ? 'impute' : `impute (${impute.neighbors} neighbors, ${impute.distance})`];
  }

  /** `Found by <mode>: <n> trials` for a run a search chose; empty for another. */
  static searchLines(search: SearchRecord | undefined): string[] {
    return search === undefined ? [] : [`Found by ${SEARCH_CAPTIONS[search.mode]}: ${search.trials} trials`];
  }

  private static metricColumn(id: MetricId): string {
    return `${METRIC_LABELS[id]} (validation)`;
  }

  private static metricTooltip(metrics: MetricsRecord | undefined, id: MetricId): string {
    const text = (values: MetricValues | undefined) => values?.[id]?.toFixed(3) ?? '-';
    return `Train ${text(metrics?.train)}, validation ${text(metrics?.validation)}`;
  }
}
