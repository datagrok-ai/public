import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {copyInto, scatteredColumn} from '../apply/apply-model';
import {Engine} from '../engines/engine';
import {ROC_THRESHOLDS, rocCurveOf} from '../metrics/metrics';
import {residualsOf, STATISTIC_DESCRIPTIONS, statisticsOf} from '../metrics/statistics';
import {releaseFrame, sharedFrame} from '../preparation/shared-frame';
import {numberColumn, readOnlyGrid, textColumn} from './data-grid';
import {MetricsSummary} from './metrics-table';
import {classesOf, PanelKey, ResultKind} from './results-layout';
import '../../css/forge.css';

const NUMERICAL = DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME;
const MAX_CONFUSION_CLASSES = 10;
const MIN_PC_COLUMNS = 3;
const MAX_PC_COLUMNS = 20;

type PanelText = Exclude<PanelKey, 'predictedVsActual' | `visualize:${number}`> | 'regressionVsActual' |
  'classifierVsActual' | 'visualize';

/** Each panel's title (`<name>` is the target's or the method's name) and what it shows, the title's tooltip. */
export const PANEL_TEXTS: {[key in PanelText]: [title: string, tooltip: string]} = {
  regressionVsActual: ['Predicted <name> vs Actual', 'The predicted value of each row against the actual one.'],
  classifierVsActual: ['Predicted <name> vs Actual',
    'The rows by two features, colored by the actual class; a cross marks a wrong prediction.'],
  rocCurve: ['ROC Curve', 'How well the predicted probabilities separate the two classes at every cutoff; the ' +
    'marked point is the current Positive class cutoff. The AUC is the one of Performance.'],
  metrics: ['Performance', 'How well the model predicts, on the training rows and in validation.'],
  confusions: ['Confusions', 'How often each actual class is predicted as each class.'],
  wrongPredictions: ['Wrong predictions', 'The rows whose predicted class differs from the actual one.'],
  insights: ['Insights & Tips', 'Notes on the data and the model, and what to check.'],
  residuals: ['Residuals', 'How far each prediction is from the actual value.'],
  distributions: ['Distributions of predicted and actual <name>',
    'How the predicted values are spread compared with the actual ones.'],
  statistics: ['Statistics', 'Summary numbers of the predictions next to those of the actual values.'],
  features: ['Features', 'Each row as a line across the features, colored by the target.'],
  correlation: ['Features correlation', 'How strongly each pair of features moves together.'],
  classBalance: ['Class balance', 'How many rows each class has.'],
  visualize: ['<name>', 'The method\'s own view of the model.'],
};
const STATISTICS_HEADERS: {[column: string]: string} = {
  'Statistic': 'Summary statistic of the column.',
  'Actual': 'The statistic of the target column.',
  'Predicted': 'The statistic of the model\'s predictions on the same rows.',
};

/** A trained or evaluated model's results on the rows of [table] that [keptRows] holds. */
export interface ResultsData {
  engine: Engine;
  /** The rows [keptRows] indexes: the Train view's table, the evaluated table. */
  table: DG.DataFrame;
  /** What the method saw (the prepared columns), kept rows only. */
  columns: DG.Column[];
  /** The selection's target (text for a classifier), kept rows only. */
  target: DG.Column;
  /** The labels or values of the kept rows: Forge's own column, renamed here. */
  prediction: DG.Column;
  /** Predict probability: the method's scores, Forge's own column. */
  score?: DG.Column;
  /** The table's rows the data holds; null for all of them. */
  keptRows: DG.BitSet | null;
  positiveClass?: string;
  /** Predict probability: the **Positive class cutoff** the scores are cut at. */
  cutoff?: number;
  blob: Uint8Array;
  /** The metrics grid's content. */
  summary: MetricsSummary;
}

/** **Wrong predictions**: the table's columns and the prediction spread over its rows, filtered to the wrong rows. */
export interface Mismatches {
  frame: DG.DataFrame;
  predicted: DG.Column;
  subs: rxjs.Subscription[];
}

/** The viewers of the shown result and the data under them. They stay while the data does (the same table, kind,
 * columns, target and kept rows): a new result of the same data only refills the model's columns. */
export interface Shown {
  table: DG.DataFrame;
  kind: ResultKind;
  /** The target's classes, its empty category aside. */
  classes: string[];
  columns: DG.Column[];
  target: DG.Column;
  keptRows: DG.BitSet | null;
  /** The table's row of each frame row, and the frame's row of each table row (-1: skipped); null for the same rows. */
  kept: Int32Array | null;
  toFrame: Int32Array | null;
  /** The viewers' frame: the columns, the target and Forge's columns below. */
  frame: DG.DataFrame;
  prediction: DG.Column;
  /** The residuals of a regression, whether a label is right for a classifier. */
  marks: DG.Column;
  score?: DG.Column;
  /** Forge's own frames of **Distributions ...** and **Statistics** (regression). */
  distributions?: DG.DataFrame;
  statistics?: DG.DataFrame;
  /** **ROC Curve** (scores only): Forge's own frame of the curve's points and the panel's title text. */
  roc?: {frame: DG.DataFrame; title: HTMLElement};
  wrong?: Mismatches;
  subs: rxjs.Subscription[];
  /** A hovered row is being passed on: the events that causes are not passed back. */
  isLinking: boolean;
  /** The viewer panels by key: every panel but the metrics, the insights and the method's own. */
  panels: Map<PanelKey, HTMLElement>;
  /** Their viewers, which show the insights about them as their descriptions. */
  viewers: Map<PanelKey, DG.Viewer>;
}

/** The viewers of [data]: one frame of the columns the method saw, the target, the prediction and Forge's own
 * columns, and the panels over it. A viewer that fails to build gives the frame and the links back. */
export function buildShown(data: ResultsData, kind: ResultKind): Shown {
  const {table, columns, target, prediction, score, keptRows} = data;
  const kept = keptRows?.getSelectedIndexes() ?? null;
  let toFrame: Int32Array | null = null;
  if (kept !== null) {
    toFrame = new Int32Array(table.rowCount).fill(-1);
    for (let i = 0; i < kept.length; i++)
      toFrame[kept[i]] = i;
  }
  const shown: Shown = {table, kind, classes: classesOf(target), columns, target, keptRows, kept, toFrame,
    prediction, score, subs: [], isLinking: false, panels: new Map(), viewers: new Map(), ...previewFrame(data)};
  linkRows(shown, shown.frame, kept, toFrame, shown.subs);
  try {
    if (kind === 'regression')
      regressionPanels(shown);
    else
      classifierPanels(shown, data.positiveClass);
    featurePanels(shown);
  } catch (e) {
    releaseShown(shown);
    throw e;
  }
  return shown;
}

/** [shown] takes the model's columns of [data] when [data] holds the rows and the values it was built on and a result
 * of the same [kind]; false, with nothing changed, otherwise. */
export function refillShown(shown: Shown, data: ResultsData, kind: ResultKind): boolean {
  return isSameData(shown, data, kind) && refill(shown, data);
}

/** The hover links of [shown] dropped, its frames given back. */
export function releaseShown(shown: Shown): void {
  releaseLinkedFrame(shown);
  if (shown.wrong !== undefined)
    releaseLinkedFrame(shown.wrong);
}

/** The hover links of [frame] dropped, the frame given back. */
export function releaseLinkedFrame({frame, subs}: {frame: DG.DataFrame; subs: rxjs.Subscription[]}): void {
  for (const sub of subs)
    sub.unsubscribe();
  releaseFrame(frame);
}

/** The frame of the columns the method saw, the target, the prediction and Forge's own columns (the residuals of a
 * regression or whether a label is right, the scores), and the residuals or the right-label column. Forge's columns
 * get names the table and the frame do not use. */
function previewFrame({table, columns, target, prediction, score}: ResultsData):
  {frame: DG.DataFrame; marks: DG.Column} {
  const taken = new Set([...table.columns.names(), ...columns.map((c) => c.name), target.name]
    .map((name) => name.toLowerCase()));
  const unused = (name: string): string => {
    let free = name;
    for (let i = 2; taken.has(free.toLowerCase()); i++)
      free = `${name} ${i}`;
    taken.add(free.toLowerCase());
    return free;
  };
  prediction.name = unused(`Predicted ${target.name}`);
  let marks: DG.Column;
  if (target.matches(NUMERICAL)) {
    marks = residualsOf(target, prediction);
    marks.name = unused('residuals');
  } else {
    marks = DG.Column.fromType(DG.COLUMN_TYPE.BOOL, unused('correctness'), target.length)
      .init((i) => !isWrong(target, prediction, i));
    marks.meta.markers.setAll({'true': DG.MARKER_TYPE.CIRCLE_BORDER, 'false': DG.MARKER_TYPE.CROSS_X_BORDER});
  }
  if (score !== undefined)
    score.name = unused('score');
  const frame = sharedFrame([...columns, target, prediction, marks, ...(score === undefined ? [] : [score])]);
  return {frame, marks};
}

/** Hovering a row of [frame] hovers its row of [shown]'s table and back; [kept] and [toFrame] map the rows (null: the
 * same rows). A value is set only when it differs, and the events a set fires are not passed back: a table row the
 * frame lacks would otherwise come back as none. */
function linkRows(shown: Shown, frame: DG.DataFrame, kept: Int32Array | null, toFrame: Int32Array | null,
  subs: rxjs.Subscription[]): void {
  const table = shown.table;
  const sync = (from: DG.DataFrame, to: DG.DataFrame, map: Int32Array | null) => {
    const i = from.mouseOverRowIdx;
    const mapped = i < 0 || map === null ? i : map[i];
    if (shown.isLinking || to.mouseOverRowIdx === mapped)
      return;
    shown.isLinking = true;
    try {
      to.mouseOverRowIdx = mapped;
    } finally {
      shown.isLinking = false;
    }
  };
  subs.push(frame.onMouseOverRowChanged.subscribe(() => sync(frame, table, kept)),
    table.onMouseOverRowChanged.subscribe(() => sync(table, frame, toFrame)));
}

/** **Wrong predictions** of [shown]: built when the model gets a row wrong, refilled while it does; when it gets
 * none, the panel goes and its frame is returned, to be given back once the rows are mounted. A grid that fails to
 * build gives its frame and links back. */
export function updateWrongPredictions(shown: Shown): Mismatches | undefined {
  const {table, prediction, marks, keptRows, kept, toFrame} = shown;
  const old = shown.wrong;
  const wrong = shown.kind === 'regression' ? null : DG.BitSet.create(marks.length, (i) => !marks.get(i));
  if (wrong === null || wrong.trueCount === 0) {
    shown.wrong = undefined;
    shown.panels.delete('wrongPredictions');
    return old;
  }
  const isWrongRow = new Uint8Array(table.rowCount);
  for (const i of wrong.getSelectedIndexes())
    isWrongRow[kept === null ? i : kept[i]] = 1;
  if (old !== undefined) {
    if (toFrame !== null) {
      old.predicted.init((i) => toFrame[i] < 0 ? null : prediction.get(toFrame[i]));
      old.predicted.fireValuesChanged();
    }
    old.frame.filter.init((i) => isWrongRow[i] === 1);
    return undefined;
  }
  // The table's own columns and the prediction spread over its rows, under the name the frame gave it.
  const predicted = keptRows === null ? prediction : scatteredColumn(prediction, keptRows, table.rowCount);
  const mismatches: Mismatches = {frame: sharedFrame([...table.columns.toList(), predicted]), predicted, subs: []};
  try {
    mismatches.frame.filter.init((i) => isWrongRow[i] === 1);
    linkRows(shown, mismatches.frame, null, null, mismatches.subs);
    const grid = readOnlyGrid(mismatches.frame, {headers: {[predicted.name]: 'The model\'s prediction for the row.'}});
    grid.props.showRowHeader = true;
    shown.panels.set('wrongPredictions', resultsPanel('wrongPredictions', grid.root));
  } catch (e) {
    releaseLinkedFrame(mismatches);
    throw e;
  }
  shown.wrong = mismatches;
  return undefined;
}

/** [data] holds the rows and the values [shown] was built on, and a result of the same kind. */
function isSameData(shown: Shown, data: ResultsData, kind: ResultKind): boolean {
  const before = [...shown.columns, shown.target];
  const now = [...data.columns, data.target];
  return shown.table === data.table && shown.kind === kind &&
    (shown.score === undefined) === (data.score === undefined) && isSameRows(shown.keptRows, data.keptRows) &&
    before.length === now.length && now.every((c, i) => isSameColumn(before[i], c));
}

function isSameRows(a: DG.BitSet | null, b: DG.BitSet | null): boolean {
  if (a === null || b === null)
    return a === b;
  return a.length === b.length && a.trueCount === b.trueCount &&
    isSamePrefix(a.getSelectedIndexes(), b.getSelectedIndexes(), a.trueCount);
}

/** The same column, or a copy with the same name, type and values. */
function isSameColumn(a: DG.Column, b: DG.Column): boolean {
  if (a === b)
    return true;
  if (a.name !== b.name || a.type !== b.type || a.length !== b.length)
    return false;
  const words = a.type === DG.COLUMN_TYPE.BOOL ? Math.ceil(a.length / 32) : a.length;
  return isSamePrefix(a.getRawData(), b.getRawData(), words) &&
    (a.type !== DG.COLUMN_TYPE.STRING || isSamePrefix(a.categories, b.categories, a.categories.length));
}

function isSamePrefix(a: ArrayLike<number | string>, b: ArrayLike<number | string>, n: number): boolean {
  if (a.length < n || b.length < n)
    return false;
  for (let i = 0; i < n; i++) {
    if (a[i] !== b[i])
      return false;
  }
  return true;
}

/** The model's columns of [shown] take the values of [data]: the prediction, the scores (unless the same, as after a
 * cutoff change), the marks, the regression frames. False, with nothing changed, when a new column cannot be
 * written into the old one. */
function refill(shown: Shown, {prediction, score, positiveClass}: ResultsData): boolean {
  if (!canCopy(shown.prediction, prediction) ||
    (shown.score !== undefined && score !== undefined && !canCopy(shown.score, score)))
    return false;
  copyValues(shown.prediction, prediction);
  if (shown.score !== undefined && score !== undefined && !isSameColumn(shown.score, score)) {
    copyValues(shown.score, score);
    if (shown.roc !== undefined)
      fillRoc(shown.roc.frame, shown.score, shown, positiveClass);
  }
  const {target, marks, distributions, statistics} = shown;
  if (shown.kind === 'regression') {
    copyValues(marks, residualsOf(target, shown.prediction));
    // A scatter plot recomputes its regression line only when the state of its inputs changes (its columns, its
    // filter's version: getRegressionState, d4 regression_line.dart), not when their values do.
    shown.frame.filter.fireChanged();
  } else {
    marks.init((i) => !isWrong(target, shown.prediction, i));
    marks.fireValuesChanged();
  }
  if (distributions !== undefined) {
    copyValues(distributions.getCol('targets'),
      DG.Column.fromFloat32Array('targets', distributionValues(target, shown.prediction)));
  }
  if (statistics !== undefined) {
    const rows = statisticsOf(target, shown.prediction);
    statistics.getCol('Actual').init((i) => rows[i].actual).fireValuesChanged();
    statistics.getCol('Predicted').init((i) => rows[i].predicted).fireValuesChanged();
  }
  return true;
}

/** [to] can hold the values of [from]: the same type and length. */
function canCopy(to: DG.Column, from: DG.Column): boolean {
  return to.type === from.type && to.length === from.length;
}

function copyValues(to: DG.Column, from: DG.Column): void {
  // `init` rebuilds a text column's categories; a value set row by row would keep the stale ones.
  if (to.type === DG.COLUMN_TYPE.STRING)
    to.init((i) => from.get(i));
  else
    copyInto(to, from, 0);
  to.fireValuesChanged();
}

/** The core's widget title (`d4-dialog-header`) of the panel [key], [name] in place of `<name>`, with what the
 * panel shows as the tooltip of its text, [span]. */
export function resultsTitle(key: PanelText, name = '',
  span = ui.span([PANEL_TEXTS[key][0].replace('<name>', () => name)])): HTMLElement {
  ui.tooltip.bind(span, PANEL_TEXTS[key][1]);
  return ui.div([span], 'd4-dialog-header d4-dialog-title forge-results-title');
}

export function resultsPanel(key: PanelText, content: HTMLElement, name = ''): HTMLElement {
  return ui.divV([resultsTitle(key, name), content], 'forge-results-panel');
}

/** The panel [key] of [viewer] in [shown]. */
function addViewer(shown: Shown, key: PanelKey & PanelText, viewer: DG.Viewer, name = ''): void {
  shown.panels.set(key, resultsPanel(key, viewer.root, name));
  shown.viewers.set(key, viewer);
}

/** **ROC Curve** (scores only; {@link rocPanel}), **Predicted <target> vs Actual** (two numerical columns),
 * **Confusions** (fewer than ten classes), **Class balance** (the target's classes over the kept rows, each its own
 * color). **Wrong predictions** follows the prediction ({@link updateWrongPredictions}). */
function classifierPanels(shown: Shown, positiveClass: string | undefined): void {
  const {frame, panels, columns, target, prediction, marks} = shown;
  rocPanel(shown, positiveClass);
  const numerical = columns.filter((c) => c.matches(NUMERICAL));
  if (numerical.length >= 2) {
    panels.set('predictedVsActual', resultsPanel('classifierVsActual', DG.Viewer.scatterPlot(frame,
      {xColumnName: numerical[0].name, yColumnName: numerical[1].name, colorColumnName: target.name,
        markersColumnName: marks.name}).root, target.name));
  }
  if (shown.classes.length < MAX_CONFUSION_CLASSES) {
    addViewer(shown, 'confusions', DG.Viewer.confusionMatrix(frame,
      {xColumnName: target.name, yColumnName: prediction.name}));
  }
  addViewer(shown, 'classBalance', DG.Viewer.barChart(frame, {splitColumnName: target.name,
    stackColumnName: target.name, valueColumnName: target.name, valueAggrType: DG.AGG.TOTAL_COUNT,
    showStackSelector: false, showCategorySelector: false, showValueSelector: false}));
}

/** **ROC Curve** of [shown]'s scores: a scatter plot over Forge's own frame of the curve's points (`Threshold`,
 * `FPR`, `TPR`, one row per threshold of {@link ROC_THRESHOLDS}) drawn as a line, with the dashed non-informative
 * baseline.
 * Nothing but the scores changes it: a retrain refills the frame ({@link fillRoc}), a cutoff moves its current row
 * ({@link markRoc}). */
function rocPanel(shown: Shown, positiveClass: string | undefined): void {
  if (shown.score === undefined)
    return;
  const frame = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Threshold', [...ROC_THRESHOLDS]),
    DG.Column.fromType(DG.COLUMN_TYPE.FLOAT, 'FPR', ROC_THRESHOLDS.length),
    DG.Column.fromType(DG.COLUMN_TYPE.FLOAT, 'TPR', ROC_THRESHOLDS.length)]);
  fillRoc(frame, shown.score, shown, positiveClass);
  frame.meta.formulaLines.addLine({title: 'Non-informative baseline', formula: `\${TPR} = \${FPR}`, width: 1,
    style: 'dashed', min: 0, max: 1});
  const viewer = DG.Viewer.scatterPlot(frame, {xColumnName: 'FPR', yColumnName: 'TPR', linesOrderColumnName: 'FPR',
    linesWidth: 5, markerType: 'dot', showCurrentPoint: true, showXSelector: false, showYSelector: false,
    showSizeSelector: false, showColorSelector: false});
  const title = ui.span([PANEL_TEXTS.rocCurve[0]]);
  shown.panels.set('rocCurve', ui.divV([resultsTitle('rocCurve', '', title), viewer.root], 'forge-results-panel'));
  shown.viewers.set('rocCurve', viewer);
  shown.roc = {frame, title};
}

/** The `FPR` and `TPR` of [frame] from [score] against [shown]'s target, in place. */
function fillRoc(frame: DG.DataFrame, score: DG.Column, {target, classes}: Shown,
  positiveClass: string | undefined): void {
  const {fpr, tpr} = rocCurveOf(target, score, positiveClass ?? classes[0]);
  frame.getCol('FPR').init((i) => fpr[i]).fireValuesChanged();
  frame.getCol('TPR').init((i) => tpr[i]).fireValuesChanged();
}

/** **ROC Curve**'s title with the AUC-ROC of the shown scores (the value of **Performance**, the trapezoid over every
 * score), and the cutoff's point as the curve's current row; each set only when it changes. */
export function markRoc({roc}: Shown, {summary, cutoff}: ResultsData): void {
  if (roc === undefined)
    return;
  const auc = summary.evaluation?.values.auc ?? summary.metrics.train.auc;
  const text = `${PANEL_TEXTS.rocCurve[0]}${auc === undefined ? '' : ` (AUC = ${auc.toFixed(3)})`}`;
  if (roc.title.textContent !== text)
    roc.title.textContent = text;
  const row = cutoff === undefined ? -1 : Math.round(cutoff * (ROC_THRESHOLDS.length - 1));
  if (roc.frame.currentRowIdx !== row)
    roc.frame.currentRowIdx = row;
}

/** **Predicted <target> vs Actual**, **Residuals**, **Distributions of predicted and actual <target>**,
 * **Statistics**. */
function regressionPanels(shown: Shown): void {
  const {frame, panels, target, prediction, marks} = shown;
  panels.set('predictedVsActual', resultsPanel('regressionVsActual', DG.Viewer.scatterPlot(frame,
    {xColumnName: target.name, yColumnName: prediction.name, showRegressionLine: true}).root, target.name));
  const range = Math.max(Math.abs(target.min), Math.abs(target.max));
  addViewer(shown, 'residuals', DG.Viewer.scatterPlot(frame, {xColumnName: prediction.name,
    yColumnName: marks.name, showRegressionLine: true, yMin: -range, yMax: range}));
  const values = distributionValues(target, prediction);
  shown.distributions = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('targets', values),
    DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'label', values.length)
      .init((i) => i < target.length ? 'TARGET' : 'PREDICTION')]);
  panels.set('distributions', resultsPanel('distributions', DG.Viewer.histogram(shown.distributions,
    {valueColumnName: 'targets', splitColumnName: 'label'}).root, target.name));
  const rows = statisticsOf(target, prediction);
  shown.statistics = DG.DataFrame.fromColumns([textColumn('Statistic', rows.map((r) => r.id)),
    numberColumn('Actual', rows.map((r) => r.actual)), numberColumn('Predicted', rows.map((r) => r.predicted))]);
  addViewer(shown, 'statistics', readOnlyGrid(shown.statistics,
    {headers: STATISTICS_HEADERS, cells: {'Statistic': (row) => STATISTIC_DESCRIPTIONS[rows[row].id]}}));
}

/** **Features** and **Features correlation**: both for 3 to 20 numerical columns. */
function featurePanels(shown: Shown): void {
  const {frame, columns, target, prediction} = shown;
  const numerical = columns.filter((c) => c.matches(NUMERICAL)).map((c) => c.name);
  if (numerical.length < MIN_PC_COLUMNS || numerical.length > MAX_PC_COLUMNS)
    return;
  const extra = target.matches(NUMERICAL) ? [target.name, prediction.name] : [];
  addViewer(shown, 'features', DG.Viewer.pcPlot(frame,
    {columnNames: [...numerical, ...extra], colorColumnName: target.name}));
  addViewer(shown, 'correlation', DG.Viewer.correlationPlot(frame,
    {xColumnNames: numerical, yColumnNames: numerical}));
}

/** The built-in tool's distribution values: the target's values then the prediction's (the `targets` column of
 * **Distributions ...**, labelled `TARGET` and `PREDICTION`). */
function distributionValues(target: DG.Column, prediction: DG.Column): Float32Array {
  const n = target.length;
  const valueOf = (col: DG.Column, i: number) => col.isNone(i) ? DG.FLOAT_NULL : col.getNumber(i);
  return Float32Array.from({length: 2 * n}, (_, i) => i < n ? valueOf(target, i) : valueOf(prediction, i - n));
}

/** Row [i] has an actual and a predicted label, and they differ (compared as text, as the metrics do). */
function isWrong(target: DG.Column, prediction: DG.Column, i: number): boolean {
  return !target.isNone(i) && !prediction.isNone(i) && String(target.get(i)) !== String(prediction.get(i));
}
