import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {copyInto, scatteredColumn} from '../apply/apply-model';
import {textColumn} from '../columns';
import {Engine} from '../engines/engine';
import {classesOf, Curves, METRIC_LABELS} from '../metrics/metrics';
import {residualsOf, STATISTIC_DESCRIPTIONS, statisticsOf} from '../metrics/statistics';
import {releaseFrame, sharedFrame} from '../preparation/shared-frame';
import {numberColumn, readOnlyGrid} from './data-grid';
import {MetricsSummary} from './metrics-table';
import {PanelKey, ResultKind} from './results-layout';
import '../../css/forge.css';

const NUMERICAL = DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME;
const MAX_CONFUSION_CLASSES = 10;
const MIN_PC_COLUMNS = 3;
const MAX_PC_COLUMNS = 20;

type PanelText = Exclude<PanelKey, 'predictedVsActual' | `visualize:${number}`> | 'regressionVsActual' |
  'classifierVsActual' | 'visualize' | 'runs';

/** Each panel's title (`<name>` is the target's or the method's name) and what it shows, the title's tooltip. */
export const PANEL_TEXTS: {[key in PanelText]: [title: string, tooltip: string]} = {
  regressionVsActual: ['Predicted <name> vs Actual', 'The predicted value of each row against the actual one.'],
  classifierVsActual: ['Predicted <name> vs Actual',
    'The rows by two features, colored by the actual class; a cross marks a wrong prediction.'],
  rocCurve: ['ROC Curve', 'Recall (sensitivity) against the false positive rate (1 − specificity) at every cutoff of ' +
    'the positive class\'s probability (or score); with Predict probability the circle is the current Positive ' +
    'class cutoff. The AUC is the one of Performance.'],
  prCurve: ['PR Curve', 'Precision against recall at every cutoff of the positive class\'s probability (or score); ' +
    'the dashed line is the positive class\'s share of the rows (a model that guesses); with Predict probability ' +
    'the circle is the current Positive class cutoff. The AP (average precision) is the one of Performance.'],
  perClass: ['Per-class', 'Each class against the rest: precision, recall (sensitivity), F1 and, with class ' +
    'probabilities, AUC-ROC, one line per class in the class\'s color.'],
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
  runs: ['Runs', 'The trainings of this view, newest first; click one to set its inputs again.'],
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
  /** The positive class's ROC and PR curves (two classes with a probability or score), from EDA's report. */
  curves?: Curves;
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
  curves?: Curves;
  /** Forge's own frames of **Distributions ...** and **Statistics** (regression). */
  distributions?: DG.DataFrame;
  statistics?: DG.DataFrame;
  /** **ROC Curve** and **PR Curve** (two classes; none otherwise), in {@link CURVES} order. */
  curvePanels: CurvePanel[];
  wrong?: Mismatches;
  /** **Per-class** (three or more classes): Forge's own frame of one row per class, the shown split's values per
   * axis. */
  perClass?: PerClassPanel;
  subs: rxjs.Subscription[];
  /** A hovered row is being passed on: the events that causes are not passed back. */
  isLinking: boolean;
  /** The viewer panels by key: every panel but the metrics, the insights and the method's own. */
  panels: Map<PanelKey, HTMLElement>;
  /** Their viewers, which show the insights about them as their descriptions. */
  viewers: Map<PanelKey, DG.Viewer>;
}

/** A curve's panel: Forge's own frame of the curve's points, the viewer, the panel's title text and the frame's row
 * marked as the cutoff (-1: none). */
interface CurvePanel { curve: Curve; frame: DG.DataFrame; viewer: DG.Viewer; title: HTMLElement; marked: number }

interface PerClassPanel { frame: DG.DataFrame; labels: string[]; axes: PerClassAxis[] }

/** The viewers of [data]: one frame of the columns the method saw, the target, the prediction and Forge's own
 * columns, and the panels over it. A viewer that fails to build gives the frame and the links back. */
export function buildShown(data: ResultsData, kind: ResultKind): Shown {
  const {table, columns, target, prediction, score, curves, keptRows} = data;
  const kept = keptRows?.getSelectedIndexes() ?? null;
  let toFrame: Int32Array | null = null;
  if (kept !== null) {
    toFrame = new Int32Array(table.rowCount).fill(-1);
    for (let i = 0; i < kept.length; i++)
      toFrame[kept[i]] = i;
  }
  const shown: Shown = {table, kind, classes: classesOf(target), columns, target, keptRows, kept, toFrame,
    prediction, score, curves, curvePanels: [], subs: [], isLinking: false, panels: new Map(), viewers: new Map(),
    ...previewFrame(data)};
  linkRows(shown, shown.frame, kept, toFrame, shown.subs);
  try {
    if (kind === 'regression')
      regressionPanels(shown);
    else
      classifierPanels(shown, data);
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
    (shown.curves === undefined) === (data.curves === undefined) &&
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
 * cutoff change), the curves (unless the same object, as after a cutoff change), the marks, the per-class values, the
 * regression frames. False, with nothing changed, when a new column cannot be written into the old one, **Per-class**
 * comes or goes, or its classes or axes differ. */
function refill(shown: Shown, {prediction, score, curves, cutoff, summary}: ResultsData): boolean {
  const shownScore = shown.score;
  const hasScores = shownScore !== undefined && score !== undefined;
  if (!canCopy(shown.prediction, prediction) || (hasScores && !canCopy(shownScore, score)))
    return false;
  const values = shown.kind === 'multiclass' ? perClassValues(summary) : undefined;
  const perClass = values !== undefined && values.labels.length > 0 && values.axes.length > 0 ? values : undefined;
  if ((perClass === undefined) !== (shown.perClass === undefined))
    return false;
  if (perClass !== undefined && (perClass.labels.join('\n') !== shown.perClass!.labels.join('\n') ||
    perClass.axes.join() !== shown.perClass!.axes.join()))
    return false;
  copyValues(shown.prediction, prediction);
  if (hasScores && !isSameColumn(shownScore, score))
    copyValues(shownScore, score);
  if (curves !== undefined && curves !== shown.curves)
    replaceCurves(shown, curves, cutoff);
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
  if (perClass !== undefined) {
    for (const id of perClass.axes)
      shown.perClass!.frame.getCol(PER_CLASS_NAMES[id]).init((i) => perClass.values.get(id)![i]).fireValuesChanged();
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

/** **ROC Curve** and **PR Curve** (two classes; {@link curvePanels}), **Predicted <target> vs Actual** (two numerical
 * columns), **Confusions** (fewer than ten classes), **Class balance** (the target's classes over the kept rows, each
 * its own color), **Per-class** (three or more classes; {@link perClassPanel}). **Wrong predictions** follows the
 * prediction ({@link updateWrongPredictions}). */
function classifierPanels(shown: Shown, {summary, cutoff}: ResultsData): void {
  const {frame, panels, columns, target, prediction, marks} = shown;
  curvePanels(shown, cutoff);
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
  if (shown.kind === 'multiclass')
    perClassPanel(shown, summary);
}

const FPR = 'FPR (1 − Specificity)';
const RECALL = METRIC_LABELS.sensitivity;
const PRECISION = METRIC_LABELS.precision;
const THRESHOLD = 'Threshold';
const MARK = 'Mark';
const POINT = 'point';
const CUTOFF = 'cutoff';
const MARK_SIZE = 8;
const BASELINE = 'Non-informative baseline';

/** A curve of **ROC Curve** / **PR Curve**: its panel, axes, points and dashed non-informative baseline, which way
 * its thresholds run and the **Performance** value its title shows. */
interface Curve {
  key: 'rocCurve' | 'prCurve';
  x: string;
  y: string;
  /** The thresholds, x and y of the points. */
  points: (curves: Curves) => [Float32Array, Float32Array, Float32Array];
  baseline: (y: Float32Array) => string;
  /** The thresholds decrease: the cutoff's row is searched from the end. */
  fromEnd: boolean;
  value: 'auc' | 'averagePrecision';
  valueName: string;
}

const CURVES: Curve[] = [
  {key: 'rocCurve', x: FPR, y: RECALL, points: ({roc}) => [roc.threshold, roc.fpr, roc.tpr],
    baseline: () => `\${${RECALL}} = \${${FPR}}`, fromEnd: true, value: 'auc', valueName: 'AUC'},
  // The first point is the lowest threshold, where every row is predicted positive: its precision is the positive
  // class's share.
  {key: 'prCurve', x: RECALL, y: PRECISION, points: ({pr}) => [pr.threshold, pr.recall, pr.precision],
    baseline: (precision) => `\${${PRECISION}} = ${precision[0].toFixed(4)}`, fromEnd: false,
    value: 'averagePrecision', valueName: 'AP'},
];

/** The axes and the markers column of [curve]'s scatter plot, set again with a new frame. */
function curveAxes(curve: Curve): {[name: string]: string} {
  return {xColumnName: curve.x, yColumnName: curve.y, linesOrderColumnName: THRESHOLD, markersColumnName: MARK};
}

/** [curve]'s frame of [curves]: `Threshold`, x, y and `Mark` - a one-pixel dot per point, a circle at the row
 * [cutoff] labels ([marked], -1 for none; {@link markCurves} moves it) - with the dashed non-informative baseline. */
function curveFrame(curve: Curve, curves: Curves, cutoff: number | undefined):
  {frame: DG.DataFrame; marked: number} {
  const [threshold, x, y] = curve.points(curves);
  const marked = thresholdRow(threshold, cutoff, curve.fromEnd);
  const mark = textColumn(MARK, Array.from({length: x.length}, (_, i) => i === marked ? CUTOFF : POINT));
  mark.meta.markers.setAll({[POINT]: DG.MARKER_TYPE.DOT, [CUTOFF]: DG.MARKER_TYPE.CIRCLE});
  const frame = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array(THRESHOLD, threshold),
    DG.Column.fromFloat32Array(curve.x, x), DG.Column.fromFloat32Array(curve.y, y), mark]);
  frame.meta.formulaLines.addLine({title: BASELINE, formula: curve.baseline(y), width: 1, style: 'dashed', min: 0,
    max: 1});
  return {frame, marked};
}

function curvePanel(shown: Shown, curve: Curve, curves: Curves, cutoff: number | undefined): CurvePanel {
  const {frame, marked} = curveFrame(curve, curves, cutoff);
  const viewer = DG.Viewer.scatterPlot(frame, {...curveAxes(curve), linesWidth: 3, markerDefaultSize: MARK_SIZE,
    showCurrentPoint: true, legendVisibility: 'Never', showXAxis: true, showYAxis: true, autoLayout: false,
    showXSelector: true, showYSelector: true, showSizeSelector: false, showColorSelector: false});
  const title = ui.span([PANEL_TEXTS[curve.key][0]]);
  shown.panels.set(curve.key, ui.divV([resultsTitle(curve.key, '', title), viewer.root], 'forge-results-panel'));
  shown.viewers.set(curve.key, viewer);
  return {curve, frame, viewer, title, marked};
}

/** **ROC Curve** and **PR Curve** of [shown]'s curves (two classes): scatter plots over Forge's own frames of the
 * curves' points (one per distinct probability, as Crux computes them), drawn as lines in threshold order, [cutoff]
 * circled. A retrain replaces the frames ({@link replaceCurves}), a cutoff moves their marks ({@link markCurves}). */
function curvePanels(shown: Shown, cutoff: number | undefined): void {
  const {curves} = shown;
  if (curves === undefined || shown.kind !== 'binary')
    return;
  shown.curvePanels = CURVES.map((curve) => curvePanel(shown, curve, curves, cutoff));
}

/** The curve viewers of [shown] take new frames of [curves], [cutoff] circled. */
function replaceCurves(shown: Shown, curves: Curves, cutoff: number | undefined): void {
  for (const panel of shown.curvePanels) {
    const {frame, marked} = curveFrame(panel.curve, curves, cutoff);
    panel.viewer.dataFrame = frame;
    panel.frame = frame;
    panel.marked = marked;
    panel.viewer.setOptions(curveAxes(panel.curve));
  }
  shown.curves = curves;
}

/** The index of the threshold of [thresholds] that [cutoff] labels, from the start ([fromEnd] false) or the end; -1
 * for none. */
function thresholdRow(thresholds: ArrayLike<number>, cutoff: number | undefined, fromEnd: boolean): number {
  if (cutoff === undefined)
    return -1;
  const n = thresholds.length;
  for (let k = 0; k < n; k++) {
    const i = fromEnd ? n - 1 - k : k;
    if (thresholds[i] >= cutoff)
      return i;
  }
  return -1;
}

/** **ROC Curve**'s and **PR Curve**'s titles with the AUC-ROC and the average precision of the shown scores (the
 * values of **Performance**), and the cutoff's point as each curve's current row and its `circle` marker (the other
 * points are dots; a row is positive when its score is at least the cutoff: the ROC's thresholds decrease, the PR's
 * increase); each set only when it changes. */
export function markCurves({curves, curvePanels}: Shown, {summary, cutoff}: ResultsData): void {
  if (curves === undefined)
    return;
  for (const panel of curvePanels) {
    const {curve, frame, title} = panel;
    const value = summary.evaluation?.values[curve.value] ?? summary.metrics.train[curve.value];
    const text = `${PANEL_TEXTS[curve.key][0]}${titleValue(curve.valueName, value)}`;
    if (title.textContent !== text)
      title.textContent = text;
    const row = thresholdRow(curve.points(curves)[0], cutoff, curve.fromEnd);
    if (panel.marked !== row) {
      const mark = frame.getCol(MARK);
      if (panel.marked >= 0)
        mark.set(panel.marked, POINT, false);
      if (row >= 0)
        mark.set(row, CUTOFF, false);
      mark.fireValuesChanged();
      panel.marked = row;
    }
    if (frame.currentRowIdx !== row)
      frame.currentRowIdx = row;
  }
}

function titleValue(name: string, value: number | undefined): string {
  return value === undefined ? '' : ` (${name} = ${value.toFixed(3)})`;
}

const CLASS = 'Class';
const PER_CLASS_AXES = ['precision', 'sensitivity', 'f1', 'auc'] as const;
type PerClassAxis = typeof PER_CLASS_AXES[number];
const PER_CLASS_NAMES: Record<PerClassAxis, string> = {precision: 'Precision', sensitivity: 'Recall', f1: 'F1',
  auc: 'AUC-ROC'};
/** The most characters of the class names (plus {@link LEGEND_ITEM_CHARS} per class for its color mark) a legend on
 * top of **Per-class** holds in one line of a third of the pane; a longer one is not shown. */
const MAX_TOP_LEGEND_CHARS = 40;
const LEGEND_ITEM_CHARS = 3;

interface PerClassValues { labels: string[]; axes: PerClassAxis[]; values: Map<PerClassAxis, (number | null)[]> }

/** The classes of [summary] (the evaluation's first) and the {@link PER_CLASS_AXES} values of the shown split - the
 * evaluation's when there is one, else the validation's; an axis whose values are all missing is left out. */
function perClassValues(summary: MetricsSummary): PerClassValues {
  const stored = summary.metrics.classes;
  const evaluated = summary.evaluation?.classes;
  const labels = [...new Set([...evaluated ?? [], ...stored ?? []].map((c) => c.label))];
  const shownOf = (label: string) => evaluated !== undefined ?
    evaluated.find((c) => c.label === label)?.values : stored?.find((c) => c.label === label)?.validation;
  const values: PerClassValues['values'] = new Map();
  for (const id of PER_CLASS_AXES) {
    const axis = labels.map((label) => shownOf(label)?.[id] ?? null);
    if (axis.some((v) => v !== null))
      values.set(id, axis);
  }
  return {labels, axes: [...values.keys()], values};
}

/** The color of each of [labels] in [target]'s rows (its user color coding kept); a label the rows lack gets the
 * categorical palette's color of its index. */
function classColors(target: DG.Column, labels: string[]): {[label: string]: number} {
  const firstRow = new Map<string, number>();
  for (let i = 0; i < target.length; i++) {
    const text = target.isNone(i) ? undefined : String(target.get(i));
    if (text !== undefined && !firstRow.has(text))
      firstRow.set(text, i);
  }
  return Object.fromEntries(labels.map((label, k) => {
    const row = firstRow.get(label);
    return [label, row === undefined ? DG.Color.getCategoricalColor(k) : DG.Color.getRowColor(target, row)];
  }));
}

/** **Per-class**: a parallel coordinates plot of one line per class of [summary] over the {@link PER_CLASS_AXES}
 * values of the shown split - the evaluation's when there is one, else the validation's - in the class's own color
 * (the target's), the legend listing the classes on top when it is short ({@link MAX_TOP_LEGEND_CHARS}), else none;
 * no in-chart filters, no color selector. */
function perClassPanel(shown: Shown, summary: MetricsSummary): void {
  const {labels, axes, values} = perClassValues(summary);
  if (labels.length === 0 || axes.length === 0)
    return;
  const classes = textColumn(CLASS, labels);
  classes.meta.colors.setCategorical(classColors(shown.target, labels));
  const names = axes.map((id) => PER_CLASS_NAMES[id]);
  const frame = DG.DataFrame.fromColumns([classes, ...axes.map((id, k) => numberColumn(names[k], values.get(id)!))]);
  const legendChars = labels.reduce((sum, label) => sum + label.length + LEGEND_ITEM_CHARS, 0);
  addViewer(shown, 'perClass', DG.Viewer.pcPlot(frame, {columnNames: names, colorColumnName: CLASS,
    showCurrentLine: false, showFilters: false, showColorSelector: false, legendPosition: 'Top',
    legendVisibility: legendChars <= MAX_TOP_LEGEND_CHARS ? 'Always' : 'Never'}));
  shown.perClass = {frame, labels, axes};
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
