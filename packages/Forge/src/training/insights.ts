import * as DG from 'datagrok-api/dg';
import {areValuesSimilar, classCounts, fixed, METRIC_DESCRIPTIONS, METRIC_LABELS, metricFindings, MetricId,
  MetricValues, Split, SplitMetrics, trivialModelOf, Underfitting, verdictOf} from '../metrics/metrics';
import {residualsOf, statisticsOf} from '../metrics/statistics';

/** The panel of the results an insight is about. */
type InsightPanel = 'classBalance' | 'statistics' | 'correlation' | 'residuals' | 'confusions' | 'metrics';

/** One line of Insights & Tips: [detail] holds the numbers behind it (one line per `\n`), [items] a list of them;
 * [panel] is the panel of the results that shows what it says. */
export interface Insight { text: string; detail?: string; items?: string[]; panel?: InsightPanel }

const NUMERICAL = DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME;
const IMBALANCE = 0.2;
const MAX_LISTED = 4;
/** The classes listed with their shares; past it, the largest ones and the rest in one line. */
const MAX_SHARES = 10;
const CORRELATION = 0.9;
const MAX_CORRELATED_COLUMNS = 99;
const UNIQUE_SHARE = 0.8;
const LOW_RATIO = 0.6;
const TAIL_FACTOR = 4;

/** The built-in tool's data checks on [columns] and [target]: unbalanced classes (the target's first), highly
 * correlated numerical columns, text columns whose values are almost all different. */
export function dataInsights(columns: DG.Column[], target: DG.Column): Insight[] {
  const all = [...columns, target];
  return [...imbalanceInsights([target, ...columns], target), ...correlationInsights(all, target),
    ...uniqueInsights(all)];
}

/** The built-in tool's checks of a trained model's predictions of [target]: one class found much better than the
 * other for a classifier ([values] are its metrics); one-sided large errors and predicted statistics far from the
 * actual ones for a regressor. */
export function resultInsights(input: {values: MetricValues; target: DG.Column; prediction: DG.Column;
  positiveClass?: string}): Insight[] {
  const {values, target, prediction} = input;
  if (target.matches(NUMERICAL))
    return [...residualInsights(target, prediction), ...statisticInsights(target, prediction)];
  const insights: Insight[] = [];
  const {sensitivity, specificity} = values;
  const classes = target.categories.filter((c) => c !== '');
  const positive = input.positiveClass ?? classes[0];
  const negative = classes.find((c) => c !== positive);
  if (sensitivity === undefined || specificity === undefined || sensitivity + specificity === 0 ||
    negative === undefined)
    return insights;
  const detail = `Recall (Sensitivity) ${sensitivity.toFixed(3)}, Specificity ${specificity.toFixed(3)}.`;
  if (sensitivity / specificity < LOW_RATIO)
    insights.push({text: `The model finds ${negative} much better than ${positive}.`, detail, panel: 'confusions'});
  else if (specificity / sensitivity < LOW_RATIO)
    insights.push({text: `The model finds ${positive} much better than ${negative}.`, detail, panel: 'confusions'});
  return insights;
}

/** The insights about **Performance** by the verdict of [metric] ({@link verdictOf}; the caller resolves it with
 * `verdictMetricOf` or the cards' **Compare by**; undefined: no verdict): an underfitting model, listing every metric
 * barely better than the trivial model's, the verdict's first; a model overfitting by [metric], listing the other
 * metrics far apart; else the metrics far apart, without the verdict word. */
export function metricInsights(metrics: SplitMetrics, metric: MetricId | undefined): Insight[] {
  const verdict = metric === undefined ? undefined : verdictOf(metrics, metric);
  const findings = metricFindings(metrics);
  const mismatches = findings.filter((f) => f.kind === 'mismatch').map((f) => f.id);
  const insights: Insight[] = [];
  const underfitting = verdict?.status === 'underfitting' ? verdict.underfitting : undefined;
  if (underfitting !== undefined) {
    const others = findings.flatMap((f) => f.underfitting !== undefined && f.id !== underfitting.id ?
      [f.underfitting] : []);
    insights.push({text: `The model is underfitting: it is barely better than ${trivialModelOf(underfitting.id)} ` +
      'even on the training rows.', panel: 'metrics', items: [underfitting, ...others].map(underfittingItem)});
  }
  const splitItem = (id: MetricId) => `${METRIC_LABELS[id]}: train ${fixed(metrics.train[id])}, ` +
    `validation ${fixed(metrics.validation[id])}`;
  if (verdict?.status === 'overfitting') {
    const others = mismatches.filter((id) => id !== verdict.metric);
    insights.push({text: `The model may be overfitting by ${splitItem(verdict.metric)}.`,
      panel: 'metrics', ...(others.length === 0 ? {} : {items: others.map(splitItem)})});
  } else if (mismatches.length > 0) {
    insights.push({text: 'Some metrics are far apart on the training and validation rows.', panel: 'metrics',
      items: mismatches.map(splitItem)});
  }
  return insights;
}

/** `<Metric>: train <t>, guessing <b>` (`predicting the mean` for R2). */
function underfittingItem({id, train, baseline}: Underfitting): string {
  return `${METRIC_LABELS[id]}: train ${fixed(train)}, ${trivialModelOf(id)} ${fixed(baseline)}`;
}

/** The insights of a trained or evaluated model: the {@link resultInsights} of [result], then the
 * {@link metricInsights} of the [stored] metrics by [metric]. */
export function evaluationInsights(result: Parameters<typeof resultInsights>[0], stored: SplitMetrics,
  metric: MetricId | undefined): Insight[] {
  return [...resultInsights(result), ...metricInsights(stored, metric)];
}

/** One line reading the model's quality when nothing else is worth a tip: Balanced accuracy for a classifier, R2 for a
 * regressor, on the evaluated rows ([evaluation]), else the validation rows, else the training rows; undefined without
 * either metric. */
export function summaryInsight(stored: Record<Split, MetricValues>, evaluation?: MetricValues): Insight | undefined {
  const metricOf = (values: MetricValues): 'balancedAccuracy' | 'r2' | undefined =>
    values.balancedAccuracy !== undefined ? 'balancedAccuracy' : values.r2 !== undefined ? 'r2' : undefined;
  const [values, where]: [MetricValues, string] = evaluation !== undefined ? [evaluation, 'on the evaluated rows'] :
    metricOf(stored.validation) !== undefined ? [stored.validation, 'on the validation rows'] :
      [stored.train, 'on the training rows'];
  const id = metricOf(values);
  const v = id === undefined ? undefined : values[id];
  if (id === undefined || v === undefined)
    return undefined;
  const reading = id === 'balancedAccuracy' ?
    `on average ${Math.round(100 * v)}% of each class is predicted correctly.` :
    v < 0 ? 'the model predicts worse than the target\'s mean.' :
      `the model explains ${Math.round(100 * v)}% of the target's variation.`;
  return {text: `${METRIC_LABELS[id]} ${v.toFixed(3)} ${where}: ${reading}`, detail: METRIC_DESCRIPTIONS[id],
    panel: 'metrics'};
}

/** A class far from the average count makes its column unbalanced (a feature with a semantic type is not checked);
 * [target]'s classes are shown in **Class balance**. */
function imbalanceInsights(columns: DG.Column[], target: DG.Column): Insight[] {
  const insights: Insight[] = [];
  for (const col of columns.filter((c) => c.isCategorical && (c === target || !c.semType))) {
    const counts = classCounts(col);
    const total = [...counts.values()].reduce((sum, n) => sum + n, 0);
    const reference = total / counts.size;
    if (![...counts.values()].some((n) => n / reference > 1 + IMBALANCE || n / reference < 1 - IMBALANCE))
      continue;
    const detail = shareLines(counts).join('\n');
    insights.push(col === target ? {text: 'The classes are unbalanced.', detail, panel: 'classBalance'} :
      {text: `The classes of ${col.name} are unbalanced.`, detail});
    if (insights.length === MAX_LISTED)
      break;
  }
  return insights;
}

/** `<class>: <share>%` for each class of [counts], one decimal, rounded so the shares add up to exactly 100 (the
 * largest remainders take the tenths left over); past {@link MAX_SHARES} classes, the largest ones and the rest in
 * one line. */
function shareLines(counts: Map<string, number>): string[] {
  let entries = [...counts];
  if (entries.length > MAX_SHARES) {
    const largest = [...entries].sort((a, b) => b[1] - a[1]);
    const kept = new Set(largest.slice(0, MAX_SHARES - 1).map(([name]) => name));
    const rest = largest.slice(MAX_SHARES - 1);
    entries = [...entries.filter(([name]) => kept.has(name)),
      [`other ${rest.length} classes`, rest.reduce((sum, [, n]) => sum + n, 0)]];
  }
  const total = entries.reduce((sum, [, n]) => sum + n, 0);
  const exact = entries.map(([, n]) => n * 1000 / total);
  const tenths = exact.map(Math.floor);
  const left = 1000 - tenths.reduce((sum, t) => sum + t, 0);
  const byRemainder = exact.map((_, i) => i).sort((a, b) => (exact[b] - tenths[b]) - (exact[a] - tenths[a]));
  for (const i of byRemainder.slice(0, left))
    tenths[i]++;
  return entries.map(([name], i) => `${name}: ${(tenths[i] / 10).toFixed(1)}%`);
}

/** One insight for the pairs of numerical columns that move together, each listed with its correlation; it points to
 * **Features correlation** when a pair of features is among them. */
function correlationInsights(columns: DG.Column[], target: DG.Column): Insight[] {
  const numerical = columns.filter((c) => c.matches(NUMERICAL));
  if (numerical.length < 2 || numerical.length > MAX_CORRELATED_COLUMNS)
    return [];
  const pairs: string[] = [];
  let hasFeaturePair = false;
  for (let a = 0; a < numerical.length - 1; a++) {
    const stats = DG.Stats.fromColumn(numerical[a]);
    for (let b = a + 1; b < numerical.length; b++) {
      const r = stats.corr(numerical[b]);
      if (r > CORRELATION) {
        pairs.push(`${numerical[a].name} and ${numerical[b].name}: ${r.toFixed(2)}`);
        hasFeaturePair ||= numerical[a] !== target && numerical[b] !== target;
      }
    }
  }
  return pairs.length === 0 ? [] :
    [{text: 'Features are highly correlated.', items: pairs, panel: hasFeaturePair ? 'correlation' : undefined}];
}

function uniqueInsights(columns: DG.Column[]): Insight[] {
  const insights: Insight[] = [];
  for (const col of columns.filter((c) => c.isCategorical && !c.semType)) {
    const share = col.categories.length / col.length;
    if (share > UNIQUE_SHARE) {
      insights.push({text: `Almost every value of ${col.name} is different; it will not help the model.`,
        detail: `${Math.round(share * 100)}% of its values are unique (above ${UNIQUE_SHARE * 100}%).`});
    }
  }
  return insights;
}

/** The residuals' tails: the top quarter (Q4 - Q3) or the bottom one (Q1 - Q0) more than four times as wide as the
 * middle half (Q3 - Q1). */
function residualInsights(target: DG.Column, prediction: DG.Column): Insight[] {
  const stats = DG.Stats.fromColumn(residualsOf(target, prediction));
  const [q0, q1, q3, q4] = [stats.min, stats.q1, stats.q3, stats.max];
  const middle = `Q3 - Q1 = ${q3.toFixed(3)} - ${q1.toFixed(3)} = ${(q3 - q1).toFixed(3)}`;
  if (q4 - q3 > TAIL_FACTOR * (q3 - q1)) {
    return [{text: 'Large errors are mostly on one side: the model over-predicts.', panel: 'residuals',
      detail: `Q4 - Q3 = ${q4.toFixed(3)} - ${q3.toFixed(3)} = ${(q4 - q3).toFixed(3)}\n${middle}`}];
  }
  if (q1 - q0 > TAIL_FACTOR * (q3 - q1)) {
    return [{text: 'Large errors are mostly on one side: the model under-predicts.', panel: 'residuals',
      detail: `Q1 - Q0 = ${q1.toFixed(3)} - ${q0.toFixed(3)} = ${(q1 - q0).toFixed(3)}\n${middle}`}];
  }
  return [];
}

/** One insight for every statistic of the predictions far from the target's. */
function statisticInsights(target: DG.Column, prediction: DG.Column): Insight[] {
  const far = statisticsOf(target, prediction).filter(({actual, predicted}) => Number.isFinite(actual) &&
    Number.isFinite(predicted) && !areValuesSimilar(actual, predicted));
  if (far.length === 0)
    return [];
  return [{text: 'The predictions\' statistics deviate significantly from the actual ones.', panel: 'statistics',
    items: far.map(({id, actual, predicted}) =>
      `${id}: predicted ${predicted.toFixed(3)}, actual ${actual.toFixed(3)}`)}];
}
