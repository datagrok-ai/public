import * as DG from 'datagrok-api/dg';
import {ForgeError} from '../forge-error';
import {BestThresholdMetric, ClassProbabilities, classesOf, isAtThreshold, probabilityOf} from '../metrics/metrics';
import {BINARY_CLASSIFICATION, CLASS_PROBABILITIES, IGNORE_MISSING, IMPUTE_MISSING, ONE_HOT, PreparationOptions,
  SKIP_UNIQUE_CATEGORIES} from './preparation-options';

/** The preparation steps chosen for a training, applied after the missing values are handled. With
 * [thresholdMetric], the threshold of two classes is the one that maximizes it on the out-of-fold probabilities
 * (`recutResult`), [threshold] when there is none. */
export interface PreparationSteps { oneHot: boolean; skipUniqueCategories: boolean; predictProbability: boolean;
  threshold: number; thresholdMetric?: BestThresholdMetric }
export interface PreparedFeatures { columns: DG.Column[]; target: DG.Column; options: PreparationOptions }

export const DEFAULT_THRESHOLD = 0.5;
export const DEFAULT_THRESHOLD_METRIC: BestThresholdMetric = 'f1';
/** The defaults a training starts with, in the Train view and in `Forge:trainModel`. */
export const DEFAULT_STEPS = {skipUniqueCategories: true, predictProbability: false} as const;
const THRESHOLD_PROBLEM = 'Enter a number from 0 to 1.';

/** Why [value] cannot be the classification threshold; null for a number from 0 to 1. */
export function thresholdProblem(value: number | null): string | null {
  return value !== null && value >= 0 && value <= 1 ? null : THRESHOLD_PROBLEM;
}
/** One-hot encoding is on by default while every text or yes/no feature it would encode has at most this many
 * categories. */
const ONE_HOT_MAX_CATEGORIES = 20;
const RECORDED_ONLY = [IGNORE_MISSING, IMPUTE_MISSING];
const COLUMN_TYPES: string[] = Object.values(DG.COLUMN_TYPE);

/** Training: applies [steps] to [columns], in this order: skip unique categories, one-hot, and records each step that
 * changed something in a copy of [options]; predict probability on a text or yes/no target with two or more classes
 * records `class-probabilities` (with the two classes, the threshold and its metric for two). Only Forge's own columns
 * are created; [columns] and [target] are never changed. */
export function prepareFeatures(columns: DG.Column[], target: DG.Column, steps: PreparationSteps,
  options: PreparationOptions): PreparedFeatures {
  const preprocessingInfo = [...options.preprocessingInfo];
  const prepared: PreparationOptions = {...options, preprocessingInfo};
  let features = columns;
  if (steps.skipUniqueCategories) {
    const kept = withoutUniqueCategories(features);
    if (kept.length < features.length) {
      preprocessingInfo.push(SKIP_UNIQUE_CATEGORIES);
      prepared.skippedColumns = features.filter((c) => !kept.includes(c)).map((c) => c.name);
    }
    features = kept;
  }
  const categorical = features.filter((c) => c.isCategorical);
  if (steps.oneHot && categorical.length > 0) {
    prepared.oneHotCategories = Object.fromEntries(categorical.map((c) => [c.name, c.categories]));
    features = oneHotEncoded(features, prepared.oneHotCategories);
    preprocessingInfo.push(ONE_HOT);
  }
  const classes = steps.predictProbability && target.isCategorical ? classesOf(target) : [];
  if (classes.length < 2)
    return {columns: features, target, options: prepared};
  return {columns: features, target, options: {...prepared,
    postprocessingInfo: [...prepared.postprocessingInfo, CLASS_PROBABILITIES],
    ...(classes.length === 2 ? {positiveClass: classes[0], negativeClass: classes[1],
      binaryClassificationThreshold: steps.threshold,
      ...(steps.thresholdMetric === undefined ? {} : {thresholdMetric: steps.thresholdMetric})} : {})}};
}

/** Application: replays `options.preprocessingInfo` on [columns]; returns [columns] itself when there is nothing to
 * replay. One-hot uses the recorded training categories, or without a record (the built-in tool's models) the
 * applied data's own; skip unique categories drops the recorded `skippedColumns`, or without one the core's rule. */
export function replayPreprocessing(columns: DG.Column[], options: PreparationOptions): DG.Column[] {
  const steps = options.preprocessingInfo.filter((id) => !RECORDED_ONLY.includes(id));
  const unknown = steps.find((id) => id !== ONE_HOT && id !== SKIP_UNIQUE_CATEGORIES);
  if (unknown !== undefined)
    throw new ForgeError(notReplayable(unknown));
  let features = columns;
  for (const id of steps) {
    features = id === ONE_HOT ? oneHotEncoded(features, options.oneHotCategories) :
      withoutUniqueCategories(features, options.skippedColumns);
  }
  return features;
}

/** The labels of a method's [prediction] after `options.postprocessingInfo`: `binary-classification` cuts the raw
 * score, `class-probabilities` takes the labels from the method's [probabilities]. */
export function replayPostprocessing(prediction: DG.Column, options: PreparationOptions,
  probabilities?: ClassProbabilities): DG.Column {
  let result = prediction;
  for (const id of options.postprocessingInfo) {
    if (id !== BINARY_CLASSIFICATION && id !== CLASS_PROBABILITIES)
      throw new ForgeError(notReplayable(id));
    result = id === BINARY_CLASSIFICATION ? binaryClasses(result, options) :
      probabilityLabels(result, options, probabilities);
  }
  return result;
}

/** The two classes of a text or yes/no target, the first one positive; null for any other target. */
export function twoClasses(target: DG.Column): [string, string] | null {
  if (!target.isCategorical)
    return null;
  const classes = target.categories.filter((c) => c !== '');
  return classes.length === 2 ? [classes[0], classes[1]] : null;
}

/** One-hot encoding's default for [features]: on while a feature is text or yes/no and every such feature it would
 * encode (the all-unique ones aside while [skipUniqueCategories]) has at most [ONE_HOT_MAX_CATEGORIES] categories. */
export function defaultOneHot(features: DG.Column[], skipUniqueCategories: boolean): boolean {
  const encoded = features.filter((c) => c.isCategorical && !(skipUniqueCategories && hasUniqueCategories(c)));
  return features.some((c) => c.isCategorical) && encoded.every((c) => c.categories.length <= ONE_HOT_MAX_CATEGORIES);
}

/** A text or yes/no column whose values are all different, such as an id: Skip unique categories leaves it out. */
export function hasUniqueCategories(col: DG.Column): boolean {
  return col.isCategorical && col.categories.length === col.length;
}

function notReplayable(id: string): string {
  return `The model uses the preparation step '${id}', which Forge cannot replay yet.`;
}

/** The labels the class [probabilities] give: two classes with the recorded keys, the positive class where its
 * probability reaches the threshold (in single precision, as the probabilities are stored), else the negative one;
 * otherwise the most probable class (ties: the first). A row without probabilities keeps [prediction]'s label; without
 * [probabilities], [prediction] itself. */
function probabilityLabels(prediction: DG.Column, options: PreparationOptions,
  probabilities?: ClassProbabilities): DG.Column {
  if (probabilities === undefined)
    return prediction;
  const {positiveClass, negativeClass, binaryClassificationThreshold: threshold} = options;
  const positive = positiveClass === undefined ? undefined : probabilityOf(probabilities, positiveClass);
  const cut = positive === undefined || negativeClass === undefined || threshold === undefined ? undefined :
    (i: number) => positive.isNone(i) ? undefined :
      isAtThreshold(positive.getNumber(i), threshold) ? positiveClass : negativeClass;
  const labels = DG.Column.fromType(prediction.type, prediction.name, prediction.length);
  for (let i = 0; i < prediction.length; i++) {
    if (prediction.isNone(i))
      continue;
    const label = cut === undefined ? mostProbableClass(probabilities, i) : cut(i);
    labels.setString(i, label ?? String(prediction.get(i)), false);
  }
  labels.compact();
  return labels;
}

/** The class of the largest probability of row [i] (ties: the first); undefined when the row has none. */
function mostProbableClass(probabilities: ClassProbabilities, i: number): string | undefined {
  let best: string | undefined;
  let max = -Infinity;
  for (const [k, col] of probabilities.columns.entries()) {
    if (!col.isNone(i) && col.getNumber(i) > max) {
      max = col.getNumber(i);
      best = probabilities.classes[k];
    }
  }
  return best;
}

/** Every text or yes/no column replaced by one 0/1 column `<column>=<category>` per category, after the other columns:
 * the [recorded] categories of the column (a value not among them is 0 in each), else its own. */
function oneHotEncoded(columns: DG.Column[], recorded?: {[column: string]: string[]}): DG.Column[] {
  const encoded: DG.Column[] = [];
  for (const col of columns.filter((c) => c.isCategorical)) {
    // A text column's raw data are the indexes of its categories (-1, a category it lacks, matches no row).
    const isText = col.type === DG.COLUMN_TYPE.STRING;
    const codes: ArrayLike<number | string> = isText ? col.getRawData().subarray(0, col.length) :
      Array.from({length: col.length}, (_, i) => col.isNone(i) ? '' : String(col.get(i)));
    for (const category of recorded?.[col.name] ?? col.categories) {
      const code = isText ? col.categories.indexOf(category) : category;
      encoded.push(DG.Column.fromInt32Array(`${col.name}=${category}`,
        Int32Array.from(codes, (value) => value === code ? 1 : 0)));
    }
  }
  return [...columns.filter((c) => !c.isCategorical), ...encoded];
}

/** [columns] without the ones training left out: the [skipped] names when recorded, else those whose values are all
 * different. */
function withoutUniqueCategories(columns: DG.Column[], skipped?: string[]): DG.Column[] {
  return columns.filter((c) => skipped !== undefined ? !skipped.includes(c.name) : !hasUniqueCategories(c));
}

function binaryClasses(prediction: DG.Column, options: PreparationOptions): DG.Column {
  const {positiveClass, negativeClass, binaryClassificationThreshold: threshold, targetType} = options;
  if (positiveClass === undefined || negativeClass === undefined || threshold === undefined ||
    !isColumnType(targetType)) {
    throw new ForgeError(`The model's preparation step '${BINARY_CLASSIFICATION}' has incomplete settings, ` +
      'so it cannot be replayed.');
  }
  const classes = DG.Column.fromType(targetType, prediction.name, prediction.length);
  for (let i = 0; i < prediction.length; i++) {
    if (!prediction.isNone(i))
      classes.setString(i, prediction.getNumber(i) >= threshold ? positiveClass : negativeClass, false);
  }
  classes.compact();
  return classes;
}

function isColumnType(type: string | undefined): type is DG.ColumnType {
  return COLUMN_TYPES.some((t) => t === type);
}
