import * as DG from 'datagrok-api/dg';
import {ForgeError} from '../forge-error';
import {BINARY_CLASSIFICATION, IGNORE_MISSING, IMPUTE_MISSING, ONE_HOT, PreparationOptions,
  SKIP_UNIQUE_CATEGORIES} from './preparation-options';

/** The preparation steps chosen for a training, applied after the missing values are handled. */
export interface PreparationSteps { oneHot: boolean; skipUniqueCategories: boolean; predictProbability: boolean;
  cutoff: number }
export interface PreparedFeatures { columns: DG.Column[]; target: DG.Column; options: PreparationOptions }

export const DEFAULT_CUTOFF = 0.5;
const RECORDED_ONLY = [IGNORE_MISSING, IMPUTE_MISSING];
const COLUMN_TYPES: string[] = Object.values(DG.COLUMN_TYPE);

/** Training: applies [steps] to [columns] and [target], in this order: skip unique categories, one-hot, predict
 * probability, and records each step that changed something in a copy of [options]. Only Forge's own columns are
 * created; [columns] and [target] are never changed. */
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
  const classes = steps.predictProbability ? twoClasses(target) : null;
  if (classes === null)
    return {columns: features, target, options: prepared};
  const [positiveClass, negativeClass] = classes;
  return {columns: features, target: probabilityTarget(target, positiveClass), options: {...prepared,
    postprocessingInfo: [...prepared.postprocessingInfo, BINARY_CLASSIFICATION], positiveClass, negativeClass,
    binaryClassificationThreshold: steps.cutoff, targetType: target.type}};
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

export function replayPostprocessing(prediction: DG.Column, options: PreparationOptions): DG.Column {
  let result = prediction;
  for (const id of options.postprocessingInfo) {
    if (id !== BINARY_CLASSIFICATION)
      throw new ForgeError(notReplayable(id));
    result = binaryClasses(result, options);
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

/** A text or yes/no column whose values are all different, such as an id: Skip unique categories leaves it out. */
export function hasUniqueCategories(col: DG.Column): boolean {
  return col.isCategorical && col.categories.length === col.length;
}

function notReplayable(id: string): string {
  return `The model uses the preparation step '${id}', which Forge cannot replay yet.`;
}

/** A float target: 1 for [positiveClass], 0 for the other class. A float, not an int: methods that round the
 * predictions of an int target (XGBoost) would return 0/1 instead of probabilities. */
function probabilityTarget(target: DG.Column, positiveClass: string): DG.Column {
  if (target.type !== DG.COLUMN_TYPE.STRING) {
    return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, target.name, Array.from({length: target.length},
      (_, i) => target.isNone(i) ? null : String(target.get(i)) === positiveClass ? 1 : 0));
  }
  // A text column's raw data are the indexes of its categories; the empty category is the missing value.
  const indexes = target.getRawData();
  const positive = target.categories.indexOf(positiveClass);
  const empty = target.categories.indexOf('');
  const values = Float32Array.from({length: target.length},
    (_, i) => indexes[i] === empty ? DG.FLOAT_NULL : indexes[i] === positive ? 1 : 0);
  return DG.Column.fromFloat32Array(target.name, values);
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
  return classes;
}

function isColumnType(type: string | undefined): type is DG.ColumnType {
  return COLUMN_TYPES.some((t) => t === type);
}
