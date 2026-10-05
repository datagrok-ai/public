import * as DG from 'datagrok-api/dg';
import {ForgeError} from '../forge-error';
import {BINARY_CLASSIFICATION, IGNORE_MISSING, IMPUTE_MISSING, ONE_HOT, PreparationOptions,
  SKIP_UNIQUE_CATEGORIES} from './preparation-options';
import {releaseFrame, sharedFrame} from './shared-frame';

const RECORDED_ONLY = [IGNORE_MISSING, IMPUTE_MISSING];
const COLUMN_TYPES: string[] = Object.values(DG.COLUMN_TYPE);

/** Replays the model's preprocessing on [features]; returns [features] itself when there is nothing to replay,
 * otherwise a new frame that may share its columns (release it with `releaseFrame`). */
export function replayPreprocessing(features: DG.DataFrame, options: PreparationOptions): DG.DataFrame {
  const steps = options.preprocessingInfo.filter((id) => !RECORDED_ONLY.includes(id));
  const unknown = steps.find((id) => id !== ONE_HOT && id !== SKIP_UNIQUE_CATEGORIES);
  if (unknown !== undefined)
    throw new ForgeError(notReplayable(unknown));
  if (steps.length === 0)
    return features;
  // A frame of the same columns: the steps add and remove columns of this frame only.
  const frame = sharedFrame(features.columns.toList());
  let isDone = false;
  try {
    for (const id of steps) {
      if (id === ONE_HOT)
        oneHotEncode(frame);
      else
        removeUniqueCategories(frame);
    }
    isDone = true;
    return frame;
  } finally {
    if (!isDone)
      releaseFrame(frame);
  }
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

function notReplayable(id: string): string {
  return `The model uses the preparation step '${id}', which Forge cannot replay yet.`;
}

/** Replaces each text or boolean column with one 0/1 column per category, named `<column>=<category>`. */
function oneHotEncode(frame: DG.DataFrame): void {
  for (const col of frame.columns.toList().filter((c) => c.isCategorical)) {
    const texts = Array.from({length: col.length}, (_, i) => col.isNone(i) ? '' : String(col.get(i)));
    frame.columns.remove(col);
    for (const category of col.categories) {
      const values = Int32Array.from(texts, (text) => text === category ? 1 : 0);
      frame.columns.add(DG.Column.fromInt32Array(`${col.name}=${category}`, values));
    }
  }
}

function removeUniqueCategories(frame: DG.DataFrame): void {
  for (const col of frame.columns.toList().filter((c) => c.isCategorical && c.categories.length === c.length))
    frame.columns.remove(col);
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
