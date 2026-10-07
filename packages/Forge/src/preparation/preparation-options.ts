import {MISSING_VALUES_MODES, MissingValuesMode} from './missing-values';

export const ONE_HOT = 'one-hot';
export const SKIP_UNIQUE_CATEGORIES = 'skip-unique-categories';
export const IGNORE_MISSING = 'ignore-missing';
export const IMPUTE_MISSING = 'impute-missing';
export const BINARY_CLASSIFICATION = 'binary-classification';

export interface MissingValuesRecord {
  mode: MissingValuesMode;
  neighbors?: number;
  distance?: string;
  skippedRows: number;
}

/** `options` of a model; the keys after `missingValues` are those of the platform's built-in models. */
export interface PreparationOptions {
  preprocessingInfo: string[];
  postprocessingInfo: string[];
  missingValues?: MissingValuesRecord;
  positiveClass?: string;
  negativeClass?: string;
  binaryClassificationThreshold?: number;
  targetType?: string;
  allowNulls?: boolean;
}

export function preparationOptionsOf(value: unknown): PreparationOptions {
  const source: Record<string, unknown> = isRecord(value) ? value : {};
  const options: PreparationOptions = {
    preprocessingInfo: stringsOf(source['preprocessingInfo']),
    postprocessingInfo: stringsOf(source['postprocessingInfo']),
  };
  const stored = source['missingValues'];
  const missingValues: Record<string, unknown> = isRecord(stored) ? stored : {};
  const {mode, neighbors, distance, skippedRows} = missingValues;
  if (isMissingValuesMode(mode)) {
    const record: MissingValuesRecord = {mode, skippedRows: typeof skippedRows === 'number' ? skippedRows : 0};
    if (typeof neighbors === 'number')
      record.neighbors = neighbors;
    if (typeof distance === 'string')
      record.distance = distance;
    options.missingValues = record;
  }
  const {positiveClass, negativeClass, binaryClassificationThreshold, targetType, allowNulls} = source;
  if (typeof positiveClass === 'string')
    options.positiveClass = positiveClass;
  if (typeof negativeClass === 'string')
    options.negativeClass = negativeClass;
  if (typeof binaryClassificationThreshold === 'number')
    options.binaryClassificationThreshold = binaryClassificationThreshold;
  if (typeof targetType === 'string')
    options.targetType = targetType;
  if (typeof allowNulls === 'boolean')
    options.allowNulls = allowNulls;
  return options;
}

export function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}

function isMissingValuesMode(value: unknown): value is MissingValuesMode {
  return MISSING_VALUES_MODES.some((mode) => mode === value);
}

/** The strings of [value] when it is an array; otherwise none. */
export function stringsOf(value: unknown): string[] {
  return Array.isArray(value) ? value.filter((v): v is string => typeof v === 'string') : [];
}
