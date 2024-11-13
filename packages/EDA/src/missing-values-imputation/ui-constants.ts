/** Error & info messages */
export enum ERROR_MSG {
  NO_DATAFRAME = 'No dataframe is opened',
  NO_MISSING_VALUES = 'No missing values',
  ONE_AVAILABLE_FEATURE = 'Not enough of feature columns to apply imputation using the KNN method',
  ONE_FEATURE_SELECTED = 'Imputation cannot be applied to',
  UNSUPPORTED_COLUMN_TYPE = 'Unsupported column type',
  UNSUPPORTED_IMPUTATION_STRATEGY = 'Unsupported imputation strategy',
  KNN_CANNOT_BE_APPLIED = 'KNN imputer cannot be applied: no columns to be used as features',
  KNN_NO_TARGET_COLUMNS = 'KNN imputer cannot be applied: no columns with missing values',
  KNN_NO_FEATURE_COLUMNS = 'KNN imputer cannot be applied: no feature columns',
  KNN_NOT_ENOUGH_OF_ROWS = 'KNN imputer cannot be applied: not enough of rows',
  KNN_IMPOSSIBLE_IMPUTATION = 'Imputation is impossible, no features can be used',
  INCORRECT_NEIGHBORS = 'Incorrect number of neighbors',
  KNN_FAILS = 'KNN IMPUTATION FAILS',
  CORE_ISSUE = 'Core issue',
  FAILED_TO_IMPUTE = 'Failed to impute',
  UNSUPPORTED_FILL_VALUE_TYPE = 'Unsupported fill value type',
  EMPTY_COLUMN = 'Column contains just null values',
  FAILS_TO_PREDICT_IMPUTATION_FAILS = 'Failed to predict imputation fails',
  WRONG_PREDICTIONS = 'wrong evaluation of KNN imputation fails',
};

/** Suffix used for column copy */
export const COPY_SUFFIX = 'copy';

/** UI titles */
export enum TITLE {
  KNN_IMPUTER = 'k-NN Imputation',
  TABLE = 'Table',
  IN_PLACE = 'In-place',
  COLUMNS = 'Impute',
  FEATURES = 'Using',
  CANCEL = 'CANCEL',
  RUN = 'RUN',
  OK = 'OK',
  NEIGHBORS = 'Neighbors',
  DISTANCE = 'Distance',
  FILL = 'Fill',
  MARK = 'Mark',
  SIMPLE_IMPUTER = 'Simple impute',
  SETTINGS = 'Settings',
  KEEP_EMPTY = 'Keep empty',
};

/** Help links */
export const KNN_IMPUTER = '/help/explore/missing-values-imputation#the-k-nn-method';

/** Tooltips */
export enum HINT {
  TARGET = 'Columns with missing values that must be filled',
  FEATURES = "Columns with features to be used for determining the 'nearest' elements in the k-NN method",
  IN_PLACE = 'Defines whether to use in-place imputation or add a new column without missing values',
  METRIC = 'Type of metric between the feature values',
  WEIGHT = 'Weight',
  NEIGHBORS = 'Neighbors count used in the KNN method',
  DISTANCE = 'Type of distance between elements with the specified features',
  METRIC_SETTINGS = 'Show additional options',
  FILL_FAILED_ITEMS = 'Impute missing values using a simple approach: mean, median or most frequent',
  MARK_FAILED_ITEMS = 'Mark missing values cells with a color',
  FILL_VALUE = 'Fill value',
  IMPUTATION_SETTINGS = 'Simple imputation settings',
  KEEP_EMPTY = 'Defines whether to keep empty missing values failed to be imputed OR fill them using simple imputation',
  RUN = 'Run imputation using the k-NN method',
};

export const MAX_INPUT_NAME_LENGTH = 15;
