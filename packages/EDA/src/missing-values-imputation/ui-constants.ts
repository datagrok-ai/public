/** Error & info messages */
export enum ERROR_MSG {
  NO_DATAFRAME = 'No dataframe is opened',
  NO_MISSING_VALUES = 'No missing values',
  UNSUPPORTED_COLUMN_TYPE = 'Unsupported column type',
  UNSUPPORTED_IMPUTATION_STRATEGY = 'Unsupported imputation strategy',
  KNN_CANNOT_BE_APPLIED = 'KNN imputer cannot be applied: no columns to be used as features',
};

/** Suffix used for column copy */
export const COPY_SUFFIX = '(copy)';

/** UI titles */
export enum TITLE {
  KNN_IMPUTER = 'Impute',
  TABLE = 'Table',
  IN_PLACE = 'In-place',
  COLUMNS = 'Columns',
  FEATURES = 'Features',
  CANCEL = 'CANCEL',
  RUN = 'RUN',  
  WEIGHT = 'metric with the weight',
  NEIGHBORS = 'Neighbors',
  DISTANCE = 'Distance',
};

/** Help links */
export enum LINK {
  KNN_IMPUTER = '/help/explore', // TODO: provide correct link
};

/** Tooltips */ 
export enum HINT {
  TARGET = 'Columns with missing values that must be filled',
  FEATURES = "Columns with features to be used for determining the 'nearest' elements in the KNN method",
  IN_PLACE = 'Defines whether to use in-place imputeation or add a new column without missing values',
  METRIC = 'Type of metric between the feature values',
  WEIGHT = 'Weight of the feature',
  NEIGHBORS = 'Neighbors count used in the KNN method',
  DISTANCE = 'Type of distance between elements with the specified features',
  METRIC_SETTINGS = 'Features metrics settings'
}
