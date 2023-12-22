export enum IMPUTATION_STRATEGY {
  MEAN = 'mean',
  MEDIAN = 'median',
};

/** */
export enum ERROR_MSG {
  NO_DATAFRAME = 'No dataframe is opened',
  NO_MISSING_VALUES = 'No missing values',
  UNSUPPORTED_COLUMN_TYPE = 'Unsupported column type',
  UNSUPPORTED_IMPUTATION_STRATEGY = 'Unsupported imputation strategy',
  KNN_CANNOT_BE_APPLIED = 'KNN imputer cannot be applied: no columns to be used as features',
};

/** */
export const COPY_SUFFIX = '(copy)';

/** */
export enum TITLE {
  KNN_IMPUTER = 'Impute',
  TABLE = 'Table',
  IN_PLACE = 'In-place',
  COLUMN = 'Column',
  FEATURES = 'Features',
  CANCEL = 'CANCEL',
  RUN = 'RUN',  
  WEIGHT = 'metric with the weight',
  NEIGHBORS = 'Neighbors'
};

/** */
export enum LINK {
  KNN_IMPUTER = '/help/explore', // TODO: provide correct link
};

/** */ 
export enum HINT { // TODO: provide tooltips
  TARGET = 'target',
  FEATURES = 'Features',
  IN_PLACE = 'in place',
  DISTANCE = 'dist',
  WEIGHT = 'weight',
  NEIGHBORS = 'neighbors',
}
