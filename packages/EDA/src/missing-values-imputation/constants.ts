export enum IMPUTATION_STRATEGY {
  MEAN = 'mean',
  MEDIAN = 'median',
};

export enum ERROR_MSG {
  NO_MISSING_VALUES = 'No missing values',
  UNSUPPORTED_COLUMN_TYPE = 'Unsupported column type',
  UNSUPPORTED_IMPUTATION_STRATEGY = 'Unsupported imputation strategy',
};

export const COPY_SUFFIX = '(copy)';
