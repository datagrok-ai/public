// PLS specific constants

/** Types of analysis using PLS */
export enum PLS_ANALYSIS {
  COMPUTE_COMPONENTS,
  PERFORM_MVA,
}

/** Errors & warnings */
export enum ERROR_MSG {
  NO_DF = 'No dataframe is opened',
  NO_COLS = 'No numeric columns without missing values',
  ONE_COL = 'No columns to be used as features (just one numeric columns without missing values)',
}

/** Widget titles */
export enum TITLE {
  PREDICT = 'Predict',
  USING = 'Using',
  COMPONENTS = 'Components',
  PLS = 'PLS',
  MVA = 'Multivariate analysis (PLS)',
  RUN = 'RUN',
  NAMES = 'Names',
  MODEL = 'Observed vs. Predicted',
  FEATURE = 'Feature',
  REGR_COEFS = 'Regression Coefficients',
  XLOADING = 'x.loading.p',
  LOADINGS = 'Loadings',
  XSCORE = 'x.score.t',
  YSCORE = 'y.score.u',
  SCORES = 'Scores',
}

/** Tooltips */
export enum HINT {
  PREDICT = 'Column with the response variable',
  FEATURES = 'Predictors (features)',
  COMPONENTS = 'Number of PLS components',
  PLS = 'Compute PLS components',
  MVA = 'Perform multivariate analysis',
  NAMES = 'Names of data samples',
}

/** Links to help */
export enum LINK {
  PLS = 'https://datagrok.ai/help/explore/multivariate-analysis/pls#pls-components',
  MVA = 'https://datagrok.ai/help/explore/multivariate-analysis/pls',
  MODEL = 'https://datagrok.ai/help/explore/multivariate-analysis/plots/predicted-vs-reference',
  COEFFS = 'https://datagrok.ai/help/explore/multivariate-analysis/plots/regression-coefficients',
  LOADINGS = 'https://datagrok.ai/help/explore/multivariate-analysis/plots/correlation-loadings',
}

/** Components consts */
export enum COMPONENTS {
  DEFAULT = 3,
  MIN = 1,
}

/**  */
export enum RESULT {
  PREFIX = 'PLS',
  SUFFIX = '(predicted)',
}

export const INT = 'Int';
export const TIMEOUT = 6;
