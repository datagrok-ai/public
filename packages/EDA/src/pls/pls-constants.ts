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
  EMPTY_DF = 'Dataframe is empty',
}

/** Widget titles */
export enum TITLE {
  PREDICT = 'Predict',
  USING = 'Using',
  COMPONENTS = 'Components',
  PLS = 'PLS',
  MVA = 'Multivariate Analysis (PLS)',
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
  EXPL_VAR = 'Explained Variance',
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
  EXPL_VARS = 'https://datagrok.ai/help/explore/multivariate-analysis/plots/explained-variance',
}

/** Components consts */
export enum COMPONENTS {
  DEFAULT = 3,
  MIN = 1,
}

/** Items used for naming results */
export enum RESULT_NAMES {
  PREFIX = 'PLS',
  SUFFIX = '(predicted)',
  COMP = 'component',
  COMPS = 'components',
}

/** Indeces of wasm-computation output */
export enum WASM_OUTPUT_IDX {
  PREDICTION = 0,
  REGR_COEFFS = 1,
  T_SCORES = 2,
  U_SCORES = 3,
  X_LOADINGS = 4,
  Y_LOADINGS = 5,
}

export const INT = 'Int';
export const TIMEOUT = 6;
