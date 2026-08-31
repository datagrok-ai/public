// PLS specific constants

/** Types of analysis using PLS */
export enum PLS_ANALYSIS {
  COMPUTE_COMPONENTS,
  PERFORM_MVA,
  DEMO,
}

/** Errors & warnings */
export enum ERROR_MSG {
  NO_DF = 'No dataframe is opened',
  NO_COLS = 'No numeric columns without missing values',
  ONE_COL = 'No columns to be used as features (just one numeric columns without missing values)',
  EMPTY_DF = 'Dataframe is empty',
  PREDICT = 'Predictors must not contain a response variable',
  ENOUGH = 'Not enough of features',
  COMP_LIN_PLS = 'Components count must be less than the number of features',
  COMP_QUA_PLS = 'Too large components count for the quadratic PLS regression',
  COMP_ROWS = 'Components count must not exceed the number of rows',
  COMPONENTS = 'Components count must be at least 1',
  INV_INP = 'Invalid inputs',
  NULL_COMPS = 'Components count is not specified',
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
  VIP = 'Variable Importance',
  XLOADING = 'x.loading.p',
  LOADINGS = 'Loadings',
  XSCORE = 'x.score.t',
  YSCORE = 'y.score.u',
  SCORES = 'Scores',
  EXPL_VAR = 'Explained Variance',
  EXPLORE = 'Explore',
  FEATURES = 'Feature names',
  BROWSE = 'Browse',
  ANALYSIS = 'Features Analysis',
  QUADRATIC = 'Quadratic',
  BIAS = 'bias',
}

/** Tooltips */
export enum HINT {
  PREDICT = 'Column with the response variable',
  FEATURES = 'Predictors (features)',
  COMPONENTS = 'Number of PLS components',
  PLS = 'Compute PLS components',
  MVA = 'Perform multivariate analysis',
  NAMES = 'Names of data samples',
  QUADRATIC = 'Specifies whether to include squared terms as additional predictors in the PLS model',
}

/** Links to help */
export enum LINK {
  PLS = '/help/explore/multivariate-analysis#pls-components',
  MVA = '/help/explore/multivariate-analysis',
  MODEL = '/help/explore/multivariate-analysis#observed-vs-predicted',
  COEFFS = '/help/explore/multivariate-analysis#regression-coefficients',
  LOADINGS = '/help/explore/multivariate-analysis#loadings',
  EXPL_VARS = '/help/explore/multivariate-analysis#explained-variance',
  SCORES = '/help/explore/multivariate-analysis#scores',
  VIPS = '/help/explore/multivariate-analysis#variable-importance',
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
  VIP = 6,
}

/** Tag of the source table that stores the fitted PLS model */
export const MVA_MODEL_TAG = 'mvaModel';

/** Function that adds the analysis results: called via FuncCall, so that the table
 * creation script replays it when a data-synced project is opened */
export const MVA_TRANSFORM_FUNC = 'multivariateAnalysisTransform';

/** Function that restores the viewer state when a project or layout is opened */
export const MVA_INIT_FUNC = 'mvaModelInitFunction';

export const INT = 'Int';
export const TIMEOUT = 6;
export const LINE_WIDTH = 1;
export const X_COORD = 200;
export const Y_COORD = 200;
export const DELAY = 2000;

export const MAX_ROWS_IN_PREDICTION_TOOLTIP = 20;

export const NUMS_AFTER_COMMA = 3;

/** Curves colors */
export enum COLOR {
  AXIS = '#838383',
  INVALID = '#EB6767',
  VALID_TEXT = '#4d5261',
  VALID_LINE = '#dbdcdf',
};

/** Hotelling's T-squared confidence ellipses drawn on the PLS scores plot */
export const ELLIPSES: {conf: number, color: string}[] = [
  {conf: 0.95, color: '#FF8C00'},
  {conf: 0.99, color: '#0000FF'},
];

/** Intro markdown for demo app */
export const DEMO_INTRO_MD = `# Data
30 cars described by 15 predictors that correlate with each other: size, weight, engine, mileage.
Multiple regression breaks down on data like this.

# Model
Predict a car's price from its other characteristics.

# Approach
[**Partial least squares (PLS)**](https://en.wikipedia.org/wiki/Partial_least_squares_regression)
regression builds a linear model on **latent factors**: combinations of the predictors that
maximize covariance with the response variable.

# Essence
Maximizing covariance balances two goals:
* summarize the variation of the predictors
* correlate with the response`;

/** Description of demo results: wizard components */
export const DEMO_RESULTS = [
  {
    caption: TITLE.MODEL,
    text: 'Closer to the line means better price prediction.',
  },
  {
    caption: TITLE.SCORES,
    text: 'Every car placed by its first two latent factors. Nearby cars are similar; clusters, trends, ' +
      'and outliers show up as patterns. ' +
      'The orange 95% and blue 99% Hotelling’s T² ellipses mark the confidence limits — ' +
      'cars outside them are outliers.',
  },
  {
    caption: TITLE.LOADINGS,
    text: 'How strongly the two latent factors describe each feature: the ones far from the origin are described ' +
      'well, and the ones grouped together correlate with each other.',
  },
  {
    caption: TITLE.VIP,
    text: 'Feature ranking by overall contribution: bars above 1 mark the predictors that matter most.',
  },
  {
    caption: TITLE.EXPL_VAR,
    text: 'The share of variance the latent factors explain, added up over the components: closer to one means ' +
      'a better fit.',
  },
];

