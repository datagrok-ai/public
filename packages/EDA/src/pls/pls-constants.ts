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
export const RADIUS = [0.49, 0.79, 0.99];
export const LINE_WIDTH = 1;
export const X_COORD = 200;
export const Y_COORD = 200;
export const DELAY = 2000;

export const MAX_ROWS_IN_PREDICTION_TOOLTIP = 20;

export const NUMS_AFTER_COMMA = 3;

/** Curves colors */
export enum COLOR {
  AXIS = '#838383',
  CIRCLE = '#0000FF',
  INVALID = '#EB6767',
  VALID_TEXT = '#4d5261',
  VALID_LINE = '#dbdcdf',
};

/** Intro markdown for demo app */
export const DEMO_INTRO_MD = `# Data
Cars have many correlated features, which makes pattern extraction difficult.

# Model
Predict a car's price from its other features.

# Approach
[**Partial Least Squares (PLS)**](https://en.wikipedia.org/wiki/Partial_least_squares_regression)
regression projects features onto latent factors that best explain the response.

# Essence
PLS finds latent factors that simultaneously:
* capture maximum variance in the features
* maximize correlation with the response`;

/** Description of demo results: wizard components */
export const DEMO_RESULTS = [
  {
    caption: TITLE.MODEL,
    text: 'Closer to the line means better price prediction.',
  },
  {
    caption: TITLE.SCORES,
    text: 'The latent factor values for each sample reflect the similarities and dissimilarities among observations.',
  },
  {
    caption: TITLE.LOADINGS,
    text: 'The impact of each feature on the latent factors: higher loading means stronger influence.',
  },
  {
    caption: TITLE.REGR_COEFS,
    text: 'Parameters of the obtained model: features make different contribution to the prediction.',
  },
  {
    caption: TITLE.EXPL_VAR,
    text: 'How well the latent components fit source data: closer to one means better fit.',
  },
];

