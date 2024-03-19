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
export const RADIUS = [0.49, 0.79, 0.99];
export const LINE_WIDTH = 1;
export const X_COORD = 200;
export const Y_COORD = 200;
export const DELAY = 2000;

/** Curves colors */
export enum COLOR {
  AXIS = '#838383',
  CIRCLE = '#0000FF',
};

/** Intro markdown for demo app */
export const DEMO_INTRO_MD = `# Data
Each car has many features - patterns extraction is complicated.

# Model
Predict car price by its other features.

# Try
Press 'RUN' to perform multivariate analysis using partial least squares
([PLS](https://en.wikipedia.org/wiki/Partial_least_squares_regression)) regression.

# Essence
The method computes the latent factors that capture the maximum variance in the features 
while maximizing correlation with the response variable.`;

/** Results markdown for demo app */
export const DEMO_RESULTS_MD = [
  `# ${TITLE.MODEL}
  Closer to the line means better price prediction.`,
  `# ${TITLE.REGR_COEFS}
  The 'diesel' feature affects the price the most.`,
  `# ${TITLE.LOADINGS}
  The impact of each feature on the latent factors: higher loading means stronger influence.`,
  `# ${TITLE.SCORES}
  Similarities & dissimilarities: alfaromeo and mercedes are different.`,
  `# ${TITLE.EXPL_VAR}
  How well the latent components fit source data: closer to one means better fit.`,
  `# Learn more
  
  * [Multivariate analysis](${LINK.MVA}),
  * [ANOVA](https://datagrok.ai/help/explore/anova)`,
];
