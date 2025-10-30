/* eslint-disable max-len */
// Specific optimization constants and definitions

export enum METHOD {
   NELDER_MEAD = 'Nelder-Mead',
   GRAD_DESC = 'Gradient descent',
 };

export const methodTooltip = new Map([
  [METHOD.NELDER_MEAD, 'The Nelder-Mead method'],
  [METHOD.GRAD_DESC, 'The gradient descent method'],
]);

export enum LOSS {
   MAD = 'MAD',
   RMSE = 'RMSE',
 };

export const lossTooltip = new Map([
  [LOSS.MAD, 'Maximum Absolute Deviation (MAD) - measures only the single largest deviation, ignoring smaller deviations'],
  [LOSS.RMSE, 'Root Mean Square Error (RMSE) - measures average deviation, penalizing larger deviations more heavily'],
]);

/** Grid elements sizes */
export enum GRID_SIZE {
  LOSS_COL_WIDTH = 60,
  ROW_HEIGHT = 280,
  LOSS_GRAPH_WIDTH = 320,
  GOF_VIEWER_WIDTH = 320,
  RADAR_WIDTH = 280,
  BAR_HEIGHT = 100,
}

/** UI titles */
export enum TITLE {
   ITER = 'Iteration',
   LOSS = 'Loss',
   LOSS_GRAPH = 'Fitting profile',
   TARGET = 'Target',
   FIT = 'Fit',
   METHOD = 'method',
   VALUE = 'Value',
   OBJECTIVE = TARGET,
   SAMPLES = 'initial points',
   ID = 'id',
   OBTAINED = 'Simulation',
   LOSS_LOW = 'loss',
   SIMILARITY = 'similarity',
 };

/** Items for the fitting help */
enum HELP_ITEMS {
  GOAL = '`Goal`',
  FIT = '`' + TITLE.FIT + '`',
  MIN = '`min`',
  MAX = '`max`',
  TARGET = '`' + TITLE.TARGET + '`',
  METHOD = '`' + TITLE.METHOD + '`',
  CONTEXT = '`Context Panel (F4)`',
  SAMPLES = '`' + TITLE.SAMPLES + '`',
  LOSS = '`' + TITLE.LOSS_LOW + '`',
  SIMILARITY = '`' + TITLE.SIMILARITY + '`',
};

export enum TIMEOUT {
  MS_TO_SLEEP = 5,
  STYLE_TIMEOUT = 3000,
  RADAR = 200,
};

/** Fitting UI constants */
export enum FITTING_UI {
  SIMILARITY_MIN = 0,
  SIMILARITY_MAX = 100, // %
  SIMILARITY_DEFAULT = 10, // %
  SAMPLES = 10,
};

export const HELP_LINK = 'https://datagrok.ai/help/compute/function-analysis#parameter-optimization';

/** Starting help markdown */
export const STARTING_HELP = `# Fitting

Use fitting to solve an inverse problem: find input conditions leading to specified output constraints.
It computes inputs minimizing deviation measured by [loss function](https://en.wikipedia.org/wiki/Loss_function).

1. In the ${HELP_ITEMS.FIT} block, use switchers to specify inputs to be found:
   * Set ${HELP_ITEMS.MIN} and ${HELP_ITEMS.MAX} values for each selected item. They define the variation range
   * Set values of all other inputs

2. Set output constraints in the ${HELP_ITEMS.TARGET} block:
   * Use switchers to specify target outputs
   * Set target value for each selected item

3. Choose ${HELP_ITEMS.METHOD}. Press <i class="grok-icon fal fa-cog"></i> to specify its settings.

4. Press the "Run" icon <i class="fas fa-play"></i> on the top panel to perform fitting. You will get a
[grid](https://datagrok.ai/help/visualize/viewers/grid) containing

   * loss function values
   * fitted inputs
   * viewers visualizing the goodness of fit:
   * [line chart](https://datagrok.ai/help/visualize/viewers/line-chart) showing the loss function minimization

5. Open ${HELP_ITEMS.CONTEXT}. You will get the function run corresponding to the selected grid row.

6. Press **?** on the top panel to learn more about
[parameters optimization](${HELP_LINK}).

# Learn more

* [Parameters optimization](${HELP_LINK})
* [Nelder-Mead method](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)`;

export enum INDICES {
  SINGLE_VAL = 0,
  SINGLE_CAT = 1,
  DOUBLE_VAL1 = 0,
  DOUBLE_VAL2 = 1,
  DOUBLE_CAT = 2,
  DIFF_STUDIO_OUTPUT = 0,
  SCALAR_VALS_COL = 1,
  OLIVE = 8,
  LIGHT_OLIVE = 16,
  GREY = 7,
  LIGHT_GREY = 15,
  BLUE = 0,
  ORANGE = 1,
};

export enum NAME {
  CATEGORY = '__Category__',
  SIZE = '__Size__',
  SIMULATION = 'Simulation',
  TARGET = 'Target',
  SCALARS = 'Scalars',
  NULL = 'Null',
};

export enum SIZE {
  LINE_CHART_LINE_WIDTH = 2,
  SIMULATION = 1,
  TARGET = 2,
  MIN_MARKER = 1,
  MAX_MARKER = 7,
  LIGHTER_PERC = 50,
};

/** Info panel content */
export const TARGET_DATAFRAME_INFO = `1. Set a dataframe with function(s) values

2. Select a column with the independent variable in the **argument** field

3. Specify one or more target columns with the dependent variable(s) in the **functions** field`;

/** Loss function line chart options */
export const LOSS_FUNC_CHART_OPTS = {
  showYAxis: true,
  showXAxis: true,
  showXSelector: true,
  showYSelectors: true,
  lineColoringType: 'custom',
  lineColor: 4288600805,
  markerColor: 4288600805,
};

export const MIN_RADAR_COLS_COUNT = 4;

/** Modulus for the linear congruential generator of random numbers */
export const RAND_MOD = 2147483647;

/** Multiplier used in the generator of random numbers */
export const RAND_MULT = 16807;

/** Reproducibility settings */
export type ReproSettings = {
  reproducible: boolean,
  seed: number,
};

export const REPRO_DEFAULT = true;
export const SEED_DEFAULT = 10;

/** Early stopping settings */
export type EarlyStoppingSettings = {
  useEarlyStopping: boolean,
  costFuncThreshold: number,
  stopAfter: number,
};

export const EARLY_STOP_DEFAULT = false;
export const COST_FUNC_THRESH = 1e-5;
export const STOP_AFTER_DEFAULT = 1;
export const STOP_AFTER_MIN = 1;
