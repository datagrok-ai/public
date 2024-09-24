/** Diff Studio application UI constants */

/** Hot keys */
export enum HOT_KEY {
  RUN = 'F5',
};

const COMPUTATION_TIME_UNITS = 'sec';

/** Tooltips messages */
export enum HINT {
  HELP = 'Open help in a new tab',
  OPEN = 'Open model',
  SAVE_MODEL = 'Save model',
  SAVE_LOC = 'Save model to local file',
  SAVE_MY = 'Save model to My files',
  LOAD = 'Load model from local file',
  BASIC = 'Basic template',
  ADV = 'Advanced template',
  EXT = 'Extended template',
  CHEM = 'Mass-action kinetics illustration',
  ROB = `Robertson's chemical reaction model`,
  FERM = 'Fermentation process simulation',
  PK = 'Pharmacokinetic model',
  PKPD = 'Pharmacokinetic-pharmacodynamic model',
  ACID = 'Gluconic acid production model',
  NIM = 'Nimotuzumab disposition model',
  BIO = 'Bioreactor simulation',
  POLL = 'The chemical reaction part of the air pollution model',
  CLEAR = 'Clear model',
  TO_JS = 'Export model to JavaScript script',
  APP = 'Export model to platform application with user interface',
  OPEN_DS = 'Open model in Diff Studio',
  SAVE = 'Save changes',
  SENS_AN = 'Run sensitivity analysis',
  FITTING = 'Run fitting inputs',
  CONTINUE = 'Continue computations',
  ABORT = 'Abort computations',
  MAX_TIME = `Max computation time, ${COMPUTATION_TIME_UNITS}.`,
  CLICK_RUN = `Click to run`,
}; // HINT

/** UI titles */
export enum TITLE {
  LOAD = 'Load...',
  FROM_FILE = 'From file...',
  SAVE_TO = 'Save to',
  LOCAL_FILE = 'Local file...',
  MY_FILES = 'My files...',
  TO_MY_FILES = 'Save to My Files...',
  AS_LOCAL = 'Save as Local File...',
  TEMPL = 'Templates',
  BASIC = 'Basic',
  ADV = 'Advanced',
  EXT = 'Extended',
  EXAMP = 'Examples',
  RECENT = 'Recent',
  CHEM = 'Chem reactions',
  ROB = `Robertson's model`,
  FERM = 'Fermentation',
  PK = 'PK',
  PKPD = 'PK-PD',
  ACID = 'Acid production',
  NIM = 'Nimotuzumab',
  BIO = 'Bioreactor',
  POLL = 'Pollution',
  CLEAR = 'Clear',
  TO_JS = 'js',
  MISC = 'Misc',
  VARY = 'Vary inputs',
  MODEL = 'Model',
  IPUTS = 'Run',
  FIT = 'Fit',
  SOLUTION = 'Solution',
  OPEN = 'Open',
  BROWSE = 'Browse',
  SAVE = 'Save',
  APPS = 'Apps',
  DIF_ST = 'Diff Studio',
  NAME = 'Name',
  TYPE = 'Type',
  INFO = 'Info',
  IS_CUST = 'Custom model',
}; // TITLE

/** Titles of template models */
export const TEMPLATE_TITLES = [TITLE.BASIC, TITLE.ADV, TITLE.EXT];

/** Titles of example models */
export const EXAMPLE_TITLES = [TITLE.CHEM, TITLE.ROB, TITLE.FERM, TITLE.PK,
  TITLE.PKPD, TITLE.ACID, TITLE.NIM, TITLE.BIO, TITLE.POLL];

/** Models' tooltips map */
export const MODEL_HINT = new Map([
  [TITLE.BASIC, HINT.BASIC],
  [TITLE.ADV, HINT.ADV],
  [TITLE.EXT, HINT.EXT],
  [TITLE.CHEM, HINT.CHEM],
  [TITLE.ROB, HINT.ROB],
  [TITLE.FERM, HINT.FERM],
  [TITLE.PK, HINT.PK],
  [TITLE.PKPD, HINT.PKPD],
  [TITLE.ACID, HINT.ACID],
  [TITLE.NIM, HINT.NIM],
  [TITLE.BIO, HINT.BIO],
  [TITLE.POLL, HINT.POLL],
]);

/** Help links */
export enum LINK {
  DIF_STUDIO_REL = '/help/compute/diff-studio',
  DIF_STUDIO = 'https://datagrok.ai/help/compute/diff-studio',
  SENS_AN = 'https://datagrok.ai/help/compute#sensitivity-analysis',
  FITTING = 'https://datagrok.ai/help/compute#input-parameter-optimization',
  CHEM_REACT = `${DIF_STUDIO_REL}#chem-reactions`,
  FERMENTATION = `${DIF_STUDIO_REL}#fermentation`,
  GA_PRODUCTION = `${DIF_STUDIO_REL}#acid-production`,
  NIMOTUZUMAB = `${DIF_STUDIO_REL}#nimotuzumab`,
  PK = `${DIF_STUDIO_REL}#pk`,
  PKPD = `${DIF_STUDIO_REL}#pk-pd`,
  ROBERTSON = `${DIF_STUDIO_REL}#robertson-model`,
  BIOREACTOR = `${DIF_STUDIO_REL}#bioreactor`,
  POLLUTION = `${DIF_STUDIO_REL}#pollution`,
};

/** Error messages */
export enum ERROR_MSG {
  SOLVING_FAILS = 'Solving fails',
  APP_CREATING_FAILS = 'Application creating fails',
  EXPORT_TO_SCRIPT_FAILS = 'Export to JavaScript script fails',
  SCRIPTING_ISSUE = 'Platform scripting issue',
  UI_ISSUE = 'UI creating issue',
  MISSING_CLOSING_BRACKET = 'Annotation: "]" is missing',
  INCORRECT_BRACES_USE = 'Annotation: incorrect use of "{}"',
  MISSING_COLON = 'Annotation: ":" is missing',
  CHECK_ARGUMENTS = ' (check the "argument" section)',
  INCORRECT_ARG_SEGM = 'Incorrect limits for the argument',
  OPEN_IN_DIF_STUD = 'To change the model, open it in Diff Studio.',
  FAILED_TO_SAVE = 'Failed to save model to the file',
  INCORRECT_MODEL = 'Incorrect model',
  SENS_AN_FAILS = 'Sensitivity Analysis fails',
  FITTING_FAILS = 'Fitting fails',
  PLATFORM_ISSUE = 'Platform issue',
};

/** Lookup table fails */
export enum LOOKUP_DF_FAIL {
  LOAD = 'Failed to load lookup table: ',
  PLATFORM = 'the platform issue',
  FUNCTION = 'incorrect function',
  NO_DF = 'no dataframe',
  INCORRECT = 'incorrect dataframe, ',
  ROWS = 'at least one row is needed.',
  NULLS = 'missing values are not allowed.',
  NUMS = 'no numerical columns.',
  CHOICES = 'first column must contain strings.',
};

/** Lookup table specification fails */
export enum LOOKUP_EXPR_FAIL {
  MISSING = 'Incorrect input: missing ',
};

/** Other UI constants */
export enum MISC {
  VIEW_DEFAULT_NAME = 'Template',
  IVP_EXT = 'ivp',
  TXT_EXT = 'txt',
  FILE_DEFAULT_NAME = `equations.${IVP_EXT}`,
  DEFAULT = 'Default',
  CHOICES = 'choices',
  NAME = 'name',
  CATEGORY = 'category',
  CAPTION = 'caption',
  TOOLTIP = 'tooltip',
};

/** Warning dialog lines */
export enum WARNING {
  TITLE = 'WARNING',
  CHECK = 'Show this warning',
  OVERWRITE_MODEL = 'Overwrite the current model?',
  OVERWRITE_FILE = 'Overwrite existing file?',
  CONTINUE = 'Continue?',
  CHECK_PERF = 'Check time',
  TIME_LIM = 'Time limit',
  UNITS = COMPUTATION_TIME_UNITS,
  PREVIEW = `Model preview is unavailble for this type. Use "${MISC.IVP_EXT}" instead`,
};

/** Code completion infos */
export enum INFO {
  NAME = 'name of the model',
  TAGS = 'scripting tags',
  DESCR = 'descritpion of the model',
  DIF_EQ = 'block of differential equation(s) specification',
  EXPR = 'block of auxiliary expressions & computations',
  ARG = 'independent variable specification',
  INITS = 'initial values of the model',
  PARAMS = 'parameters of the model',
  CONSTS = 'constants definition',
  TOL = 'tolerance of numerical solution',
  EPS_SCALE = 'scale used when computing Jacobian',
  LOOP = 'loop feature',
  UPDATE = 'update model feature',
  OUTPUT = 'output specification',
  COMMENT = 'block with comments',
  SOLVER = 'solver options',
  INPUS = 'path to table with inputs',
};

/** Demo app help info */
export const demoInfo = `# Try
Modify formulas and go to the **${TITLE.IPUTS}** tab.

# No-code
Define equations in a declarative form.

# Interactivity
Play with model inputs on the **${TITLE.IPUTS}** tab.

# Examples
Click <i class="fas fa-folder-open"></i> **Open** icon and explore **Examples**.

# Scripting
Click **JS** button and export model to JavaScript script.

# Analysis
* Click <i class="fas fa-analytics"></i> **Run sensitivity analysis** icon to explore model
* Click <i class="fas fa-chart-line"></i> **Run fitting inputs** icon to optimize inputs

# Learn more
* [Diff Studio](${LINK.DIF_STUDIO})
* [Sensitivity Analysis](${LINK.SENS_AN})
* [Parameter Optimization](${LINK.FITTING})`;

/** Inputs types */
export enum INPUT_TYPE {
  FLOAT = 'Float',
  INT = 'Int',
};

/** Path related consts */
export enum PATH {
  APP_DATA_DS = '/files/system.appdata/diffstudio',
  APPS_DS = '/apps/DiffStudio',
  MODEL = `?model=`,
  CUSTOM = `${MODEL}custom`,
  EMPTY = `${MODEL}empty`,
  EQ = '=',
  AND = '&',
  PARAM = `?params:`,
  BROWSE = 'browse',
  RECENT = 'diff-studio-recent.d42',
  MY_FILES = 'Myfiles',
  HOME = 'Home',
  SYSTEM = 'System',
};

/** UI time constants */
export enum UI_TIME {
  DOCK_EDITOR_TIMEOUT = 100,
  PREVIEW_RUN_SOLVING = 105,
  APP_RUN_SOLVING = 100,
  SOLV_DEFAULT_TIME_SEC = 5,
  SOLV_TIME_MIN_SEC = 1,
  BROWSING = APP_RUN_SOLVING + 300,
  SWITCH_TO_FOLDER = 1,
};

/** Numerical methods names */
export enum METHOD {
  MRT = 'mrt',
  ROS3PRw = 'ros3prw',
  ROS34PRw = 'ros34prw',
};

export const DOCK_RATIO = 0.3;

export const MAX_RECENT_COUNT = 10;

export const CUSTOM_MODEL_IMAGE_LINK = 'images/custom.png';

/** Model image link */
export const modelImageLink = new Map([
  [TITLE.BASIC, 'images/basic.png'],
  [TITLE.ADV, 'images/advanced.png'],
  [TITLE.EXT, 'images/extended.png'],
  [TITLE.CHEM, 'images/chem-react.png'],
  [TITLE.ROB, 'images/robertson.png'],
  [TITLE.FERM, 'images/fermentation.png'],
  [TITLE.PK, 'images/pk.png'],
  [TITLE.PKPD, 'images/pk-pd.png'],
  [TITLE.ACID, 'images/ga-production.png'],
  [TITLE.NIM, 'images/nimotuzumab.png'],
  [TITLE.BIO, 'images/bioreactor.png'],
  [TITLE.POLL, 'images/pollution.png'],
]);

/** Inputs table constants */
export enum INPUTS_DF {
  MIN_ROWS_COUNT = 1,
  INP_NAMES_IDX = 0,
  INPUT_SETS_COL_IDX = 0,
};
