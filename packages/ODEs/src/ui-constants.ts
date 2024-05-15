/** Diff Studio application UI constants */

/** Hot keys */
export enum HOT_KEY {
  RUN = 'F5',
};

/** Tooltips messages */
export enum HINT {
  HELP = 'Open help in a new tab',
  OPEN = 'Open model',
  SAVE_LOC = 'Save model to local file',
  LOAD = 'Load model from local file',
  BASIC = 'Open basic template',
  ADV = 'Open advanced template',
  EXT = 'Open extended template',
  CHEM = 'Mass-action kinetics illustration',
  ROB = `Robertson's chemical reaction model`,
  FERM = 'Fermentation process simulation',
  PK = 'Pharmacokinetic model',
  PKPD = 'Pharmacokinetic-pharmacodynamic model',
  ACID = 'Gluconic acid production model',
  NIM = 'Nimotuzumab disposition model',
  BIO = 'Bioreactor simulation',
  CLEAR = 'Clear model',
  TO_JS = 'Export model to JavaScript script',
  APP = 'Export model to platform application with user interface',
  OPEN_DS = 'Open model in Diff Studio',
  SAVE = 'Save changes',
  SENS_AN = 'Run sensitivity analysis',
  FITTING = 'Run fitting inputs',
};

/** UI titles */
export enum TITLE {
  SAVE_DOTS = 'Save...',
  LOAD = 'Load...',
  FROM_FILE = 'From file...',
  TEMPL = 'Templates',
  BASIC = 'Basic',
  ADV = 'Advanced',
  EXT = 'Extended',
  CASES = 'Examples',
  CHEM = 'Chem reactions',
  ROB = `Robertson's model`,
  FERM = 'Fermentation',
  PK = 'PK',
  PKPD = 'PK-PD',
  ACID = 'Acid production',
  NIM = 'Nimotuzumab',
  BIO = 'Bioreactor',
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
};

/** Help links */
export enum LINK {
  DIF_STUDIO_REL = '/help/compute/diff-studio',
  DIF_STUDIO = 'https://datagrok.ai/help/compute/diff-studio',
  SENS_AN = 'https://datagrok.ai/help/compute#sensitivity-analysis',
  CHEM_REACT = `${DIF_STUDIO_REL}#chem-reactions`,
  FERMENTATION = `${DIF_STUDIO_REL}#fermentation`,
  GA_PRODUCTION = `${DIF_STUDIO_REL}#acid-production`,
  NIMOTUZUMAB = `${DIF_STUDIO_REL}#nimotuzumab`,
  PK = `${DIF_STUDIO_REL}#pk`,
  PKPD = `${DIF_STUDIO_REL}#pk-pd`,
  ROBERTSON = `${DIF_STUDIO_REL}#robertson-model`,
  BIOREACTOR = `${DIF_STUDIO_REL}#bioreactor`,
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
};

/** Warning dialog lines */
export enum WARNING {
  TITLE = 'WARNING',
  CHECK = 'Show this warning',
  MES = 'Overwrite the current model?',
};

/** Other UI constants */
export enum MISC {
  VIEW_DEFAULT_NAME = 'Template',
  FILE_EXT = 'ivp',
  FILE_DEFAULT_NAME = `equations.${FILE_EXT}`,
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
  LOOP = 'loop feature',
  UPDATE = 'update model feature',
  OUTPUT = 'output specification',
  COMMENT = 'block with comments',
};

/** Demo app help info */
export const demoInfo = `# Try
Modify formulas and go to the **${TITLE.IPUTS}** tab.

# No-code
Define equations in a declarative form.

# Interactivity
Play with model inputs on the **${TITLE.IPUTS}** tab.

# Examples
Press <i class="fas fa-folder-open"></i> **Open** icon and explore **Examples**.

# Scripting
Press **JS** button and export model to JavaScript script.

# Analysis
Press <i class="fas fa-analytics"></i> **Run sensitivity analysis** icon to explore model

# Learn more
* [Diff Studio](${LINK.DIF_STUDIO})
* [Sensitivity Analysis](${LINK.SENS_AN})`;

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
};

/** UI timeouts */
export enum TIMEOUT {
  PREVIEW_DOCK_EDITOR = 1000,
  PREVIEW_RUN_SOLVING = 1100,
  APP_RUN_SOLVING = 1100,
};
