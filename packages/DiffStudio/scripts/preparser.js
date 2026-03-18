
//name: ivpLanguagePreparser
//language: nodejs
//input: string code
//output: string result
;// ./src/constants.ts
// Control constants
const CONTROL_TAG = '#';
const CONTROL_TAG_LEN = CONTROL_TAG.length;
const constants_DF_NAME = 'df';
const META = (/* unused pure expression or super */ null && (`${CONTROL_TAG}meta`));
const MAX_LINE_CHART = 4;
/** Control expressions for the problem specifying */
var CONTROL_EXPR;
(function (CONTROL_EXPR) {
    CONTROL_EXPR["NAME"] = "#name";
    CONTROL_EXPR["TAGS"] = "#tags";
    CONTROL_EXPR["DESCR"] = "#description";
    CONTROL_EXPR["DIF_EQ"] = "#equations";
    CONTROL_EXPR["EXPR"] = "#expressions";
    CONTROL_EXPR["ARG"] = "#argument";
    CONTROL_EXPR["INITS"] = "#inits";
    CONTROL_EXPR["CONSTS"] = "#constants";
    CONTROL_EXPR["PARAMS"] = "#parameters";
    CONTROL_EXPR["TOL"] = "#tolerance";
    CONTROL_EXPR["LOOP"] = "#loop";
    CONTROL_EXPR["UPDATE"] = "#update";
    CONTROL_EXPR["RUN_ON_OPEN"] = "#meta.runOnOpen";
    CONTROL_EXPR["RUN_ON_INPUT"] = "#meta.runOnInput";
    CONTROL_EXPR["OUTPUT"] = "#output";
    CONTROL_EXPR["COMMENT"] = "#comment";
    CONTROL_EXPR["SOLVER"] = "#meta.solver";
    CONTROL_EXPR["INPUTS"] = "#meta.inputs";
})(CONTROL_EXPR || (CONTROL_EXPR = {}));
;
/** Loop consts */
var constants_LOOP;
(function (LOOP) {
    LOOP[LOOP["MIN_LINES_COUNT"] = 1] = "MIN_LINES_COUNT";
    LOOP[LOOP["COUNT_IDX"] = 0] = "COUNT_IDX";
    LOOP["COUNT_NAME"] = "_count";
    LOOP[LOOP["MIN_COUNT"] = 1] = "MIN_COUNT";
})(constants_LOOP || (constants_LOOP = {}));
;
/** UPDATE consts */
var constants_UPDATE;
(function (UPDATE) {
    UPDATE[UPDATE["MIN_LINES_COUNT"] = 1] = "MIN_LINES_COUNT";
    UPDATE[UPDATE["DURATION_IDX"] = 0] = "DURATION_IDX";
    UPDATE["DURATION"] = "_duration";
})(constants_UPDATE || (constants_UPDATE = {}));
;
/** Ranges of the solver options */
const SOLVER_OPTIONS_RANGES = new Map([
    ['maxTime', { min: 1, max: 10000 }],
    ['scale', { min: 0.5, max: 1 }],
]);
const constants_TINY = 0.0001;
const constants_STEP_RATIO = 0.5;

;// ./src/ui-constants.ts
/** Diff Studio application UI constants */
/** Hot keys */
var HOT_KEY;
(function (HOT_KEY) {
    HOT_KEY["RUN"] = "F5";
})(HOT_KEY || (HOT_KEY = {}));
;
const COMPUTATION_TIME_UNITS = 'sec';
/** Tooltips messages */
var HINT;
(function (HINT) {
    HINT["HELP"] = "Open help in a new tab";
    HINT["OPEN"] = "Open model";
    HINT["SAVE_MODEL"] = "Save model";
    HINT["SAVE_LOC"] = "Save model to local file";
    HINT["SAVE_MY"] = "Save model to My files";
    HINT["LOAD"] = "Load model from local file";
    HINT["BASIC"] = "Basic template";
    HINT["ADV"] = "Advanced template";
    HINT["EXT"] = "Extended template";
    HINT["CHEM"] = "Mass-action kinetics illustration";
    HINT["ROB"] = "Robertson's chemical reaction model";
    HINT["FERM"] = "Fermentation process simulation";
    HINT["PK"] = "Pharmacokinetic model";
    HINT["PKPD"] = "Pharmacokinetic-pharmacodynamic model";
    HINT["ACID"] = "Gluconic acid production model";
    HINT["NIM"] = "Nimotuzumab disposition model";
    HINT["BIO"] = "Bioreactor simulation";
    HINT["POLL"] = "The chemical reaction part of the air pollution model";
    HINT["CLEAR"] = "Clear model";
    HINT["TO_JS"] = "Open in script editor";
    HINT["APP"] = "Export model to platform application with user interface";
    HINT["OPEN_DS"] = "Open model in Diff Studio";
    HINT["SAVE"] = "Save changes";
    HINT["SENS_AN"] = "Run sensitivity analysis";
    HINT["FITTING"] = "Run fitting inputs";
    HINT["CONTINUE"] = "Continue computations";
    HINT["ABORT"] = "Abort computations";
    HINT["MAX_TIME"] = "Max computation time, sec.";
    HINT["CLICK_RUN"] = "Click to run";
    HINT["SOLVE"] = "Solve equations (F5)";
    HINT["NO_MY_MODELS"] = "No models in My files";
    HINT["FAILED_TO_LOAD_MY_MODELS"] = "Failed to load models from My files: file system error";
    HINT["NO_RECENT_MODELS"] = "No recent models";
    HINT["FAILED_TO_LOAD_RECENT_MODELS"] = "Failed to load recent models: file system error";
    HINT["CORRUPTED_DATA_FILE"] = "Failed to load recent models: corrupted data file";
    HINT["EDIT"] = "Edit";
    HINT["UPDATE"] = "Apply all changes to the current file";
})(HINT || (HINT = {}));
; // HINT
/** UI titles */
var TITLE;
(function (TITLE) {
    TITLE["LOAD"] = "Load...";
    TITLE["IMPORT"] = "Import...";
    TITLE["SAVE_TO_MY_FILES"] = "Save to My files";
    TITLE["LOCAL_FILE"] = "Local file...";
    TITLE["MY_FILES"] = "My files...";
    TITLE["TO_MY_FILES"] = "Save to My Files...";
    TITLE["AS_LOCAL"] = "Save as Local File...";
    TITLE["TEMPL"] = "Templates";
    TITLE["BASIC"] = "Basic";
    TITLE["ADV"] = "Advanced";
    TITLE["EXT"] = "Extended";
    TITLE["LIBRARY"] = "Library";
    TITLE["RECENT"] = "Recent";
    TITLE["CHEM"] = "Chem reactions";
    TITLE["ROB"] = "Robertson's model";
    TITLE["FERM"] = "Fermentation";
    TITLE["PK"] = "PK";
    TITLE["PKPD"] = "PK-PD";
    TITLE["ACID"] = "Acid production";
    TITLE["NIM"] = "Nimotuzumab";
    TITLE["BIO"] = "Bioreactor";
    TITLE["POLL"] = "Pollution";
    TITLE["CLEAR"] = "Clear";
    TITLE["TO_JS"] = "js";
    TITLE["MISC"] = "Misc";
    TITLE["VARY"] = "Vary inputs";
    TITLE["EDIT"] = "Edit";
    TITLE["SOLVE"] = "Solve";
    TITLE["FIT"] = "Fit";
    TITLE["SOLUTION"] = "Solution";
    TITLE["OPEN"] = "Open";
    TITLE["BROWSE"] = "Browse";
    TITLE["SAVE"] = "Save";
    TITLE["APPS"] = "Apps";
    TITLE["COMP"] = "Compute";
    TITLE["DIF_ST"] = "Diff Studio";
    TITLE["NAME"] = "Name";
    TITLE["TYPE"] = "Type";
    TITLE["INFO"] = "Info";
    TITLE["IS_CUST"] = "Custom model";
    TITLE["MY_MODELS"] = "My Models";
    TITLE["NO_MODELS"] = "None";
    TITLE["MULTI_AXIS"] = "Multiaxis";
    TITLE["FACET"] = "Facet";
    TITLE["UPDATE"] = "Update";
    TITLE["CONTROLS"] = "Model controls";
})(TITLE || (TITLE = {}));
; // TITLE
/** Titles of template models */
const TEMPLATE_TITLES = [TITLE.BASIC, TITLE.ADV, TITLE.EXT];
/** Titles of example models */
const EXAMPLE_TITLES = [TITLE.CHEM, TITLE.ROB, TITLE.FERM, TITLE.PK,
    TITLE.PKPD, TITLE.ACID, TITLE.NIM, TITLE.BIO, TITLE.POLL];
/** Models' tooltips map */
const MODEL_HINT = new Map([
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
var LINK;
(function (LINK) {
    LINK["DIF_STUDIO_REL"] = "/help/compute/diff-studio";
    LINK["MODELS"] = "/help/compute/models";
    LINK["DIF_STUDIO"] = "https://datagrok.ai/help/compute/diff-studio";
    LINK["SENS_AN"] = "/help/compute/function-analysis#sensitivity-analysis";
    LINK["FITTING"] = "/help/compute/function-analysis#parameter-optimization";
    LINK["CHEM_REACT"] = "/help/compute/models#chem-reactions";
    LINK["FERMENTATION"] = "/help/compute/models#fermentation";
    LINK["GA_PRODUCTION"] = "/help/compute/models#acid-production";
    LINK["NIMOTUZUMAB"] = "/help/compute/models#nimotuzumab";
    LINK["PK"] = "/help/compute/models#pk";
    LINK["PKPD"] = "/help/compute/models#pk-pd";
    LINK["ROBERTSON"] = "/help/compute/models#robertson-model";
    LINK["BIOREACTOR"] = "/help/compute/models#bioreactor";
    LINK["POLLUTION"] = "/help/compute/models#pollution";
    LINK["COMPUTE"] = "https://datagrok.ai/help/compute";
    LINK["LOAD_SAVE"] = "/help/compute/diff-studio#working-with-models";
    LINK["MODEL_COMPONENTS"] = "/help/compute/diff-studio#model-components-and-syntax";
    LINK["INTERFACE"] = "/help/compute/diff-studio#user-interface-options";
})(LINK || (LINK = {}));
;
/** Error messages */
var ERROR_MSG;
(function (ERROR_MSG) {
    ERROR_MSG["SOLVING_FAILS"] = "Solving fails";
    ERROR_MSG["APP_CREATING_FAILS"] = "Application creating fails";
    ERROR_MSG["EXPORT_TO_SCRIPT_FAILS"] = "Export to JavaScript script fails";
    ERROR_MSG["SCRIPTING_ISSUE"] = "Platform scripting issue";
    ERROR_MSG["UI_ISSUE"] = "UI creating issue";
    ERROR_MSG["MISSING_CLOSING_BRACKET"] = "Annotation error: **\"]\"** is missing";
    ERROR_MSG["INCORRECT_BRACES_USE"] = "Annotation error: incorrect use of **\"{}\"**";
    ERROR_MSG["MISSING_COLON"] = "Annotation error: **\":\"** is missing";
    ERROR_MSG["CHECK_ARGUMENTS"] = " (check the **#argument** block)";
    ERROR_MSG["INCORRECT_ARG_SEGM"] = "Incorrect limits for the argument";
    ERROR_MSG["OPEN_IN_DIF_STUD"] = "To change the model, open it in Diff Studio.";
    ERROR_MSG["FAILED_TO_SAVE"] = "Failed to save model to the file";
    ERROR_MSG["INCORRECT_MODEL"] = "Incorrect model";
    ERROR_MSG["SENS_AN_FAILS"] = "Sensitivity Analysis fails";
    ERROR_MSG["FITTING_FAILS"] = "Fitting fails";
    ERROR_MSG["PLATFORM_ISSUE"] = "Platform issue";
    ERROR_MSG["INCORTECT_INPUT"] = "Incorrect input caption";
    ERROR_MSG["NANS_OBTAINED"] = "Failed computations: obtained NaN-s. Check inputs and formulas.";
})(ERROR_MSG || (ERROR_MSG = {}));
;
/** Lookup table fails */
var LOOKUP_DF_FAIL;
(function (LOOKUP_DF_FAIL) {
    LOOKUP_DF_FAIL["LOAD"] = "Failed to load lookup table: ";
    LOOKUP_DF_FAIL["PLATFORM"] = "the platform issue";
    LOOKUP_DF_FAIL["FUNCTION"] = "incorrect function";
    LOOKUP_DF_FAIL["NO_DF"] = "no dataframe";
    LOOKUP_DF_FAIL["INCORRECT"] = "incorrect dataframe, ";
    LOOKUP_DF_FAIL["ROWS"] = "at least one row is needed.";
    LOOKUP_DF_FAIL["NULLS"] = "missing values are not allowed.";
    LOOKUP_DF_FAIL["NUMS"] = "no numerical columns.";
    LOOKUP_DF_FAIL["CHOICES"] = "first column must contain strings.";
})(LOOKUP_DF_FAIL || (LOOKUP_DF_FAIL = {}));
;
/** Lookup table specification fails */
var LOOKUP_EXPR_FAIL;
(function (LOOKUP_EXPR_FAIL) {
    LOOKUP_EXPR_FAIL["MISSING"] = "Incorrect input: missing ";
})(LOOKUP_EXPR_FAIL || (LOOKUP_EXPR_FAIL = {}));
;
/** Other UI constants */
var MISC;
(function (MISC) {
    MISC["VIEW_DEFAULT_NAME"] = "Template";
    MISC["MODEL_FILE_EXT"] = "ivp";
    MISC["FILE_DEFAULT_NAME"] = "equations.ivp";
    MISC["DEFAULT"] = "Default";
    MISC["CHOICES"] = "choices";
    MISC["NAME"] = "name";
    MISC["CATEGORY"] = "category";
    MISC["CAPTION"] = "caption";
    MISC["TOOLTIP"] = "tooltip";
    MISC["IS_NOT_DEF"] = "is not defined";
    MISC["UNEXPECTED"] = "Unexpected identifier";
    MISC["PROP_OF_NULL"] = "Cannot set properties of null";
})(MISC || (MISC = {}));
;
/** Warning dialog lines */
var WARNING;
(function (WARNING) {
    WARNING["TITLE"] = "WARNING";
    WARNING["CHECK"] = "Show this warning";
    WARNING["OVERWRITE_MODEL"] = "Overwrite the current model?";
    WARNING["OVERWRITE_FILE"] = "Overwrite existing file?";
    WARNING["CONTINUE"] = "Continue?";
    WARNING["CHECK_PERF"] = "Check time";
    WARNING["TIME_LIM"] = "Time limit";
    WARNING["UNITS"] = "sec";
    WARNING["PREVIEW"] = "Model preview is unavailble for this type. Use \"ivp\" instead";
})(WARNING || (WARNING = {}));
;
/** Code completion infos */
var INFO;
(function (INFO) {
    INFO["NAME"] = "name of the model";
    INFO["TAGS"] = "scripting tags";
    INFO["DESCR"] = "descritpion of the model";
    INFO["DIF_EQ"] = "block of differential equation(s) specification";
    INFO["EXPR"] = "block of auxiliary expressions & computations";
    INFO["ARG"] = "independent variable specification";
    INFO["INITS"] = "initial values of the model";
    INFO["PARAMS"] = "parameters of the model";
    INFO["CONSTS"] = "constants definition";
    INFO["TOL"] = "tolerance of numerical solution";
    INFO["EPS_SCALE"] = "scale used when computing Jacobian";
    INFO["LOOP"] = "loop feature";
    INFO["UPDATE"] = "update model feature";
    INFO["OUTPUT"] = "output specification";
    INFO["COMMENT"] = "block with comments";
    INFO["SOLVER"] = "solver options";
    INFO["INPUS"] = "path to table with inputs";
})(INFO || (INFO = {}));
;
/** Demo app help info */
const demoInfo = (/* unused pure expression or super */ null && (`# Diff Studio
In-browser solver of ordinary differential equations 
([ODEs](https://en.wikipedia.org/wiki/Ordinary_differential_equation))

# Interactivity
Play with the model inputs. Move sliders to explore its behavior.

# Model
Turn on the **${TITLE.EDIT}** toggle on the top panel to modify formulas.

# No-code
Define equations in a declarative form.

# Examples
Click <i class="fas fa-folder-open d4-combo-popup" style="min-width: 0px; cursor: default"></i> icon and
explore **Library**.

# Scripting
Click **</>** icon to export model to JavaScript script.

# Analysis
Turn off the **${TITLE.EDIT}** toggle, and perform analysis:
* Click the **Fit** icon on the top panel to [optimize inputs](${LINK.FITTING}).
* Click the **Sensitivity** icon to run [sensitivity analysis](${LINK.SENS_AN}).

# Catalog
Click <i class="fas fa-layer-plus"></i> icon to save model to **Model Hub**.

# Learn more
* [Diff Studio](${LINK.DIF_STUDIO})
* [Compute](${LINK.COMPUTE})`));
/** Inputs types */
var INPUT_TYPE;
(function (INPUT_TYPE) {
    INPUT_TYPE["FLOAT"] = "Float";
    INPUT_TYPE["INT"] = "Int";
})(INPUT_TYPE || (INPUT_TYPE = {}));
;
/** Path related consts */
var PATH;
(function (PATH) {
    PATH["APP_DATA_DS"] = "/files/system.appdata/diffstudio";
    PATH["APPS_DS"] = "/apps/DiffStudio";
    PATH["MODEL"] = "?model=";
    PATH["CUSTOM"] = "?model=custom";
    PATH["EMPTY"] = "?model=empty";
    PATH["EQ"] = "=";
    PATH["AND"] = "&";
    PATH["PARAM"] = "?params:";
    PATH["BROWSE"] = "browse";
    PATH["RECENT"] = "diff-studio-recent.d42";
    PATH["MY_FILES"] = "Myfiles";
    PATH["HOME"] = "Home";
    PATH["SYSTEM"] = "System";
    PATH["FILE"] = "file";
    PATH["SLASH"] = "/";
})(PATH || (PATH = {}));
;
/** UI time constants */
var UI_TIME;
(function (UI_TIME) {
    UI_TIME[UI_TIME["DOCK_EDITOR_TIMEOUT"] = 100] = "DOCK_EDITOR_TIMEOUT";
    UI_TIME[UI_TIME["PREVIEW_RUN_SOLVING"] = 105] = "PREVIEW_RUN_SOLVING";
    UI_TIME[UI_TIME["APP_RUN_SOLVING"] = 100] = "APP_RUN_SOLVING";
    UI_TIME[UI_TIME["SOLV_DEFAULT_TIME_SEC"] = 5] = "SOLV_DEFAULT_TIME_SEC";
    UI_TIME[UI_TIME["SOLV_TIME_MIN_SEC"] = 1] = "SOLV_TIME_MIN_SEC";
    UI_TIME[UI_TIME["BROWSING"] = 600] = "BROWSING";
    UI_TIME[UI_TIME["SWITCH_TO_FOLDER"] = 100] = "SWITCH_TO_FOLDER";
    UI_TIME[UI_TIME["WGT_CLICK"] = 10] = "WGT_CLICK";
    UI_TIME[UI_TIME["FACET_DOCKING"] = 100] = "FACET_DOCKING";
    UI_TIME[UI_TIME["TITLE_REMOVING"] = 500] = "TITLE_REMOVING";
})(UI_TIME || (UI_TIME = {}));
;
/** Numerical methods names */
var METHOD;
(function (METHOD) {
    METHOD["MRT"] = "mrt";
    METHOD["ROS3PRw"] = "ros3prw";
    METHOD["ROS34PRw"] = "ros34prw";
})(METHOD || (METHOD = {}));
;
/** Dock ratios */
var DOCK_RATIO;
(function (DOCK_RATIO) {
    DOCK_RATIO[DOCK_RATIO["INPUTS_TAB"] = 0.2] = "INPUTS_TAB";
    DOCK_RATIO[DOCK_RATIO["GRAPHS"] = 0.8] = "GRAPHS";
})(DOCK_RATIO || (DOCK_RATIO = {}));
;
const MAX_RECENT_COUNT = 10;
const CUSTOM_MODEL_IMAGE_LINK = 'images/custom.png';
/** Model image link */
const modelImageLink = new Map([
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
var INPUTS_DF;
(function (INPUTS_DF) {
    INPUTS_DF[INPUTS_DF["MIN_ROWS_COUNT"] = 1] = "MIN_ROWS_COUNT";
    INPUTS_DF[INPUTS_DF["INP_NAMES_IDX"] = 0] = "INP_NAMES_IDX";
    INPUTS_DF[INPUTS_DF["INPUT_SETS_COL_IDX"] = 0] = "INPUT_SETS_COL_IDX";
})(INPUTS_DF || (INPUTS_DF = {}));
;
/** Max number of facet grid graphs count */
const MAX_FACET_GRAPHS_COUNT = 20;

;// ./src/shared-utils.ts



/** Diff Studio model error */
class ModelError extends Error {
    constructor(message, helpUrl, toHighlight) {
        super(message);
        this.toHighlight = undefined;
        this.helpUrl = helpUrl;
        this.toHighlight = toHighlight;
    }
    getHelpUrl() {
        return this.helpUrl;
    }
    getToHighlight() {
        return this.toHighlight;
    }
}
; // ModelError
/**  String-to-value */
const strToVal = (s) => {
    const num = Number(s);
    return !isNaN(num) ? num : s === 'true' ? true : s === 'false' ? false : s;
};
/** Copy of getLookupsInfo without UI warnings */
function getLookupsInfoData(inputsLookup) {
    var _a, _b, _c, _d, _e;
    const info = new Map();
    const braceOpenIdx = inputsLookup.indexOf(BRACE_OPEN);
    const braceCloseIdx = inputsLookup.indexOf(BRACE_CLOSE);
    if (braceOpenIdx < 0) {
        return null;
    }
    if (braceCloseIdx < 0) {
        return null;
    }
    // extract name
    info.set('name', inputsLookup.slice(0, braceOpenIdx).replaceAll(' ', ''));
    // extract features
    const options = inputsLookup.slice(braceOpenIdx + 1, braceCloseIdx).split(ANNOT_SEPAR);
    let sepIdx;
    for (const opt of options) {
        sepIdx = opt.indexOf(CONTROL_SEP);
        if (sepIdx < 0) {
            return null;
        }
        info.set(opt.slice(0, sepIdx).trim(), opt.slice(sepIdx + 1).trim());
    }
    // extract tooltip
    const bracketOpenIdx = inputsLookup.indexOf(BRACKET_OPEN);
    if (bracketOpenIdx > 0) {
        const bracketCloseIdx = inputsLookup.indexOf(BRACKET_CLOSE);
        if (bracketCloseIdx < 0) {
            return null;
        }
        info.set(MISC.TOOLTIP, inputsLookup.slice(bracketOpenIdx + 1, bracketCloseIdx));
    }
    if (info.get(MISC.CHOICES) === undefined) {
        return null;
    }
    return {
        name: (_a = info.get(MISC.NAME)) !== null && _a !== void 0 ? _a : '',
        caption: (_b = info.get(MISC.CAPTION)) !== null && _b !== void 0 ? _b : ((_c = info.get(MISC.NAME)) !== null && _c !== void 0 ? _c : ''),
        category: (_d = info.get(MISC.CATEGORY)) !== null && _d !== void 0 ? _d : TITLE.MISC,
        tooltip: (_e = info.get(MISC.TOOLTIP)) !== null && _e !== void 0 ? _e : '',
        choices: info.get(MISC.CHOICES),
    };
} // getLookupsInfoData
/** Return options with respect to the model input specification */
function getOptions(name, modelInput, modelBlock, startingInputs) {
    var _a, _b;
    const options = {
        name: name,
        defaultValue: modelInput.value,
        type: 'double',
        inputType: 'Float',
    };
    if (modelInput.annot !== null) {
        let annot = modelInput.annot;
        let descr = undefined;
        let posOpen = annot.indexOf(BRACKET_OPEN);
        let posClose = annot.indexOf(BRACKET_CLOSE);
        if (posOpen !== -1) {
            if (posClose === -1) {
                throw new ModelError(`${ERROR_MSG.MISSING_CLOSING_BRACKET}. Correct annotation in the **${modelBlock}** block.`, LINK.INTERFACE, annot);
            }
            descr = annot.slice(posOpen + 1, posClose);
            annot = annot.slice(0, posOpen);
        }
        posOpen = annot.indexOf(BRACE_OPEN);
        posClose = annot.indexOf(BRACE_CLOSE);
        if (posOpen >= posClose) {
            throw new ModelError(`${ERROR_MSG.INCORRECT_BRACES_USE}. Correct annotation in the ***${modelBlock}** block.`, LINK.INTERFACE, annot);
        }
        let pos;
        let key;
        let val;
        annot.slice(posOpen + 1, posClose).split(ANNOT_SEPAR).forEach((str) => {
            pos = str.indexOf(CONTROL_SEP);
            if (pos === -1) {
                throw new ModelError(`${ERROR_MSG.MISSING_COLON}. Correct annotation in the **${modelBlock}** block.`, LINK.INTERFACE, annot);
            }
            key = str.slice(0, pos).trim();
            val = str.slice(pos + 1).trim();
            options[key !== 'caption' ? key : 'friendlyName'] = strToVal(val);
        });
        options.description = descr !== null && descr !== void 0 ? descr : '';
        options.friendlyName = (_a = options.friendlyName) !== null && _a !== void 0 ? _a : options.name;
        options.caption = options.friendlyName;
    }
    if (modelBlock === CONTROL_EXPR.ARG) {
        options.friendlyName = options.name;
        options.caption = options.name;
        options.name = ARG_INPUT_KEYS_MAPPING[options.name];
    }
    if (startingInputs) {
        options.defaultValue = (_b = startingInputs
            .get(options.name.replace(' ', '').toLowerCase())) !== null && _b !== void 0 ? _b : options.defaultValue;
        modelInput.value = options.defaultValue;
    }
    return options;
}
; // getOptions

;// ./src/scripting-tools.ts
/* eslint-disable max-len */
/* Scripting tools for the Initial Value Problem (IVP) solver:
     - parser of formulas defining IVP;
     - JS-script generator.

   The parser processes IVP formulas given in the special form (see "Project structure" in README.md).

   JS-script generator creates DATAGROK JavaScript script: annotation & code.
*/


// Scripting specific constants
const CONTROL_SEP = ':';
const COMMA = ',';
const EQUAL_SIGN = '=';
const DIV_SIGN = '/';
const SERVICE = '_';
const BRACE_OPEN = '{';
const BRACE_CLOSE = '}';
const BRACKET_OPEN = '[';
const BRACKET_CLOSE = ']';
const ANNOT_SEPAR = ';';
const DEFAULT_TOL = '0.00005';
const DEFAULT_SOLVER_SETTINGS = '{}';
const COLUMNS = (/* unused pure expression or super */ null && (`${SERVICE}columns`));
const COMMENT_SEQ = '//';
const STAGE_COL_NAME = `${SERVICE}Stage`;
const INCEPTION = 'Inception';
/** Solver package name */
const PACKAGE_NAME = 'DiffStudio';
/** Numerical solver function */
const SOLVER_FUNC = 'solveEquations';
/** Elementary math tools */
const MATH_FUNCS = ['pow', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sqrt', 'exp', 'log', 'sinh', 'cosh', 'tanh'];
const POW_IDX = MATH_FUNCS.indexOf('pow');
const MATH_CONSTS = ['PI', 'E'];
/** Default meta */
const defaultMetas = (/* unused pure expression or super */ null && (`//meta.runOnOpen: true\n//meta.runOnInput: true\n//meta.features: {"sens-analysis": true, "fitting": true}\n`));
/** Input keys of Arg */
const ARG_INPUT_KEYS = ['initial', 'final', 'step'];
const ARG_INPUT_KEYS_MAPPING = {
    'initial': '_t0',
    'final': '_t1',
    'step': '_h'
};
/** Scripting specific constants */
var SCRIPTING;
(function (SCRIPTING) {
    SCRIPTING["ARG_NAME"] = "name";
    SCRIPTING["COUNT"] = "count";
    SCRIPTING["DURATION"] = "duration";
})(SCRIPTING || (SCRIPTING = {}));
;
/** Help links for model errors */
var ERROR_LINK;
(function (ERROR_LINK) {
    ERROR_LINK["MAIN_DOCS"] = "/help/compute/diff-studio";
    ERROR_LINK["CORE_BLOCKS"] = "/help/compute/diff-studio#core-blocks";
    ERROR_LINK["COMPS_SYNTAX"] = "/help/compute/diff-studio#model-components-and-syntax";
    ERROR_LINK["LOOP"] = "/help/compute/diff-studio#cyclic-processes";
    ERROR_LINK["UPDATE"] = "/help/compute/diff-studio#multistage-processes";
    ERROR_LINK["UI_OPTS"] = "/help/compute/diff-studio#user-interface-options";
    ERROR_LINK["LOOP_VS_UPDATE"] = "/help/compute/diff-studio#advanced-features";
    ERROR_LINK["SOLVER_CONFIG"] = "/help/compute/diff-studio#solver-configuration";
    ERROR_LINK["MODEL_PARAMS"] = "/help/compute/diff-studio#model-parameters";
})(ERROR_LINK || (ERROR_LINK = {}));
;
/** Specific error messages */
var scripting_tools_ERROR_MSG;
(function (ERROR_MSG) {
    ERROR_MSG["CTRL_EXPR"] = "Unsupported control expression with the tag **\"#\"**";
    ERROR_MSG["ARG"] = "'The **#argument** block must consist of 3 lines specifying initial and final time, and solution grid step.";
    ERROR_MSG["LOOP"] = "The **#loop** block must contain at least one line.";
    ERROR_MSG["COUNT"] = "Incorrect loop count";
    ERROR_MSG["LOOP_VS_UPDATE"] = "The **#loop** and **#update** blocks cannot be used simultaneously.";
    ERROR_MSG["UPDATE_LINES_COUNT"] = "The **#update** block must contain at least one line.";
    ERROR_MSG["DURATION"] = "Incorrect update duration";
    ERROR_MSG["BRACES"] = " Missing one of the braces (**{**, **}**).";
    ERROR_MSG["COLON"] = "Incorrect position of **\":\"**.";
    ERROR_MSG["CASE_INSENS"] = "Non-unique name (case-insensitive): use different caption for ";
    ERROR_MSG["MISSING_INIT"] = "Correct the **#inits** block.";
    ERROR_MSG["UNDEF_NAME"] = "Model name missing. Specify the model name in the **#name** block.";
    ERROR_MSG["UNDEF_DEQS"] = "Differential equation(s) are required for this model. Add equation(s) under the **#equations** block.";
    ERROR_MSG["UNDEF_INITS"] = "Initial conditions are required for this model. Add initial conditions under the **#inits** block.";
    ERROR_MSG["UNDEF_ARG"] = "Argument specification is required for this model. Specify an argument, its range, and a grid step in the **#argument** block.";
    ERROR_MSG["CORRECT_ARG_LIM"] = "Correct limits in the **#argument** block.";
    ERROR_MSG["INTERVAL"] = "Incorrect range for";
    ERROR_MSG["NEGATIVE_STEP"] = "Solution grid step must be positive. Correct the **#argument** block.";
    ERROR_MSG["INCOR_STEP"] = "Grid step must less than the length of solution interval. Correct the **#argument** block.";
    ERROR_MSG["MISS_COLON"] = "Missing **\":\"**";
    ERROR_MSG["NAN"] = "is not a valid number. Correct the line";
    ERROR_MSG["SERVICE_START"] = "Variable names must not begin with **\"_\"**.";
    ERROR_MSG["REUSE_NAME"] = "Variable reuse (case-insensitive): rename ";
    ERROR_MSG["SOLVER"] = "Incorrect solver options. Correct the **#meta.solver** line.";
})(scripting_tools_ERROR_MSG || (scripting_tools_ERROR_MSG = {})); // ERROR_MSG
/** Datagrok annotations */
var ANNOT;
(function (ANNOT) {
    ANNOT["NAME"] = "//name:";
    ANNOT["DESCR"] = "//description:";
    ANNOT["TAGS"] = "//tags:";
    ANNOT["LANG"] = "//language: javascript";
    ANNOT["DOUBLE_INPUT"] = "//input: double";
    ANNOT["INT_INPUT"] = "//input: int";
    ANNOT["OUTPUT"] = "//output: dataframe df";
    ANNOT["EDITOR"] = "//editor: Compute2:RichFunctionViewEditor";
    ANNOT["CAPTION"] = "caption:";
    ANNOT["ARG_INIT"] = "{caption: Initial; category: Argument}";
    ANNOT["ARG_FIN"] = "{caption: Final; category: Argument}";
    ANNOT["ARG_STEP"] = "{caption: Step; category: Argument}";
    ANNOT["INITS"] = "category: Initial values";
    ANNOT["PARAMS"] = "category: Parameters";
})(ANNOT || (ANNOT = {}));
/** JS-scripting components */
var SCRIPT;
(function (SCRIPT) {
    SCRIPT["CONSTS"] = "// constants";
    SCRIPT["ODE_COM"] = "// the problem definition";
    SCRIPT["ODE"] = "let odes = {";
    SCRIPT["SOLVER_COM"] = "// solve the problem";
    SCRIPT["SOLVER"] = "const solver = await grok.functions.eval('DiffStudio:solveEquations');";
    SCRIPT["PREPARE"] = "let call = solver.prepare({problem: odes, options: opts});";
    SCRIPT["CALL"] = "await call.call();";
    SCRIPT["FUNCTION_CALL"] = "let df = await grok.functions.call('DiffStudio:solveEquations', {problem: odes, options: opts});";
    SCRIPT["OUTPUT"] = "let df = call.getParamValue('df');";
    SCRIPT["SPACE2"] = "  ";
    SCRIPT["SPACE4"] = "    ";
    SCRIPT["SPACE6"] = "      ";
    SCRIPT["SPACE8"] = "        ";
    SCRIPT["FUNC_VALS"] = "// extract function values";
    SCRIPT["EVAL_EXPR"] = "// evaluate expressions";
    SCRIPT["COMP_OUT"] = "// compute output";
    SCRIPT["MATH_FUNC_COM"] = "// used Math-functions";
    SCRIPT["MATH_CONST_COM"] = "// used Math-constants";
    SCRIPT["ONE_STAGE_COM"] = "\n// one stage solution";
    SCRIPT["ONE_STAGE_BEGIN"] = "let _oneStage = async (";
    SCRIPT["ONE_STAGE_END"] = ") => {";
    SCRIPT["ASYNC_OUTPUT"] = "let df = await _oneStage(";
    SCRIPT["RETURN"] = "return df;";
    SCRIPT["RETURN_OUTPUT"] = "return call.getParamValue('df');";
    SCRIPT["EMPTY_OUTPUT"] = "let df = DG.DataFrame.create();";
    SCRIPT["APPEND"] = "df.append(";
    SCRIPT["SOLUTION_DF_COM"] = "// solution dataframe";
    SCRIPT["LOOP_INTERVAL_COM"] = "// loop interval";
    SCRIPT["LOOP_INTERVAL"] = "_interval";
    SCRIPT["LAST_IDX"] = "_lastIdx";
    SCRIPT["UPDATE_COM"] = "// update ";
    SCRIPT["CUSTOM_OUTPUT_COM"] = "// create custom output";
    SCRIPT["CUSTOM_COLUMNS"] = "let _columns = [";
    SCRIPT["ONE_STAGE"] = "await _oneStage(";
    SCRIPT["SEGMENT_COM"] = "// add segment category";
})(SCRIPT || (SCRIPT = {}));
/** Get start of the problem skipping note-lines */
function getStartOfProblemDef(lines) {
    const linesCount = lines.length;
    let idx = 0;
    while (!lines[idx].startsWith(CONTROL_TAG)) {
        ++idx;
        if (idx === linesCount)
            throw new ModelError(scripting_tools_ERROR_MSG.UNDEF_NAME, ERROR_LINK.CORE_BLOCKS);
    }
    return idx;
}
/** Get limits of blocks specifying IVP */
function getBlocks(lines, start) {
    const linesCount = lines.length;
    let beg = start;
    let idx = start + 1;
    const blocks = [];
    while (true) {
        if (idx === linesCount) {
            blocks.push({ begin: beg, end: idx });
            break;
        }
        if (lines[idx].startsWith(CONTROL_TAG)) {
            blocks.push({ begin: beg, end: idx });
            beg = idx;
        }
        ++idx;
    }
    return blocks;
}
/** Concatenate multi-line formulas & expressions */
function concatMultilineFormulas(source) {
    const res = [source[0]];
    let idxRes = 0;
    let idxSource = 1;
    const size = source.length;
    while (idxSource < size) {
        const str = source[idxSource];
        if (!str.includes(EQUAL_SIGN))
            res[idxRes] += str;
        else {
            res.push(str);
            ++idxRes;
        }
        ++idxSource;
    }
    return res;
}
/** Get indeces of math functions used in the current text */
function getUsedMathIds(text, mathIds) {
    const res = [];
    const size = mathIds.length;
    for (let i = 0; i < size; ++i) {
        if (text.includes(`${mathIds[i]}`))
            res.push(i);
    }
    return res;
}
/** Get differential equations */
function getDifEquations(lines) {
    const deqs = new Map();
    const names = [];
    let divIdx = 0;
    let eqIdx = 0;
    for (const line of lines) {
        if (line === undefined)
            throw new ModelError(scripting_tools_ERROR_MSG.UNDEF_DEQS, ERROR_LINK.CORE_BLOCKS, CONTROL_EXPR.DIF_EQ);
        divIdx = line.indexOf(DIV_SIGN);
        eqIdx = line.indexOf(EQUAL_SIGN);
        const name = line.slice(line.indexOf('d') + 1, divIdx).replace(/[() ]/g, '');
        names.push(name);
        deqs.set(name, line.slice(eqIdx + 1).trim());
    }
    return {
        equations: deqs,
        solutionNames: names,
    };
}
/** Get expressions of IVP */
function getExpressions(lines) {
    const exprs = new Map();
    let eqIdx = 0;
    for (const line of lines) {
        if (line === undefined)
            break;
        eqIdx = line.indexOf(EQUAL_SIGN);
        exprs.set(line.slice(0, eqIdx).replaceAll(' ', ''), line.slice(eqIdx + 1).trim());
    }
    return exprs;
}
/** Get input specification */
function getInput(line) {
    const str = line.slice(line.indexOf(EQUAL_SIGN) + 1).trim();
    const braceIdx = str.indexOf(BRACE_OPEN);
    let val;
    if (braceIdx === -1) { // no annotation
        val = Number(str);
        // Check right-hand side
        if (isNaN(val))
            throw new ModelError(`**${str}** ${scripting_tools_ERROR_MSG.NAN}\n **${line}**.`, ERROR_LINK.COMPS_SYNTAX, line);
        return {
            value: val,
            annot: null,
        };
    }
    // There is an annotation
    val = Number(str.slice(0, braceIdx));
    // Check right-hand side
    if (isNaN(val))
        throw new ModelError(`**${str.slice(0, braceIdx)}** ${scripting_tools_ERROR_MSG.NAN}: **${line}**`, ERROR_LINK.COMPS_SYNTAX, line);
    return {
        value: val,
        annot: str.slice(braceIdx),
    };
}
/** Get argument (independent variable) of IVP */
function getArg(lines) {
    if (lines.length !== 4)
        throw new ModelError(scripting_tools_ERROR_MSG.ARG, ERROR_LINK.CORE_BLOCKS, CONTROL_EXPR.ARG);
    const sepIdx = lines[0].indexOf(CONTROL_SEP);
    if (sepIdx === -1)
        throw new ModelError(`${scripting_tools_ERROR_MSG.MISS_COLON}. Correct the line **${lines[0]}**`, ERROR_LINK.CORE_BLOCKS, lines[0]);
    const commaIdx = lines[0].indexOf(COMMA);
    return {
        name: ((commaIdx > -1) ? lines[0].slice(sepIdx + 1, commaIdx) : lines[0].slice(sepIdx + 1)).trim(),
        initial: getInput(lines[1]),
        final: getInput(lines[2]),
        step: getInput(lines[3]),
        updateName: (commaIdx > -1) ? lines[0].slice(commaIdx + 1).trim() : undefined,
    };
}
/** Get equalities specification */
function getEqualities(lines, begin, end) {
    const source = concatMultilineFormulas(lines.slice(begin, end));
    const eqs = new Map();
    let eqIdx = 0;
    for (const line of source) {
        if (line === undefined)
            break;
        eqIdx = line.indexOf(EQUAL_SIGN);
        eqs.set(line.slice(0, eqIdx).replaceAll(' ', ''), getInput(line));
    }
    return eqs;
}
/** Get loop specification */
function getLoop(lines, begin, end) {
    if (begin >= end)
        throw new ModelError(scripting_tools_ERROR_MSG.LOOP, ERROR_LINK.LOOP, CONTROL_EXPR.LOOP);
    const source = concatMultilineFormulas(lines.slice(begin, end));
    const size = source.length;
    if (size < constants_LOOP.MIN_LINES_COUNT)
        throw new ModelError(scripting_tools_ERROR_MSG.LOOP, ERROR_LINK.LOOP, CONTROL_EXPR.LOOP);
    const count = getInput(source[constants_LOOP.COUNT_IDX]);
    if (count.value < constants_LOOP.MIN_COUNT)
        throw new ModelError(`${scripting_tools_ERROR_MSG.COUNT}: **${count.value}**.`, ERROR_LINK.LOOP, source[constants_LOOP.COUNT_IDX]);
    return { count: count, updates: source.slice(constants_LOOP.COUNT_IDX + 1) };
}
/** Get update specification */
function getUpdate(lines, begin, end) {
    const colonIdx = lines[begin].indexOf(CONTROL_SEP);
    if (colonIdx === -1) {
        throw new ModelError(`${scripting_tools_ERROR_MSG.MISS_COLON}. Correct the line: **${lines[begin]}**`, ERROR_LINK.UPDATE, lines[begin]);
    }
    const source = concatMultilineFormulas(lines.slice(begin + 1, end));
    const size = source.length;
    if (size < constants_UPDATE.MIN_LINES_COUNT)
        throw new ModelError(scripting_tools_ERROR_MSG.UPDATE_LINES_COUNT, ERROR_LINK.UPDATE, CONTROL_EXPR.UPDATE);
    if (source[constants_UPDATE.DURATION_IDX] == undefined)
        throw new ModelError(scripting_tools_ERROR_MSG.UPDATE_LINES_COUNT, ERROR_LINK.UPDATE, CONTROL_EXPR.UPDATE);
    const eqIdx = source[constants_UPDATE.DURATION_IDX].indexOf(EQUAL_SIGN);
    if (eqIdx === -1) {
        throw new ModelError(`Missing **${EQUAL_SIGN}**. Correct the line: **${source[constants_UPDATE.DURATION_IDX]}**`, ERROR_LINK.UPDATE, source[constants_UPDATE.DURATION_IDX]);
    }
    return {
        name: lines[begin].slice(colonIdx + 1),
        durationFormula: source[constants_UPDATE.DURATION_IDX].slice(eqIdx + 1).trim(),
        updates: source.slice(constants_UPDATE.DURATION_IDX + 1),
    };
} // getUpdate
/** Get output specification */
function getOutput(lines, begin, end) {
    const res = new Map();
    let line;
    let token;
    let eqIdx;
    let openIdx;
    let closeIdx;
    let colonIdx;
    for (let i = begin; i < end; ++i) {
        line = lines[i].trim();
        if (line.length < 1)
            continue;
        eqIdx = line.indexOf(EQUAL_SIGN);
        openIdx = line.indexOf(BRACE_OPEN);
        closeIdx = line.indexOf(BRACE_CLOSE);
        colonIdx = line.indexOf(CONTROL_SEP);
        // check braces
        if (openIdx * closeIdx <= 0) {
            throw new ModelError(`${scripting_tools_ERROR_MSG.BRACES} Correct the line **${line}** in the **${CONTROL_EXPR.OUTPUT}** block.`, ERROR_LINK.UI_OPTS, line);
        }
        // check ":"
        if (openIdx * colonIdx <= 0) {
            throw new ModelError(`${scripting_tools_ERROR_MSG.COLON} Correct the line **${line}** in the **#output** block.`, ERROR_LINK.UI_OPTS, line);
        }
        if (eqIdx === -1) { // no formula
            if (openIdx === -1) { // no annotation
                token = line.trim();
                res.set(token, {
                    caption: token,
                    formula: null,
                });
            }
            else { // there is an annotation
                token = line.slice(0, openIdx).trim();
                res.set(token, {
                    caption: line.slice(colonIdx + 1, closeIdx).trim(),
                    formula: null,
                });
            }
        }
        else // there is a formula
         if (openIdx === -1) { // no annotation
            token = line.slice(0, eqIdx).trim();
            res.set(token, {
                caption: token,
                formula: line.slice(eqIdx + 1),
            });
        }
        else { // there is an annotation
            token = line.slice(0, eqIdx).trim();
            res.set(token, {
                caption: line.slice(colonIdx + 1, closeIdx),
                formula: line.slice(eqIdx + 1, openIdx),
            });
        }
    } // for i
    return res;
} // getOutput
/** Get initial value problem specification given in the text */
function getIVP(text) {
    // The current Initial Value Problem (IVP) specification
    let name;
    let tags = null;
    let descr = null;
    let deqs;
    let exprs = null;
    let arg;
    let inits;
    let consts = null;
    let params = null;
    let tolerance = DEFAULT_TOL;
    let loop = null;
    const updates = [];
    const metas = [];
    let outputs = null;
    let solverSettings = DEFAULT_SOLVER_SETTINGS;
    let inputsLookup = null;
    // 0. Split text into lines & remove comments
    const lines = text.replaceAll('\t', ' ').split('\n')
        .map((s) => {
        const idx = s.indexOf(COMMENT_SEQ);
        return s.slice(0, (idx >= 0) ? idx : undefined).trimStart();
    })
        .filter((s) => s !== '');
    // 1. Skip first lines without the control tag
    const start = getStartOfProblemDef(lines);
    // 2. Get blocks limits
    const blocks = getBlocks(lines, start);
    // 3. Process blocks
    for (const block of blocks) {
        const firstLine = lines[block.begin];
        if (firstLine.startsWith(CONTROL_EXPR.NAME)) { // the 'name' block
            name = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
        }
        else if (firstLine.startsWith(CONTROL_EXPR.TAGS)) { // the 'tags' block
            tags = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
        }
        else if (firstLine.startsWith(CONTROL_EXPR.DESCR)) { // the 'description' block
            descr = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
        }
        else if (firstLine.startsWith(CONTROL_EXPR.DIF_EQ)) { // the 'differential equations' block
            deqs = getDifEquations(concatMultilineFormulas(lines.slice(block.begin + 1, block.end)));
        }
        else if (firstLine.startsWith(CONTROL_EXPR.EXPR)) { // the 'expressions' block
            exprs = getExpressions(concatMultilineFormulas(lines.slice(block.begin + 1, block.end)));
        }
        else if (firstLine.startsWith(CONTROL_EXPR.ARG)) { // the 'argument' block
            arg = getArg(concatMultilineFormulas(lines.slice(block.begin, block.end)));
        }
        else if (firstLine.startsWith(CONTROL_EXPR.INITS)) { // the 'initial values' block
            inits = getEqualities(lines, block.begin + 1, block.end);
        }
        else if (firstLine.startsWith(CONTROL_EXPR.CONSTS)) { // the 'constants' block
            consts = getEqualities(lines, block.begin + 1, block.end);
        }
        else if (firstLine.startsWith(CONTROL_EXPR.PARAMS)) { // the 'parameters' block
            params = getEqualities(lines, block.begin + 1, block.end);
        }
        else if (firstLine.startsWith(CONTROL_EXPR.TOL)) { // the 'tolerance' block
            tolerance = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
        }
        else if (firstLine.startsWith(CONTROL_EXPR.SOLVER)) { // the 'solver settings' block
            solverSettings = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
        }
        else if (firstLine.startsWith(CONTROL_EXPR.LOOP)) { // the 'loop' block
            loop = getLoop(lines, block.begin + 1, block.end);
        }
        else if (firstLine.startsWith(CONTROL_EXPR.UPDATE)) { // the 'update' block
            updates.push(getUpdate(lines, block.begin, block.end));
        }
        else if (firstLine.startsWith(CONTROL_EXPR.OUTPUT)) { // the 'output' block
            outputs = getOutput(lines, block.begin + 1, block.end);
        }
        else if (firstLine.startsWith(CONTROL_EXPR.INPUTS)) { // the 'inputs' block
            inputsLookup = firstLine.slice(firstLine.indexOf(CONTROL_SEP) + 1).trim();
        }
        else if (firstLine.startsWith(CONTROL_EXPR.COMMENT)) { // the 'comment' block
            // just skip it
        }
        else
            metas.push(firstLine.slice(CONTROL_TAG_LEN));
    } // for
    // check loop- & update-features: just one is supported
    if ((loop !== null) && (updates.length > 0))
        throw new ModelError(scripting_tools_ERROR_MSG.LOOP_VS_UPDATE, ERROR_LINK.LOOP_VS_UPDATE, CONTROL_EXPR.LOOP);
    // obtained model
    const ivp = {
        name: name,
        tags: tags,
        descr: descr,
        deqs: deqs,
        exprs: exprs,
        arg: arg,
        inits: inits,
        consts: consts,
        params: params,
        tolerance: tolerance,
        usedMathFuncs: getUsedMathIds(text, MATH_FUNCS),
        usedMathConsts: getUsedMathIds(text, MATH_CONSTS),
        loop: loop,
        updates: (updates.length === 0) ? null : updates,
        metas: metas,
        outputs: outputs,
        solverSettings: solverSettings,
        inputsLookup: inputsLookup,
    };
    checkCorrectness(ivp);
    return ivp;
} // getIVP
/** Return input specification, required for annotation generating */
function getInputSpec(inp) {
    if (inp.annot)
        return `${inp.value} ${inp.annot}`;
    return `${inp.value}`;
}
function getViewersSpec(ivp) {
    const outputColsCount = (ivp.outputs) ? ivp.outputs.size : ivp.inits.size;
    const multiAxis = (outputColsCount > MAX_LINE_CHART - 1) ? 'true' : 'false';
    const segments = (ivp.updates) ? ` segmentColumnName: "${STAGE_COL_NAME}",` : '';
    // eslint-disable-next-line max-len
    return `Line chart(block: 100, multiAxis: "${multiAxis}",${segments} multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | Grid(block: 100)`;
}
/** Return annotation line specifying viewers */
function getViewersLine(ivp) {
    const spec = getViewersSpec(ivp);
    return `viewer: ${spec}`;
}
/** Generate annotation lines */
function getAnnot(ivp, toAddViewers = true, toAddEditor = false) {
    const res = [];
    // the 'name' line
    res.push(`${ANNOT.NAME} ${ivp.name}`);
    // the 'tags' line
    if (ivp.tags)
        res.push(`${ANNOT.TAGS} ${ivp.tags}`);
    // the 'description' line
    if (ivp.descr)
        res.push(`${ANNOT.DESCR} ${ivp.descr}`);
    // the 'language' line
    res.push(ANNOT.LANG);
    // the 'loop' lines
    if (ivp.loop)
        res.push(`${ANNOT.INT_INPUT} ${LOOP.COUNT_NAME} = ${getInputSpec(ivp.loop.count)}`);
    // argument lines
    const arg = ivp.arg;
    const t0 = `${SERVICE}${arg.name}0`;
    const t1 = `${SERVICE}${arg.name}1`;
    const h = `${SERVICE}h`;
    res.push(`${ANNOT.DOUBLE_INPUT} ${t0} = ${getInputSpec(arg.initial)}`);
    res.push(`${ANNOT.DOUBLE_INPUT} ${t1} = ${getInputSpec(arg.final)}`);
    res.push(`${ANNOT.DOUBLE_INPUT} ${h} = ${getInputSpec(arg.step)}`);
    // initial values lines
    ivp.inits.forEach((val, key) => res.push(`${ANNOT.DOUBLE_INPUT} ${key} = ${getInputSpec(val)}`));
    // parameters lines
    if (ivp.params !== null)
        ivp.params.forEach((val, key) => res.push(`${ANNOT.DOUBLE_INPUT} ${key} = ${getInputSpec(val)}`));
    // the 'output' line
    if (toAddViewers)
        res.push(`${ANNOT.OUTPUT} {${ANNOT.CAPTION} ${ivp.name}; ${getViewersLine(ivp)}}`);
    else
        res.push(`${ANNOT.OUTPUT} {${ANNOT.CAPTION} ${ivp.name}`);
    // the 'editor' line
    if (toAddEditor) {
        res.push(ANNOT.EDITOR);
    }
    if (ivp.metas.length > 0)
        ivp.metas.forEach((line) => res.push(`//${line}`));
    else
        res.push(defaultMetas);
    return res;
} // getAnnot
/** Returns math functions arguments string */
function getMathArg(funcIdx) {
    return (funcIdx > POW_IDX) ? '(x)' : '(x, y)';
}
/** Return custom output lines: no expressions */
function getCustomOutputLinesNoExpressions(name, outputs, toAddUpdateCol) {
    const res = [''];
    res.push(SCRIPT.CUSTOM_OUTPUT_COM);
    res.push(SCRIPT.CUSTOM_COLUMNS);
    outputs.forEach((val, key) => {
        if (!val.formula)
            res.push(`${SCRIPT.SPACE2}${DF_NAME}.col('${key}'),`);
    });
    if (toAddUpdateCol)
        res.push(`${SCRIPT.SPACE2}${DF_NAME}.col('${STAGE_COL_NAME}'),`);
    res.push('];');
    outputs.forEach((val, key) => {
        if (!val.formula)
            res.push(`${DF_NAME}.col('${key}').name = '${val.caption}';`);
    });
    res.push(`${DF_NAME} = DG.DataFrame.fromColumns(${COLUMNS});`);
    res.push(`${DF_NAME}.name = '${name}';`);
    return res;
} // getCustomOutputLinesNoExpressions
/** Return custom output lines: with expressions */
function getCustomOutputLinesWithExpressions(ivp) {
    const res = [''];
    res.push(SCRIPT.CUSTOM_OUTPUT_COM);
    // 1. Solution raw data
    res.push(`const ${ivp.arg.name}RawData = ${DF_NAME}.col('${ivp.arg.name}').getRawData();`);
    res.push(`let ${ivp.arg.name};`);
    res.push(`const len = ${DF_NAME}.rowCount;\n`);
    ivp.inits.forEach((val, key) => {
        res.push(`const ${key}RawData = ${DF_NAME}.col('${key}').getRawData();`);
    });
    res.push('');
    // 2. Expressions raw data & variables
    ivp.exprs.forEach((val, key) => {
        if (ivp.outputs.has(key))
            res.push(`const ${key}RawData = new Float64Array(len);`);
        res.push(`let ${key};`);
    });
    res.push('');
    // 3. Computations
    res.push('for (let i = 0; i < len; ++i) {');
    res.push(`${SCRIPT.SPACE2}${ivp.arg.name} = ${ivp.arg.name}RawData[i];`);
    ivp.inits.forEach((val, key) => {
        res.push(`${SCRIPT.SPACE2}${key} = ${key}RawData[i];`);
    });
    ivp.exprs.forEach((val, key) => {
        res.push(`${SCRIPT.SPACE2}${key} = ${val};`);
        if (ivp.outputs.has(key))
            res.push(`${SCRIPT.SPACE2}${key}RawData[i] = ${key};`);
    });
    res.push('}\n');
    // 4. Form output
    res.push(`${DF_NAME} = DG.DataFrame.fromColumns([`);
    ivp.outputs.forEach((val, key) => {
        if (!val.formula)
            res.push(`${SCRIPT.SPACE2}DG.Column.fromFloat64Array('${val.caption}', ${key}RawData.slice(0, len)),`);
    });
    if (ivp.updates !== null)
        res.push(`${SCRIPT.SPACE2}${DF_NAME}.col('${STAGE_COL_NAME}'),`);
    res.push(']);');
    res.push(`${DF_NAME}.name = '${ivp.name}';`);
    return res;
} // getCustomOutputLinesWithExpressions
/** Return custom output lines */
function getCustomOutputLines(ivp) {
    var _a;
    if (ivp.exprs !== null) {
        for (const key of ivp.exprs.keys()) {
            if ((_a = ivp.outputs) === null || _a === void 0 ? void 0 : _a.has(key))
                return getCustomOutputLinesWithExpressions(ivp);
        }
    }
    return getCustomOutputLinesNoExpressions(ivp.name, ivp.outputs, (ivp.updates !== null));
} // getCustomOutputLinesWithExpressions
/** Return main body of JS-script: basic variant */
function getScriptMainBodyBasic(ivp) {
    const res = [];
    // 1. Constants lines
    if (ivp.consts !== null) {
        res.push('');
        res.push(SCRIPT.CONSTS);
        ivp.consts.forEach((val, key) => res.push(`const ${key} = ${val.value};`));
    }
    // 2. The problem definition lines
    res.push('');
    res.push(SCRIPT.ODE_COM);
    res.push(SCRIPT.ODE);
    res.push(`${SCRIPT.SPACE4}name: '${ivp.name}',`);
    // 2.1) argument
    const t = ivp.arg.name;
    const t0 = `${SERVICE}${t}0`;
    const t1 = `${SERVICE}${t}1`;
    const h = `${SERVICE}h`;
    res.push(`${SCRIPT.SPACE4}arg: {name: '${t}', start: ${t0}, finish: ${t1}, step: ${h}},`);
    const names = ivp.deqs.solutionNames;
    // 2.2) initial values
    res.push(`${SCRIPT.SPACE4}initial: [${names.join(', ')}],`);
    // 2.3) the right-hand side of the problem
    res.push(`${SCRIPT.SPACE4}func: (${t}, ${SERVICE}y, ${SERVICE}output) => {`);
    res.push(`${SCRIPT.SPACE6}${SCRIPT.FUNC_VALS}`);
    names.forEach((name, idx) => res.push(`${SCRIPT.SPACE6}const ${name} = ${SERVICE}y[${idx}];`));
    if (ivp.exprs !== null) {
        res.push(`\n${SCRIPT.SPACE6}${SCRIPT.EVAL_EXPR}`);
        ivp.exprs.forEach((val, key) => res.push(`${SCRIPT.SPACE6}const ${key} = ${val};`));
    }
    res.push(`\n${SCRIPT.SPACE6}${SCRIPT.COMP_OUT}`);
    names.forEach((name, idx) => res.push(`${SCRIPT.SPACE6}${SERVICE}output[${idx}] = ${ivp.deqs.equations.get(name)};`));
    res.push(`${SCRIPT.SPACE4}},`);
    // 2.4) final lines of the problem specification
    res.push(`${SCRIPT.SPACE4}tolerance: ${ivp.tolerance},`);
    //res.push(`${SCRIPT.SPACE6}solverOptions: ${ivp.solverSettings.replaceAll(ANNOT_SEPAR, COMMA)},`);
    res.push(`${SCRIPT.SPACE4}solutionColNames: [${names.map((key) => `'${key}'`).join(', ')}]`);
    res.push('};');
    // 2.5) solver options
    res.push('');
    res.push(`let opts = ${ivp.solverSettings.replaceAll(';', ',')};`);
    // 3. Math functions
    if (ivp.usedMathFuncs.length > 0) {
        res.push('');
        res.push(SCRIPT.MATH_FUNC_COM);
        ivp.usedMathFuncs.forEach((i) => res.push(`const ${MATH_FUNCS[i]} = ${getMathArg(i)} => Math.${MATH_FUNCS[i]}${getMathArg(i)};`));
    }
    // 4. Math constants
    if (ivp.usedMathConsts.length > 0) {
        res.push('');
        res.push(SCRIPT.MATH_CONST_COM);
        ivp.usedMathConsts.forEach((i) => res.push(`const ${MATH_CONSTS[i]} = Math.${MATH_CONSTS[i]};`));
    }
    // 5. The 'call solver' lines
    res.push('');
    res.push(SCRIPT.SOLVER_COM);
    res.push(SCRIPT.FUNCTION_CALL);
    // res.push(SCRIPT.SOLVER);
    // res.push(SCRIPT.PREPARE);
    // res.push(SCRIPT.CALL);
    // res.push(SCRIPT.OUTPUT);
    return res;
} // getScriptMainBodyBasic
/** Return function for JS-script */
function getScriptFunc(ivp, funcParamsNames) {
    const res = [];
    // 0. Function declaration
    res.push(SCRIPT.ONE_STAGE_COM);
    res.push(`${SCRIPT.ONE_STAGE_BEGIN}${funcParamsNames}${SCRIPT.ONE_STAGE_END}`);
    // 1. Constants lines
    if (ivp.consts !== null) {
        res.push('');
        res.push(SCRIPT.SPACE2 + SCRIPT.CONSTS);
        ivp.consts.forEach((val, key) => res.push(`${SCRIPT.SPACE2}const ${key} = ${val.value};`));
    }
    // 2. The problem definition lines
    res.push(SCRIPT.SPACE2 + SCRIPT.ODE_COM);
    res.push(SCRIPT.SPACE2 + SCRIPT.ODE);
    res.push(`${SCRIPT.SPACE6}name: '${ivp.name}',`);
    // 2.1) argument
    const t = ivp.arg.name;
    const t0 = `${SERVICE}${t}0`;
    const t1 = `${SERVICE}${t}1`;
    const h = `${SERVICE}h`;
    res.push(`${SCRIPT.SPACE6}arg: {name: '${t}', start: ${t0}, finish: ${t1}, step: ${h}},`);
    const names = ivp.deqs.solutionNames;
    // 2.2) initial values
    res.push(`${SCRIPT.SPACE6}initial: [${names.join(', ')}],`);
    // 2.3) the right-hand side of the problem
    res.push(`${SCRIPT.SPACE6}func: (${t}, ${SERVICE}y, ${SERVICE}output) => {`);
    res.push(`${SCRIPT.SPACE8}${SCRIPT.FUNC_VALS}`);
    names.forEach((name, idx) => res.push(`${SCRIPT.SPACE8}const ${name} = ${SERVICE}y[${idx}];`));
    if (ivp.exprs !== null) {
        res.push(`\n${SCRIPT.SPACE8}${SCRIPT.EVAL_EXPR}`);
        ivp.exprs.forEach((val, key) => res.push(`${SCRIPT.SPACE8}const ${key} = ${val};`));
    }
    res.push(`\n${SCRIPT.SPACE8}${SCRIPT.COMP_OUT}`);
    names.forEach((name, idx) => res.push(`${SCRIPT.SPACE8}${SERVICE}output[${idx}] = ${ivp.deqs.equations.get(name)};`));
    res.push(`${SCRIPT.SPACE6}},`);
    // 2.4) final lines of the problem specification
    res.push(`${SCRIPT.SPACE6}tolerance: ${ivp.tolerance},`);
    //res.push(`${SCRIPT.SPACE6}solverOptions: ${ivp.solverSettings.replaceAll(ANNOT_SEPAR, COMMA)},`);
    res.push(`${SCRIPT.SPACE6}solutionColNames: [${names.map((key) => `'${key}'`).join(', ')}]`);
    res.push(`${SCRIPT.SPACE2}};`);
    // 2.5) solver options
    res.push('');
    res.push(`${SCRIPT.SPACE2}let opts = ${ivp.solverSettings.replaceAll(';', ',')};`);
    // 3. Math functions
    if (ivp.usedMathFuncs.length > 0) {
        res.push('');
        res.push(SCRIPT.SPACE2 + SCRIPT.MATH_FUNC_COM);
        // eslint-disable-next-line max-len
        ivp.usedMathFuncs.forEach((i) => res.push(`${SCRIPT.SPACE2}const ${MATH_FUNCS[i]} = ${getMathArg(i)} => Math.${MATH_FUNCS[i]}${getMathArg(i)};`));
    }
    // 4. Math constants
    if (ivp.usedMathConsts.length > 0) {
        res.push('');
        res.push(SCRIPT.SPACE2 + SCRIPT.MATH_CONST_COM);
        ivp.usedMathConsts.forEach((i) => res.push(`${SCRIPT.SPACE2}const ${MATH_CONSTS[i]} = Math.${MATH_CONSTS[i]};`));
    }
    // 5. The 'call solver' lines
    res.push('');
    res.push(SCRIPT.SPACE2 + SCRIPT.SOLVER_COM);
    res.push(SCRIPT.SPACE2 + SCRIPT.FUNCTION_CALL);
    res.push(SCRIPT.SPACE2 + SCRIPT.RETURN);
    // res.push(SCRIPT.SPACE2 + SCRIPT.SOLVER);
    // res.push(SCRIPT.SPACE2 + SCRIPT.PREPARE);
    // res.push(SCRIPT.SPACE2 + SCRIPT.CALL);
    res.push(SCRIPT.SPACE2 + SCRIPT.RETURN_OUTPUT);
    // 6. Close the function
    res.push('};');
    return res;
} // getScriptFunc
/** Return main body of JS-script: loop case */
function getScriptMainBodyLoopCase(ivp) {
    const funcParamsNames = getFuncParamsNames(ivp);
    const res = getScriptFunc(ivp, funcParamsNames);
    res.push('');
    //res.push(`${SCRIPT.ASYNC_OUTPUT}${funcParamsNames});`);
    res.push(SCRIPT.SOLUTION_DF_COM);
    const dfNames = getSolutionDfColsNames(ivp);
    res.push(`let ${DF_NAME} = DG.DataFrame.fromColumns([`);
    dfNames.forEach((name) => res.push(`${SCRIPT.SPACE2}DG.Column.fromFloat64Array('${name}', []),`));
    res.push(`]);`);
    res.push(`${DF_NAME}.name = '${ivp.name}';`);
    res.push('');
    res.push(SCRIPT.LOOP_INTERVAL_COM);
    res.push(`const ${SCRIPT.LOOP_INTERVAL} = ${SERVICE}${ivp.arg.name}1 - ${SERVICE}${ivp.arg.name}0;`);
    res.push('');
    res.push(`let ${SCRIPT.LAST_IDX} = 0;\n`);
    res.push(SCRIPT.SOLVER_COM);
    res.push(`for (let ${SERVICE}idx = 0; ${SERVICE}idx < ${LOOP.COUNT_NAME}; ++${SERVICE}idx) {`);
    ivp.loop.updates.forEach((upd) => res.push(`${SCRIPT.SPACE2}${upd};`));
    res.push(`${SCRIPT.SPACE2}${SCRIPT.APPEND}${SCRIPT.ONE_STAGE}${funcParamsNames}), true);`);
    res.push(`${SCRIPT.SPACE2}${SERVICE}${ivp.arg.name}0 = ${SERVICE}${ivp.arg.name}1;`);
    res.push(`${SCRIPT.SPACE2}${SERVICE}${ivp.arg.name}1 += ${SCRIPT.LOOP_INTERVAL};`);
    res.push(`${SCRIPT.SPACE2}${SCRIPT.LAST_IDX} = ${DF_NAME}.rowCount - 1;`);
    dfNames.forEach((name, idx) => {
        if (idx !== 0)
            res.push(`${SCRIPT.SPACE2}${name} = ${DF_NAME}.get('${name}', ${SCRIPT.LAST_IDX});`);
    });
    // eslint-disable-next-line max-len
    res.push(`${SCRIPT.SPACE2}${DF_NAME}.set('${ivp.arg.name}', ${SCRIPT.LAST_IDX}, ${SERVICE}${ivp.arg.name}0 - Math.min(${SERVICE}h * ${STEP_RATIO}, ${TINY}));`);
    res.push('};');
    return res;
} // getScriptMainBodyLoopCase
/** Return main body of JS-script: update case */
function getScriptMainBodyUpdateCase(ivp) {
    var _a;
    const funcParamsNames = getFuncParamsNames(ivp);
    const res = getScriptFunc(ivp, funcParamsNames);
    res.push('');
    res.push(SCRIPT.SOLUTION_DF_COM);
    const dfNames = getSolutionDfColsNames(ivp);
    res.push(`${SCRIPT.ASYNC_OUTPUT}${funcParamsNames});`);
    res.push('');
    res.push(SCRIPT.SEGMENT_COM);
    // eslint-disable-next-line max-len
    res.push(`${DF_NAME}.columns.add(DG.Column.fromList('string', '${STAGE_COL_NAME}', new Array(${DF_NAME}.rowCount).fill('${(_a = ivp.arg.updateName) !== null && _a !== void 0 ? _a : INCEPTION}')));`);
    res.push('');
    res.push(`let ${SCRIPT.LAST_IDX} = 0;`);
    ivp.updates.forEach((upd, idx) => {
        res.push('');
        res.push(`${SCRIPT.UPDATE_COM} ${idx + 1}`);
        res.push(`const ${UPDATE.DURATION}${idx + 1} = ${upd.durationFormula};`);
        res.push(`${SCRIPT.LAST_IDX} = ${DF_NAME}.rowCount - 1;`);
        // eslint-disable-next-line max-len
        res.push(`${DF_NAME}.set('${ivp.arg.name}', ${SCRIPT.LAST_IDX}, ${SERVICE}${ivp.arg.name}1 - Math.min(${SERVICE}h * ${STEP_RATIO}, ${TINY}));`);
        dfNames.forEach((name, idx) => {
            if (idx !== 0)
                res.push(`${name} = ${DF_NAME}.get('${name}', ${SCRIPT.LAST_IDX});`);
        });
        upd.updates.forEach((upd) => res.push(`${upd};`));
        res.push(`${SERVICE}${ivp.arg.name}0 = ${SERVICE}${ivp.arg.name}1;`);
        res.push(`${SERVICE}${ivp.arg.name}1 += ${UPDATE.DURATION}${idx + 1};`);
        res.push(`const ${SERVICE}DF${idx + 1} = ${SCRIPT.ONE_STAGE}${funcParamsNames});`);
        // eslint-disable-next-line max-len
        res.push(`${SERVICE}DF${idx + 1}.columns.add(DG.Column.fromList('string', '${STAGE_COL_NAME}', new Array(${SERVICE}DF${idx + 1}.rowCount).fill('${upd.name}')));`);
        res.push(`${SCRIPT.APPEND}${SERVICE}DF${idx + 1}, true);`);
    });
    return res;
} // getScriptMainBodyUpdateCase
/** Return main body of JS-script */
function getScriptMainBody(ivp) {
    if (ivp.loop)
        return getScriptMainBodyLoopCase(ivp);
    if (ivp.updates)
        return getScriptMainBodyUpdateCase(ivp);
    return getScriptMainBodyBasic(ivp);
}
/** Return JS-script lines */
function getScriptLines(ivp, toAddViewers = true, toAddEditor = false) {
    const res = getAnnot(ivp, toAddViewers, toAddEditor).concat(getScriptMainBody(ivp));
    if (ivp.outputs)
        return res.concat(getCustomOutputLines(ivp));
    return res;
}
/** Return parameters of JS-script */
function getScriptParams(ivp) {
    const res = {};
    if (ivp.loop)
        res[LOOP.COUNT_NAME] = ivp.loop.count.value;
    const arg = ivp.arg;
    res[`${SERVICE}${arg.name}0`] = arg.initial.value;
    res[`${SERVICE}${arg.name}1`] = arg.final.value;
    res[`${SERVICE}h`] = arg.step.value;
    ivp.inits.forEach((val, key) => res[key] = val.value);
    if (ivp.params)
        ivp.params.forEach((val, key) => res[key] = val.value);
    return res;
}
/** Return func parameters names string */
function getFuncParamsNames(ivp) {
    const names = [];
    const arg = ivp.arg.name;
    names.push(`${SERVICE}${arg}0`);
    names.push(`${SERVICE}${arg}1`);
    names.push(`${SERVICE}h`);
    ivp.inits.forEach((val, key) => names.push(key));
    if (ivp.params)
        ivp.params.forEach((val, key) => names.push(key));
    return names.join(', ');
}
/** Return solution dataframe columns names*/
function getSolutionDfColsNames(ivp) {
    const res = [];
    res.push(ivp.arg.name);
    ivp.inits.forEach((val, key) => res.push(key));
    return res;
}
/** Check solver settings */
function checkSolverSettings(line) {
    const settings = new Map();
    const openBraceIdx = line.indexOf(BRACE_OPEN);
    const closeBraceIdx = line.indexOf(BRACE_CLOSE);
    let sepIdx;
    if ((openBraceIdx < 0) || (closeBraceIdx < 0))
        throw new ModelError(`${scripting_tools_ERROR_MSG.BRACES}. Correct the line **${line}**.`, ERROR_LINK.SOLVER_CONFIG, line);
    for (const item of line.slice(openBraceIdx + 1, closeBraceIdx).split(ANNOT_SEPAR)) {
        sepIdx = item.indexOf(CONTROL_SEP);
        if (sepIdx > 1)
            settings.set(item.slice(0, sepIdx).trim(), item.slice(sepIdx + 1).trim());
    }
    SOLVER_OPTIONS_RANGES.forEach((range, opt) => {
        if (settings.has(opt)) {
            const val = Number(settings.get(opt));
            if ((val < range.min) || (val > range.max)) {
                throw new ModelError(`${scripting_tools_ERROR_MSG.SOLVER}: **${opt}** must be in the range **${range.min}..${range.max}**.`, ERROR_LINK.SOLVER_CONFIG, line);
            }
        }
    });
} // checkSolverSettings
/** Check IVP correctness */
function checkCorrectness(ivp) {
    // 0. Check basic elements
    if (ivp.name === undefined)
        throw new ModelError(scripting_tools_ERROR_MSG.UNDEF_NAME, ERROR_LINK.CORE_BLOCKS);
    if ((ivp.deqs === undefined) || (ivp.deqs.equations.size === 0))
        throw new ModelError(scripting_tools_ERROR_MSG.UNDEF_DEQS, ERROR_LINK.CORE_BLOCKS);
    if ((ivp.inits === undefined) || (ivp.inits.size === 0))
        throw new ModelError(scripting_tools_ERROR_MSG.UNDEF_INITS, ERROR_LINK.CORE_BLOCKS);
    if (ivp.arg === undefined)
        throw new ModelError(scripting_tools_ERROR_MSG.UNDEF_ARG, ERROR_LINK.CORE_BLOCKS);
    // 1. Check initial values
    ivp.deqs.equations.forEach((ignore, name) => {
        if (!ivp.inits.has(name)) {
            throw new ModelError(`Initial value for **${name}** is missing. ${scripting_tools_ERROR_MSG.MISSING_INIT}`, ERROR_LINK.CORE_BLOCKS, CONTROL_EXPR.INITS);
        }
    });
    // 2. Check names of output columns
    const usedNames = [];
    let lowCase;
    if (ivp.outputs !== null) {
        ivp.outputs.forEach((val) => {
            lowCase = val.caption.toLowerCase();
            if (usedNames.includes(lowCase))
                throw new ModelError(`${scripting_tools_ERROR_MSG.CASE_INSENS}**${val.caption}**.`, ERROR_LINK.UI_OPTS, CONTROL_EXPR.OUTPUT);
            else
                usedNames.push(lowCase);
        });
    }
    else {
        const usedNames = [ivp.arg.name];
        ivp.deqs.solutionNames.forEach((name) => {
            lowCase = name.toLowerCase();
            if (usedNames.includes(lowCase))
                throw new ModelError(`${scripting_tools_ERROR_MSG.CASE_INSENS}**${name}**.`, ERROR_LINK.UI_OPTS);
            else
                usedNames.push(lowCase);
        });
    }
    // 3. Check argument
    if (ivp.arg.initial.value >= ivp.arg.final.value) {
        throw new ModelError(`${scripting_tools_ERROR_MSG.INTERVAL} **${ivp.arg.name}**. ${scripting_tools_ERROR_MSG.CORRECT_ARG_LIM}`, ERROR_LINK.CORE_BLOCKS, CONTROL_EXPR.ARG);
    }
    if (ivp.arg.step.value <= 0) {
        throw new ModelError(`${scripting_tools_ERROR_MSG.NEGATIVE_STEP}`, ERROR_LINK.CORE_BLOCKS);
    }
    if (ivp.arg.step.value > ivp.arg.final.value - ivp.arg.initial.value)
        throw new ModelError(scripting_tools_ERROR_MSG.INCOR_STEP, ERROR_LINK.CORE_BLOCKS);
    // 4. Check script inputs, due to https://reddata.atlassian.net/browse/GROK-15152
    const scriptInputs = [];
    let current;
    // 4.1) initial values
    ivp.inits.forEach((_, key) => {
        if (key[0] === SERVICE) {
            throw new ModelError(`${scripting_tools_ERROR_MSG.SERVICE_START} Correct **"${key}"** in the **${CONTROL_EXPR.INITS}** block.`, ERROR_LINK.CORE_BLOCKS);
        }
        current = key.toLocaleLowerCase();
        if (scriptInputs.includes(current)) {
            throw new ModelError(`${scripting_tools_ERROR_MSG.REUSE_NAME} **"${key}"** in the **${CONTROL_EXPR.INITS}** block.`, ERROR_LINK.CORE_BLOCKS);
        }
        scriptInputs.push(current);
    });
    // 4.2) parameters
    if (ivp.params !== null) {
        ivp.params.forEach((_, key) => {
            if (key[0] === SERVICE) {
                throw new ModelError(`${scripting_tools_ERROR_MSG.SERVICE_START}: correct **"${key}"** in the **${CONTROL_EXPR.PARAMS}** block.`, ERROR_LINK.MODEL_PARAMS);
            }
            current = key.toLocaleLowerCase();
            if (scriptInputs.includes(current)) {
                throw new ModelError(`${scripting_tools_ERROR_MSG.REUSE_NAME} **"${key}"** in the **${CONTROL_EXPR.PARAMS}** block.`, ERROR_LINK.MODEL_PARAMS);
            }
            scriptInputs.push(current);
        });
    }
    // 5. Check solver settings
    checkSolverSettings(ivp.solverSettings);
} // checkCorrectness

;// ./src/preparser.ts



const ivp = getIVP(code);
const lookupsOptions = ivp.inputsLookup ? [getLookupsInfoData(ivp.inputsLookup)] : [];
const argOptions = ARG_INPUT_KEYS.map((key) => getOptions(key, ivp.arg[key], CONTROL_EXPR.ARG));
const initsOptions = [...ivp.inits.entries()].map(([key, val]) => getOptions(key, val, CONTROL_EXPR.INITS));
const paramsOptions = ivp.params ? [...ivp.params.entries()].map(([key, val]) => getOptions(key, val, CONTROL_EXPR.PARAMS)) : [];
const loopOptions = ivp.loop ? [getOptions(constants_LOOP.COUNT_NAME, ivp.loop.count, CONTROL_EXPR.LOOP)] : [];
const inputs = ([...lookupsOptions, ...argOptions, ...initsOptions, ...paramsOptions, ...loopOptions]);
const outputs = [{
        name: 'df',
        type: 'dataframe',
        options: {
            caption: ivp.name,
            viewer: getViewersSpec(ivp),
        }
    }];
result = JSON.stringify({ inputs, outputs });

