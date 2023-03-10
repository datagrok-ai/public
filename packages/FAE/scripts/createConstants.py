""" constants.py
    Definition of constants required for parsing initial value problem (IVP)
    specification file.
"""

# Names of special categories specifying IVP
NAME = 'name'
METHOD = 'method'
DIF_EQUATIONS = 'differential equations'
EXPRESSIONS = 'expressions'
ARGUMENT = 'argument'
INIT_VALS = 'initial values'
CONSTANTS = 'constants'
PARAMS = 'parameters'
DESCRIPTION = 'description'
TAGS = 'tag'

# Solving methods constants
IVP_IMPLICIT_METHOD = 'implicit'
IVP_EXPLICIT_METHOD = 'explicit'

# IVP specification constants
IVP_SPEC_SPACE = ' '
IVP_SPEC_UNITS_OPEN = '('
IVP_SPEC_UNITS_CLOSE = ' )'
IVP_SPEC_STRIP_SYMBS = ' ()'
IVP_SPEC_EQUAL = '='

# tag that defines a start of category specification
CONTROL_TAG = '#'

# separator of category name and data
SEPARATOR = ':'

# separator between left & right parts
EQ_SIGN = '='

# Argumnent specification constants
ARG_SPEC_NAME = 'name'
ARG_SPEC_TITLE = 'title'
ARG_SPEC_INITIAL = 'initial'
ARG_SPEC_FINAL = 'final'
ARG_SPEC_STEP = 'step'
ARG_SPEC_VALUE = 'value'

# Quantities specification constants
QUAN_SPEC_VALUE = 'value'
QUAN_SPEC_UNITS = 'units'

# Fraction separator: it separates numerator and denominator
FRAC_SEP = '/'

# Derivative symbol: 'd' for the case dy/dt
DERIV_SYMBOL = 'd'


# name of file with C++-to-wasm settings
EXPORT_SETTINGS_FILE = 'createSettings.json'

# ODEs solver lib
ODES_SOLVER_LIB = 'odes.h'
ODES_SOLVER_NAMESPACE = 'odes'

# C++ formatting
SPACE = ' ' * 4
SUBSPACE = ' ' * 8
SUBSUBSPACE = ' ' * 12

# Types
DATA_TYPE = 'float' # dataframe type
ARG_TYPE = 'double' # solver operating type
VEC_TYPE = 'VectorXd' # solver operating vector type
MATRIX_TYPE = 'MatrixXd' # solver operating matrix type

# ODEs solver objects
SOLVER_PREFIX = 'solve'
IN_WEBWORKER_SUFFIX = 'InWebWorker'
TOLLERANCE_NAME = '_tol'
TOLLERANCE_VALUE = '0.00005f'
DIMENSION_NAME = 'DIM'
PARAMETERS_INIT_VALUE = '0.0'
ODES_RIGHT_PART_NAME = '_f'
VEC_ARG_NAME = '_y'
VEC_OUTPUT_NAME = '_res'
JACOBIAN_NAME = '_J'
T_DERIVATIVE_NAME = '_T'
EPS_NAME = '_eps'
MAT_OUTPUT_NAME = '_res'
VAL_NAME = '_val'
Y_DER_NAME = '_yDer'
SOLVER_RESULT_NAME = '_solution'
SOLVER_RESULT_ROW_COUNT_NAME = '_solutionRowCount'
SOLVER_RESULT_COL_COUNT_NAME = '_solutionColCount'
INIT_VALS_ARR_NAME = '_initialVals'

# RichFunctionView annotation constants
ANNOT_ARG_TYPE = 'double'
EDITOR_LINE = '//editor: Compute:RichFunctionViewEditor'
DF_SOLUTION_NAME = 'dfSolution'
DF_OUTPUT_ANNOT = '{caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }'

# DATAGROK package file constants
JS_SPACE = ' ' * 2
JS_SUBSPACE = ' ' * 4

PACKAGE_FILE_FIRST_LINES = ''' // THIS FILE IS GENERATED AUTOMATICALLY. DO NOT CHANGE ANYTHING!

/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

// Imports for call wasm runtime-system: in the main stream and in webworkers
import { callWasm } from '../wasm/callWasm';
import { getCppInput, getResult } from '../wasm/callWasmForWebWorker';\n
'''

CUSTOM_UI_LINES = '''    let view = grok.shell.addTableView(output);
    view.lineChart({ markerType: 'dot', sharex: 'true', multiAxis: 'true'});
'''

# File operating constants
WRITE_MODE = 'w'
APPEND_MODE = 'a'
READ_MODE = 'r'
REWRITE = 'rewrite'
