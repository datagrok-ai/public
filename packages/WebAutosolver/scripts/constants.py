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

# tag that defines a start of category specification
CONTROL_TAG = '#'

# name of file with C++-to-wasm settings
EXPORT_SETTINGS_FILE = 'exportSettings.json'

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

