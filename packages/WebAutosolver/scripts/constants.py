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

# Types
DATA_TYPE = 'float' # dataframe type
ARG_TYPE = 'double' # solver operating type
VEC_TYPE = 'VectorXd' # sollver operating vector type 

SOLVER_PREFIX = 'solve'
