""" create.py
    This script parses a file with initial value problem (IVP)
    and creates webautosolver Datagrok function.
"""

import os
import json
import shutil
import re

from createConstants import *
from settingsConstants import *
import export

def loadSettings():
    """
    Load export settings.
    """
    with open(EXPORT_SETTINGS_FILE, READ_MODE) as file:
        return json.load(file)

def getRawIVPfromFile(ivProblemFileName):
    """"
    Returns raw specification of IVP extracted from the file.
    """
    # dictionary with raw specifiaction of IVP
    ivpRaw = {NAME: None, METHOD: None, DESCRIPTION: None, TAGS: None,
              DIF_EQUATIONS: [], EXPRESSIONS: [],
              ARGUMENT: [], INIT_VALS: [], CONSTANTS: [], PARAMS: [] }

    # process file with IVP
    with open(ivProblemFileName, READ_MODE) as file:
        listOfLines = file.readlines()

        category = None

        # process each line of the IVP file
        for line in listOfLines:
            if line.startswith(CONTROL_TAG):
                category, ignore, rest = line.partition(SEPARATOR) # extract category data
                category = category.lstrip(CONTROL_TAG).strip() # extract name of category

                # process specific categories
                if category in {NAME, METHOD, DESCRIPTION, TAGS}:
                    ivpRaw[category] = rest.strip()
                elif category == ARGUMENT:
                    ivpRaw[ARGUMENT].append(rest.strip())

            else:
                line = line.strip()
                if(line != ''):
                    ivpRaw[category].append(line)

    return ivpRaw

def getSpecificationFromRawIVP(ivProblemRaw):
    """
    Processes raw description of IVP and return its specification.
    """

    def extractArgSpecification(argData):
        """
        Processes raw description of argument and returns it specification.
        """
        argSpecification = {}

        # 1) get name of argument & its title
        argName, ignore, argTitle = argData[0].partition(IVP_SPEC_SPACE) # 0-th element contains arg's data
        argName = argName.strip()
        argTitle = argTitle.strip(IVP_SPEC_STRIP_SYMBS)
        argSpecification[ARG_SPEC_NAME] = argName
        argSpecification[ARG_SPEC_TITLE] = argTitle

        # 2) get initial value data
        name, ignore, value = argData[1].partition(EQ_SIGN) # 1-st element contains arg's initial data
        name = name.strip()
        value = value.strip()
        argSpecification[ARG_SPEC_INITIAL] = {ARG_SPEC_NAME: name, ARG_SPEC_VALUE: value}

        # 3) get final value data
        name, ignore, value = argData[2].partition(EQ_SIGN) # 2-nd element contains arg's final data
        name = name.strip()
        value = value.strip()
        argSpecification[ARG_SPEC_FINAL] = {ARG_SPEC_NAME: name, ARG_SPEC_VALUE: value}

        # 4) get step data
        name, ignore, value = argData[3].partition(EQ_SIGN) # 3-rd element contains arg's step data
        name = name.strip()
        value = value.strip()
        argSpecification[ARG_SPEC_STEP] = {ARG_SPEC_NAME: name, ARG_SPEC_VALUE: value}

        return argSpecification

    def extractQuantitiesSpecification(rawDescription):
        """
        Processes quantities raw description and returns their specification.
        """
        specification = {}

        for line in rawDescription:

            # extract name
            name, ignore, rest = line.partition(EQ_SIGN)
            name = name.strip()

            # extract value and units
            value, ignore, units = rest.partition(IVP_SPEC_UNITS_OPEN)
            value = value.strip()
            units = units.strip(IVP_SPEC_UNITS_CLOSE)

            # store data
            specification[name] = {QUAN_SPEC_VALUE: value, QUAN_SPEC_UNITS: units}

        return specification

    def joinMultilineFormulas(formulas):
        """
        Returns joined multiline formulas.
        """
        # nothing to do if there are no formulas
        if len(formulas) == 0:
            return formulas

        joinedFormulas = []

        currentFormula = formulas[0]

        # each line without the symbol '=' is appended to the previous one
        for i in range(1, len(formulas)):
            if IVP_SPEC_EQUAL not in formulas[i]:
                currentFormula += IVP_SPEC_SPACE + formulas[i]
            else:
                joinedFormulas.append(currentFormula)
                currentFormula = formulas[i]

        joinedFormulas.append(currentFormula)

        return joinedFormulas

    def extractExpressionsSpecification(expressions):
        """
        Process raw expressions lines and returns their specification.
        """
        specification = {}

        for line in expressions:
            name, ignore, formula = line.partition(EQ_SIGN)
            specification[name.strip()] = formula.strip() # spaces are removed

        return specification

    def extractDifEqsSpecification(expressions):
        """
        Process raw differential equations' lines and returns their specification.
        """
        specification = {}

        for line in expressions:

            # get left & right parts
            left, ignore, right = line.partition(EQ_SIGN)

            # get expression with the name
            numerator, ignore, denominator = left.partition(FRAC_SEP)

            # get name: from 'dy' we extract 'y'
            ignore1, ignore, name = numerator.partition(DERIV_SYMBOL)

            # store results: spaces and parenthesis are removed
            specification[name.strip(IVP_SPEC_STRIP_SYMBS)] = right.strip()

        return specification

    ivProblem = {NAME: ivProblemRaw[NAME],
                 METHOD: ivProblemRaw[METHOD],
                 DESCRIPTION: ivProblemRaw[DESCRIPTION],
                 TAGS: ivProblemRaw[TAGS]}

    # 1. Process argument's description
    ivProblem[ARGUMENT] = extractArgSpecification(ivProblemRaw[ARGUMENT])

    # 2. Process initial values' description
    ivProblem[INIT_VALS] = extractQuantitiesSpecification(ivProblemRaw[INIT_VALS])

    # 3. Process constants' description
    ivProblem[CONSTANTS] = extractQuantitiesSpecification(ivProblemRaw[CONSTANTS])

    # 4. Process parameters' description
    ivProblem[PARAMS] = extractQuantitiesSpecification(ivProblemRaw[PARAMS])

    # 5. Process expressions' description
    ivProblem[EXPRESSIONS] = extractExpressionsSpecification(joinMultilineFormulas(ivProblemRaw[EXPRESSIONS]))

    # 6. Process differential equations' description
    ivProblem[DIF_EQUATIONS] = extractDifEqsSpecification(joinMultilineFormulas(ivProblemRaw[DIF_EQUATIONS]))


    # 7. Check results TODO: this section should be improved!

    # 7.1) check differential variables
    if set(ivProblem[DIF_EQUATIONS].keys()) != set(ivProblem[INIT_VALS].keys()):
        raise Warning

    return ivProblem

def generateCppCode(ivProblem):
    """
    Generates C++-code that solves the problem.
    """

    def replaceIntsByDoubles(expression):
        """
        Replace each integer value by the corresponding double number.
        """
        # add '.0' after each number
        bufferLine1 = re.sub(r'([^.A-Za-z_0-9]\d+)', r'\1.0', f' {expression}')

        # put spaces before/after artihmetic operations symbols
        bufferLine2 = re.sub(r'(\+|\*|-|/)', r' \1 ', bufferLine1)

        # remove extra '.0.'-s that may occur and double spaces
        result = bufferLine2.replace('.0.', '.').replace('  ', ' ')

        return result

    def putExternBlock(put, ivProblem):
        """
        Put the 'extern "C"' block required for Emscripten compiling.
        """
        put('\nextern "C" {\n')

        # name of solving function
        put(f'{SPACE} int {SOLVER_PREFIX}{ivProblem[NAME]}(')

        # arguments
        put(f'{DATA_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_INITIAL][ARG_SPEC_NAME]},')
        put(f' {DATA_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_FINAL][ARG_SPEC_NAME]},')
        put(f' {DATA_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_STEP][ARG_SPEC_NAME]},\n')

        # initial values
        line = SUBSPACE

        for name in ivProblem[INIT_VALS].keys():
            line += f'{DATA_TYPE} _{name}Initial, '

        put(line + '\n')

        # parameters
        line = SUBSPACE

        for name in ivProblem[PARAMS].keys():
            line += f'{DATA_TYPE} _{name}Val, '

        put(line + '\n')

        put(f'{SUBSPACE}int _{ivProblem[ARGUMENT][ARG_SPEC_NAME]}Count, ')
        put('int _varsCount,\n')

        # dataframe
        put(f'{SUBSPACE}{DATA_TYPE} * {SOLVER_RESULT_NAME}, ')
        put(f'int {SOLVER_RESULT_ROW_COUNT_NAME}, ')
        put(f'int {SOLVER_RESULT_COL_COUNT_NAME}) noexcept;\n')

        put('}\n')

    def putIVPcode(put, ivProblem):
        """
        Put C++-code specifying IVP.
        """
        # 1. PUT THE RIGHT-HAND SIDE OF ODEs

        put(f'\n\nnamespace {ivProblem[NAME]}\n')
        put('{\n')

        put(f'{SPACE}// tollerance\n')
        put(f'{SPACE}{DATA_TYPE} {TOLLERANCE_NAME} = {TOLLERANCE_VALUE};\n\n')

        put(f'{SPACE}// dimension of solution\n')
        put(f'{SPACE}const int {DIMENSION_NAME} = {len(ivProblem[INIT_VALS])};\n\n')

        put(f'{SPACE}// parameters\n')

        # put parameters declaration
        for name in ivProblem[PARAMS]:
            put(f'{SPACE}{ARG_TYPE} {name} = {PARAMETERS_INIT_VALUE};\n')

        # operating name of the argument, for example, t
        argName = ivProblem[ARGUMENT][ARG_SPEC_NAME]

        put(f'\n{SPACE}//the right part of the ODEs\n')

        # put ODEs right-hand side header, for example, f(...)
        put(SPACE)
        put(f'{VEC_TYPE} {ODES_RIGHT_PART_NAME}({ARG_TYPE} {argName}, ')
        put(f'{VEC_TYPE} & {VEC_ARG_NAME}) noexcept\n')

        # put body of the function
        put(SPACE + '{\n')

        put(f'{SUBSPACE}{VEC_TYPE} {VEC_OUTPUT_NAME}({DIMENSION_NAME});\n')

        if len(ivProblem[CONSTANTS]) > 0:
            put(f'\n{SUBSPACE}// constants\n')

        # put constants
        for name in ivProblem[CONSTANTS]:
            put(f'{SUBSPACE}const {ARG_TYPE} {name} = {ivProblem[CONSTANTS][name][QUAN_SPEC_VALUE]};\n')

        index = 0
        put(f'\n{SUBSPACE}// extract variables\n')

        # transform input vector coordinates to the current names
        for name in ivProblem[DIF_EQUATIONS]:
            put(f'{SUBSPACE}{ARG_TYPE} {name} = {VEC_ARG_NAME}({index});\n')
            index += 1

        if len(ivProblem[EXPRESSIONS]) > 0:
            put(f'\n{SUBSPACE}// expressions\n')

        # put expressions
        for name in ivProblem[EXPRESSIONS]:
            expression = replaceIntsByDoubles(ivProblem[EXPRESSIONS][name])
            put(f'{SUBSPACE}{ARG_TYPE} {name} = {expression};\n')
            #put(f'{SUBSPACE}{ARG_TYPE} {name} = {ivProblem[EXPRESSIONS][name]};\n')

        put(f'\n{SUBSPACE}// output computation\n')
        index = 0

        # computation lines
        for name in ivProblem[DIF_EQUATIONS]:
            expression = replaceIntsByDoubles(ivProblem[DIF_EQUATIONS][name])
            put(f'{SUBSPACE}{VEC_OUTPUT_NAME}({index}) = {expression};\n')
            #put(f'{SUBSPACE}{VEC_OUTPUT_NAME}({index}) = {ivProblem[DIF_EQUATIONS][name]};\n')
            index += 1

        put(f'\n{SUBSPACE}return {VEC_OUTPUT_NAME};\n')

        put(SPACE + '}' + f' // {ODES_RIGHT_PART_NAME}\n\n')

        # 1. PUT JACOBIAN

        put(f'{SPACE}// Jacobian (it is required, when applying implicit method)\n')

        # header
        put(f'{SPACE}{MATRIX_TYPE} {JACOBIAN_NAME}(')
        put(f'{ARG_TYPE} {argName}, ')
        put(f'{VEC_TYPE} & {VEC_ARG_NAME}, ')
        put(f'{ARG_TYPE} {EPS_NAME}) noexcept\n')

        # body
        put(SPACE + '{\n')

        put(f'{SUBSPACE}{MATRIX_TYPE} {MAT_OUTPUT_NAME}({DIMENSION_NAME}, {DIMENSION_NAME});\n')

        put(f'{SUBSPACE}{VEC_TYPE} {VAL_NAME} = {ODES_RIGHT_PART_NAME}')
        put(f'({argName}, {VEC_ARG_NAME});\n')

        put(f'{SUBSPACE}{VEC_TYPE} {Y_DER_NAME} = {VEC_ARG_NAME};\n\n')

        # put loop-block
        put(f'{SUBSPACE}for (int i = 0; i < {DIMENSION_NAME}; i++)' + ' {\n')
        put(f'{SUBSUBSPACE}{Y_DER_NAME}(i) += {EPS_NAME};\n')
        put(f'{SUBSUBSPACE}{MAT_OUTPUT_NAME}.col(i) = (')
        put(f'{ODES_RIGHT_PART_NAME}({argName}, {Y_DER_NAME}) - ')
        put(f'{VAL_NAME}) / {EPS_NAME};\n')
        put(f'{SUBSUBSPACE}{Y_DER_NAME}(i) -= {EPS_NAME};\n')

        put(SUBSPACE + '}\n\n')

        put(f'{SUBSPACE}return {MAT_OUTPUT_NAME};\n')

        put(SPACE + '}' + f' // {JACOBIAN_NAME}\n\n')

        # 3. PUT T-DERIVATIVE BLOCK

        put(f'{SPACE}// Derivative with respect to t (it is required, when applying implicit method)\n')

        # header
        put(f'{SPACE}{VEC_TYPE} {T_DERIVATIVE_NAME}(')
        put(f'{ARG_TYPE} {argName}, ')
        put(f'{VEC_TYPE} & {VEC_ARG_NAME}, ')
        put(f'{ARG_TYPE} {EPS_NAME}) noexcept\n')

        # body
        put(SPACE + '{\n')

        put(f'{SUBSPACE}return ({ODES_RIGHT_PART_NAME}(')
        put(f'{argName} + {EPS_NAME}, {VEC_ARG_NAME}) - ')
        put(f'{ODES_RIGHT_PART_NAME}({argName}, {VEC_ARG_NAME}))')
        put(f' / {EPS_NAME};\n')

        put(SPACE + '}' + f' // {T_DERIVATIVE_NAME}\n\n')

        put('}; ' + f'// {ivProblem[NAME]}\n')

    def putAnnotation(put, ivProblem):
        """
        Put ReachFunctionView annotation.
        """
        put(f'\n//name: {SOLVER_PREFIX + ivProblem[NAME]}\n')

        # 1. PUT ANNOTATION CONCERNING ARGUMENT

        # initial value
        put(f'//input: {ANNOT_ARG_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_INITIAL][ARG_SPEC_NAME]}')
        put(f' = {ivProblem[ARGUMENT][ARG_SPEC_INITIAL][ARG_SPEC_VALUE]}')
        put(' {caption: ' + ivProblem[ARGUMENT][ARG_SPEC_INITIAL][ARG_SPEC_NAME] + ';')
        put(f' category: {ivProblem[ARGUMENT][ARG_SPEC_TITLE]}' + '}\n')

        # final value
        put(f'//input: {ANNOT_ARG_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_FINAL][ARG_SPEC_NAME]}')
        put(f' = {ivProblem[ARGUMENT][ARG_SPEC_FINAL][ARG_SPEC_VALUE]}')
        put(' {caption: ' + ivProblem[ARGUMENT][ARG_SPEC_FINAL][ARG_SPEC_NAME] + ';')
        put(f' category: {ivProblem[ARGUMENT][ARG_SPEC_TITLE]}' + '}\n')

         # step value
        put(f'//input: {ANNOT_ARG_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_STEP][ARG_SPEC_NAME]}')
        put(f' = {ivProblem[ARGUMENT][ARG_SPEC_STEP][ARG_SPEC_VALUE]}')
        put(' {caption: ' + ivProblem[ARGUMENT][ARG_SPEC_STEP][ARG_SPEC_NAME] + ';')
        put(f' category: {ivProblem[ARGUMENT][ARG_SPEC_TITLE]}' + '}\n')

        # 2. INITIAL VALUES

        # dict with initial values specification
        data = ivProblem[INIT_VALS]

        # put annotation
        for name in data.keys():
            put(f'//input: {ANNOT_ARG_TYPE} _{name}Initial = {data[name][QUAN_SPEC_VALUE]}')
            put(' {' + f'units: {data[name][QUAN_SPEC_UNITS]};')
            put(f' caption: {name}; category: {INIT_VALS}' + '}\n')

        # 3. PARAMETERS

        # dict with parameters specification
        data = ivProblem[PARAMS]

        # put annotation
        for name in data.keys():
            put(f'//input: {ANNOT_ARG_TYPE} _{name}Val = {data[name][QUAN_SPEC_VALUE]}')
            put(' {' + f'units: {data[name][QUAN_SPEC_UNITS]};')
            put(f' caption: {name}; category: {PARAMS}' + '}\n')

        # 4. SPECIAL Input
        put(f'//input: int _{ivProblem[ARGUMENT][ARG_SPEC_NAME]}Count\n')
        put(f'//input: int _varsCount\n')

        # 5. OUTPUT COLUMNS
        put(f'//output: column_list {SOLVER_RESULT_NAME} ')
        put(f'[new(_{ivProblem[ARGUMENT][ARG_SPEC_NAME]}Count, _varsCount)')

        # create new columns names
        line = f"; '{ivProblem[ARGUMENT][ARG_SPEC_NAME]}'"

        for name in ivProblem[DIF_EQUATIONS]:
            line += f"; '{name}({ivProblem[ARGUMENT][ARG_SPEC_NAME]})'"

        put(line + ']\n')

        # 6. OUTPUT DATAFRAME
        put(f'//output: dataframe {DF_SOLUTION_NAME} [{SOLVER_RESULT_NAME}] {DF_OUTPUT_ANNOT}\n')

        # 7. EDITOR
        put(EDITOR_LINE + '\n')

    def putSolverHeader(put, ivProblem):
        """
        Put IVP solving function header.
        """
        put('EMSCRIPTEN_KEEPALIVE\n')

        # name of solving function
        put(f'int {SOLVER_PREFIX}{ivProblem[NAME]}(')

        # arguments
        put(f'{DATA_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_INITIAL][ARG_SPEC_NAME]},')
        put(f' {DATA_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_FINAL][ARG_SPEC_NAME]},')
        put(f' {DATA_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_STEP][ARG_SPEC_NAME]},\n')

        # initial values
        line = SPACE

        for name in ivProblem[INIT_VALS].keys():
            line += f'{DATA_TYPE} _{name}Initial, '

        put(line + '\n')

        # parameters
        line = SPACE

        for name in ivProblem[PARAMS].keys():
            line += f'{DATA_TYPE} _{name}Val, '

        put(line + '\n')

        put(f'{SPACE}int _{ivProblem[ARGUMENT][ARG_SPEC_NAME]}Count, ')
        put('int _varsCount,\n')

        # dataframe
        put(f'{SPACE}{DATA_TYPE} * {SOLVER_RESULT_NAME}, ')
        put(f'int {SOLVER_RESULT_ROW_COUNT_NAME}, ')
        put(f'int {SOLVER_RESULT_COL_COUNT_NAME}) noexcept\n')

    def putSolverBody(put, ivProblem):
        """
        Put IVP solving function body.
        """
        put('{\n')

        put(f'{SPACE}using namespace {ivProblem[NAME]};\n\n')

        put(f'{SPACE}{DATA_TYPE} {INIT_VALS_ARR_NAME}[DIM];\n\n')

        put(f'{SPACE}// initial values\n')

        index = 0

        # put initial values
        #for name in ivProblem[INIT_VALS]:
        for name in ivProblem[DIF_EQUATIONS]:
            put(f'{SPACE}{INIT_VALS_ARR_NAME}[{index}] = _{name}Initial;\n')
            index += 1

        put(f'\n{SPACE}// parameters\n')

        # put parameters
        for name in ivProblem[PARAMS]:
            put(f'{SPACE}{name} = _{name}Val;\n')

        # creating return-line
        line = f'\n{SPACE}return solveODE({ODES_RIGHT_PART_NAME}, '

        if ivProblem[METHOD] == IVP_IMPLICIT_METHOD:
            line += f'{T_DERIVATIVE_NAME}, {JACOBIAN_NAME}, '

        line += f'{ivProblem[ARGUMENT][ARG_SPEC_INITIAL][ARG_SPEC_NAME]}, '
        line += f'{ivProblem[ARGUMENT][ARG_SPEC_FINAL][ARG_SPEC_NAME]}, '
        line += f'{ivProblem[ARGUMENT][ARG_SPEC_STEP][ARG_SPEC_NAME]}, '
        line += f'{INIT_VALS_ARR_NAME}, {TOLLERANCE_NAME},\n' + SUBSUBSPACE
        line += f'{SOLVER_RESULT_NAME}, _{ivProblem[ARGUMENT][ARG_SPEC_NAME]}Count, _varsCount);\n'

        put(line)

        put('} //' + SOLVER_PREFIX + ivProblem[NAME])

    settings = {}

    # load settings
    with open(EXPORT_SETTINGS_FILE, READ_MODE) as file:
         settings = json.load(file)

    # generate C++-code
    with open(f'{settings[SET_FOLDER]}/{ivProblem[NAME]}.cpp', WRITE_MODE) as file:
        put = file.write

        put(f'// {ivProblem[NAME]}.cpp\n')

        # STL includes TODO: consider other libs
        put('\n#include <cmath>\n')
        put('using namespace std;\n')

        # Emscripten tools
        put(f'\n#include <{settings[SET_EMSCRIPTEN]}>\n')

        # the 'extern "C"' block - it's required for Emscripten
        putExternBlock(put, ivProblem)

        # Eigen tools
        put(f'\n#include "{settings[SET_EIGEN]}"\n')
        put('using namespace Eigen;\n')

        # ODEs solver tools
        put(f'\n#include "{ODES_SOLVER_LIB}"\n')
        put(f'using namespace {ODES_SOLVER_NAMESPACE};')

        # IVP code
        putIVPcode(put, ivProblem)

        # RichFunctionView annotation
        putAnnotation(put, ivProblem)

        # solving function's header
        putSolverHeader(put, ivProblem)

        # solving function's body
        putSolverBody(put, ivProblem)

def updateExportSettings(settings, ivProblem):
    """
    Update C++-to-wasm export settings.
    """
    settings[SET_NAME] = ivProblem[NAME]
    settings[SET_SOURCE] = [f'{ivProblem[NAME]}.cpp']

def cppToWasmExport(settings):
    """
    C++-to-Wasm export.
    """
    # parse C++-files and get exported functions data
    functionsData = export.getFunctionsFromListOfFiles(settings)

    # get command for the Emscripten tool
    command = export.getCommand(settings, functionsData)

    # execute command by Emscripten
    os.system(command)

    # create JS-wrapper for wasm-functions run in webworkers
    export.createJSwrapperForWasmInWebWorkerRun(settings)

    # create worker files
    export.createWorkerFiles(settings, functionsData)

    # complete JS-file created by Emscripten
    export.completeJsWasmfile(settings, functionsData)

    # update the file package.json
    export.updatePackageJsonFile(settings)

def addSolverToPackageFile(settings, ivProblem):
    """
    Add IVP solver function to DATAGROK package file.
    """

    def putAnnotation(put, ivProblem, inMainStreamComutation=True):
        """
        Put RichFunctionView annotation.
        """
        # 0. PUT FUNCTION NAME ANNOTATION
        put(f'\n//name: {SOLVER_PREFIX + ivProblem[NAME]}')
        put('\n' if inMainStreamComutation else f'{IN_WEBWORKER_SUFFIX}\n')

        # put description line
        if ivProblem[DESCRIPTION] is not None:
            put(f'//description: {ivProblem[DESCRIPTION]}\n')

        # put tags
        if ivProblem[TAGS] is not None and inMainStreamComutation:
            put(f'//tags: {ivProblem[TAGS]}\n')

        # 1. PUT ANNOTATION CONCERNING ARGUMENTS

        # initial value
        put(f'//input: {ANNOT_ARG_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_INITIAL][ARG_SPEC_NAME]}')
        put(f' = {ivProblem[ARGUMENT][ARG_SPEC_INITIAL][ARG_SPEC_VALUE]}')
        put(' {caption: ' + ivProblem[ARGUMENT][ARG_SPEC_INITIAL][ARG_SPEC_NAME] + ';')
        put(f' category: {ivProblem[ARGUMENT][ARG_SPEC_TITLE]}' + '}\n')

        # final value
        put(f'//input: {ANNOT_ARG_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_FINAL][ARG_SPEC_NAME]}')
        put(f' = {ivProblem[ARGUMENT][ARG_SPEC_FINAL][ARG_SPEC_VALUE]}')
        put(' {caption: ' + ivProblem[ARGUMENT][ARG_SPEC_FINAL][ARG_SPEC_NAME] + ';')
        put(f' category: {ivProblem[ARGUMENT][ARG_SPEC_TITLE]}' + '}\n')

         # step value
        put(f'//input: {ANNOT_ARG_TYPE} {ivProblem[ARGUMENT][ARG_SPEC_STEP][ARG_SPEC_NAME]}')
        put(f' = {ivProblem[ARGUMENT][ARG_SPEC_STEP][ARG_SPEC_VALUE]}')
        put(' {caption: ' + ivProblem[ARGUMENT][ARG_SPEC_STEP][ARG_SPEC_NAME] + ';')
        put(f' category: {ivProblem[ARGUMENT][ARG_SPEC_TITLE]}' + '}\n')

        # 2. INITIAL VALUES

        # dict with initial values specification
        data = ivProblem[INIT_VALS]

        # put annotation
        for name in data.keys():
            put(f'//input: {ANNOT_ARG_TYPE} _{name}Initial = {data[name][QUAN_SPEC_VALUE]}')
            put(' {' + f'units: {data[name][QUAN_SPEC_UNITS]};')
            put(f' caption: {name}; category: {INIT_VALS}' + '}\n')

        # 3. PARAMETERS

        # dict with parameters specification
        data = ivProblem[PARAMS]

        # put annotation
        for name in data.keys():
            put(f'//input: {ANNOT_ARG_TYPE} _{name}Val = {data[name][QUAN_SPEC_VALUE]}')
            put(' {' + f'units: {data[name][QUAN_SPEC_UNITS]};')
            put(f' caption: {name}; category: {PARAMS}' + '}\n')

        # 4. OUTPUT DATAFRAME & EDITOR if running computations in the main stream
        if inMainStreamComutation:
            put(f'//output: dataframe {DF_SOLUTION_NAME} {DF_OUTPUT_ANNOT}\n')
            put(EDITOR_LINE + '\n')

    def putHeader(put, ivProblem, inMainStreamComutation=True):
        """
        Put header of the function.
        """
        # name
        put(f'export async function {SOLVER_PREFIX + ivProblem[NAME]}')
        put('(' if inMainStreamComutation else f'{IN_WEBWORKER_SUFFIX}(')

        # argument
        arg = ivProblem[ARGUMENT]
        put(f'{arg[ARG_SPEC_INITIAL][ARG_SPEC_NAME]}, ')
        put(f'{arg[ARG_SPEC_FINAL][ARG_SPEC_NAME]}, ')
        put(f'{arg[ARG_SPEC_STEP][ARG_SPEC_NAME]},\n')

        # initial values
        line = ''
        for name in ivProblem[INIT_VALS].keys():
            line += f'_{name}Initial' + ', '
        put(f'{JS_SPACE}{line}\n')

        # parameters
        line = ''
        for name in ivProblem[PARAMS].keys():
            line += f'_{name}Val' + ', '

        # remove last comma
        line = line.rstrip(', ')
        put(f'{JS_SPACE}{line})\n')

    def putBodyMainStreamComputations(put, ivProblem):
        """
        Put body of the function for computations in the main stream.
        """
        put('{\n')

        # put computation of dataframe size: rows count
        arg = ivProblem[ARGUMENT]
        put(f'{JS_SPACE}let _{ivProblem[ARGUMENT][ARG_SPEC_NAME]}Count = ')
        put(f'Math.trunc(({arg[ARG_SPEC_FINAL][ARG_SPEC_NAME]} ')
        put(f'- {arg[ARG_SPEC_INITIAL][ARG_SPEC_NAME]}) ')
        put(f'/ {arg[ARG_SPEC_STEP][ARG_SPEC_NAME]}) + 1;\n')

        # put computation of dataframe size: columns count
        put(f'{JS_SPACE}let _varsCount = ')
        put(f'{len(ivProblem[INIT_VALS].keys()) + 1};\n\n')

        # put return line
        put(f"{JS_SPACE}return callWasm({ivProblem[NAME]}, '{SOLVER_PREFIX + ivProblem[NAME]}',\n")

        # put argument
        put(f'{JS_SUBSPACE}[')
        arg = ivProblem[ARGUMENT]
        put(f' {arg[ARG_SPEC_INITIAL][ARG_SPEC_NAME]}, ')
        put(f'{arg[ARG_SPEC_FINAL][ARG_SPEC_NAME]}, ')
        put(f'{arg[ARG_SPEC_STEP][ARG_SPEC_NAME]},\n')

        # put initial values
        line = ''
        for name in ivProblem[INIT_VALS].keys():
            line += f' _{name}Initial,'
        put(JS_SUBSPACE + line + '\n')

        # put parameters
        line = ''
        for name in ivProblem[PARAMS].keys():
            line += f' _{name}Val,'
        put(JS_SUBSPACE + line + '\n')

        # put row- and col- counts
        put(f'{JS_SUBSPACE} _{ivProblem[ARGUMENT][ARG_SPEC_NAME]}Count, _varsCount ] );\n')

        put('}\n')

    def putBodyWebworkerComputations(put, ivProblem):
        """
        Put body of the function for computations in the main stream.
        """
        put('{\n')

        # put computation of dataframe size: rows count
        arg = ivProblem[ARGUMENT]
        put(f'{JS_SPACE}let _{ivProblem[ARGUMENT][ARG_SPEC_NAME]}Count = ')
        put(f'Math.trunc(({arg[ARG_SPEC_FINAL][ARG_SPEC_NAME]} ')
        put(f'- {arg[ARG_SPEC_INITIAL][ARG_SPEC_NAME]}) ')
        put(f'/ {arg[ARG_SPEC_STEP][ARG_SPEC_NAME]}) + 1;\n')

        # put computation of dataframe size: columns count
        put(f'{JS_SPACE}let _varsCount = ')
        put(f'{len(ivProblem[INIT_VALS].keys()) + 1};\n\n')

        funcName = f'{SOLVER_PREFIX}{ivProblem[NAME]}'

        workerName = f'{settings[export.FOLDER]}/{funcName}{export.WORKER_SUFFIX}{export.WORKER_EXTENSION}'
        put(f"{JS_SPACE}var worker = new Worker(new URL('{workerName}', import.meta.url));\n\n")

        put(f"{JS_SPACE}worker.postMessage(getCppInput({settings[export.NAME]}['{funcName}'].arguments,\n")

        # put argument
        put(f'{JS_SUBSPACE}[')
        arg = ivProblem[ARGUMENT]
        put(f' {arg[ARG_SPEC_INITIAL][ARG_SPEC_NAME]}, ')
        put(f'{arg[ARG_SPEC_FINAL][ARG_SPEC_NAME]}, ')
        put(f'{arg[ARG_SPEC_STEP][ARG_SPEC_NAME]},\n')

        # put initial values
        line = ''
        for name in ivProblem[INIT_VALS].keys():
            line += f' _{name}Initial,'
        put(JS_SUBSPACE + line + '\n')

        # put parameters
        line = ''
        for name in ivProblem[PARAMS].keys():
            line += f' _{name}Val,'
        put(JS_SUBSPACE + line + '\n')

        # put row- and col- counts
        put(f'{JS_SUBSPACE} _{ivProblem[ARGUMENT][ARG_SPEC_NAME]}Count, _varsCount ] ));\n\n')

        put(JS_SPACE + 'worker.onmessage = function(e) {\n')

        put(f"{JS_SUBSPACE}let output = getResult({settings[export.NAME]}['{funcName}'], e.data);\n")

        put(f'\n{CUSTOM_UI_LINES}')

        put(JS_SPACE + '}\n')

        put('}\n\n')


    # set mode for openning file
    openMode = WRITE_MODE if settings[SET_PCKG_FILE_UPD_MODE] == REWRITE else APPEND_MODE

    # put code to the package file
    with open(settings[SET_PACKAGE_FILE], openMode) as file:
        put = file.write

        # put first lines if package file is opened in the write mode
        if openMode == WRITE_MODE:
            put(PACKAGE_FILE_FIRST_LINES)

        # put init-function
        put('//tags: init\nexport async function init() {\n')
        put(f'{JS_SPACE}await init{ivProblem[NAME]}();\n' + '}\n')

        # 1. Generation of code for solving IVP in the main stream
        putAnnotation(put, ivProblem)
        putHeader(put, ivProblem)
        putBodyMainStreamComputations(put, ivProblem)

        # 2. Generation of code for solving IVP in webworker
        putAnnotation(put, ivProblem, False)
        putHeader(put, ivProblem, False)
        putBodyWebworkerComputations(put, ivProblem)

def main():
    """
    The main script:
        1) load settings;
        2) get initial value problem from file;
        3) get specification of IVP from its raw description;
        4) generate C++ code;
        5) update settings for C++-to-wasm export;
        6) C++-to-wasm export;
        7) add IVP solver function to DATAGROK package file;
        8) delete of the __pycache__ folder.
    """
    try:
        # 1) load settings
        settings = loadSettings()

        # 2) get initial value problem from file
        ivProblemRaw = getRawIVPfromFile(settings[SET_IVP_FILE])

        # save raw IVP data, it's useful when debugging
        #with open('rawIVP.json', WRITE_MODE) as file:
        #    json.dump(ivProblemRaw, file)

        # 3) get specification of IVP from its raw description
        ivProblem = getSpecificationFromRawIVP(ivProblemRaw)

        # save parsed IVP specification, it's useful when debugging
        #with open('IVP.json', WRITE_MODE) as file:
        #    json.dump(ivProblem, file)

        # 4) generate C++ code
        generateCppCode(ivProblem)

        # 5) update settings for C++-to-wasm export
        updateExportSettings(settings, ivProblem)

        # 6) C++-to-wasm export
        cppToWasmExport(settings)

        # 7) add IVP solver function to DATAGROK package file
        addSolverToPackageFile(settings, ivProblem)

        # 8) delete of the __pycache__ folder, otherwise it must be done further manually
        shutil.rmtree('__pycache__')

    except OSError:
        print("OS ERROR: check paths and file names!")

    except KeyError:
        print("PARSING ERROR: check specificaition of the initial value problem!")

    except IndexError:
        print("PARSING ERROR: check C/C++-code!")

    except Warning:
        print("IVP SPECIFICATION ERROR: check names ODEs and initial values!")

if __name__ == '__main__':
    main()
