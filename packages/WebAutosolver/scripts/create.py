""" create.py
    This script parses a file with initial value problem (IVP) 
    and creates webautosolver Datagrok function.
"""

import os
import json

from constants import *

def getRawIVPfromFile(ivpFileName):
    """"
    Returns raw specification of IVP extracted from the file.
    """
    # dictionary with raw specifiaction of IVP
    ivpRaw = {NAME: None, METHOD: None, DIF_EQUATIONS: [], EXPRESSIONS: [],
    ARGUMENT: [], INIT_VALS: [], CONSTANTS: [], PARAMS: [] }

    # process file with IVP 
    with open(ivpFileName, 'r') as file:
        listOfLines = file.readlines()

        category = None

        # process each line of the IVP file
        for line in listOfLines:
            if line.startswith(CONTROL_TAG):                
                category, ignore, rest = line.partition(':') # extract category data
                category = category.lstrip('#').strip() # extract name of category

                # process specific categories
                if category in {NAME, METHOD}:                
                    ivpRaw[category] = rest.strip()
                elif category == ARGUMENT:
                    ivpRaw[ARGUMENT].append(rest.strip())                        

            else:
                line = line.strip()           
                if(line != ''):
                    ivpRaw[category].append(line)        

    return ivpRaw

def getSpecificationFromRawIVP(ivpRaw):
    """
    Processes raw description of IVP and return its specification.
    """

    def extractArgSpecification(argData):
        """ 
        Processes raw description of argument and returns it specification.
        """

        argSpecification = {}

        # 1) get name of argument & its title
        argName, ignore, argTitle = argData[0].partition(' ') # 0-th element contains arg's data
        argName = argName.strip()
        argTitle = argTitle.strip(' ()')
        argSpecification['name'] = argName
        argSpecification['title'] = argTitle

        # 2) get initial value data
        name, ignore, value = argData[1].partition('=') # 1-st element contains arg's initial data
        name = name.strip()
        value = value.strip()
        argSpecification['initial'] = {'name': name, 'value': value} 

        # 3) get final value data
        name, ignore, value = argData[2].partition('=') # 2-nd element contains arg's final data
        name = name.strip()
        value = value.strip()
        argSpecification['final'] = {'name': name, 'value': value}

        # 4) get step data
        name, ignore, value = argData[3].partition('=') # 3-rd element contains arg's step data
        name = name.strip()
        value = value.strip()
        argSpecification['step'] = {'name': name, 'value': value}

        return argSpecification

    def extractQuantitiesSpecification(rawDescription):
        """
        Processes quantities raw description and returns their specification.
        """
        specification = {}
        
        for line in rawDescription:
            
            # extract name
            name, ignore, rest = line.partition('=')
            name = name.strip()

            # extract value and units
            value, ignore, units = rest.partition('(')
            value = value.strip()
            units = units.strip(' )')

            # store data
            specification[name] = {'value': value, 'units': units}

        return specification

    def joinMultilineFormulas(formulas):
        """
        Returns joined multiline formulas.
        """
        joinedFormulas = []        

        currentFormula = formulas[0]

        # each line without the symbol '=' is appended to the previous one
        for i in range(1, len(formulas)):
            if '=' not in formulas[i]:
                currentFormula += ' ' + formulas[i]
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
            name, ignore, formula = line.partition('=')
            specification[name.strip()] = formula.strip() # spaces are removed

        return specification

    def extractDifEqsSpecification(expressions):
        """
        Process raw differential equations' lines and returns their specification.
        """
        specification = {}

        for line in expressions:

            # get left & right parts
            left, ignore, right = line.partition('=')

            # get expression with the name
            numerator, ignore, denominator = left.partition('/')

            # get name: from 'dy' we extract 'y'
            ignore1, ignore, name = numerator.partition('d')           

            # store results: spaces and parenthesis are removed
            specification[name.strip(' ()')] = right.strip()

        return specification

    ivp = {NAME: ivpRaw[NAME], METHOD: ivpRaw[METHOD] }

    # 1. Process argument's description
    ivp[ARGUMENT] = extractArgSpecification(ivpRaw[ARGUMENT])

    # 2. Process initial values' description
    ivp[INIT_VALS] = extractQuantitiesSpecification(ivpRaw[INIT_VALS])

    # 3. Process constants' description
    ivp[CONSTANTS] = extractQuantitiesSpecification(ivpRaw[CONSTANTS])

    # 4. Process parameters' description
    ivp[PARAMS] = extractQuantitiesSpecification(ivpRaw[PARAMS])

    # 5. Process expressions' description
    ivp[EXPRESSIONS] = extractExpressionsSpecification(joinMultilineFormulas(ivpRaw[EXPRESSIONS]))

    # 6. Process differential equations' description
    ivp[DIF_EQUATIONS] = extractDifEqsSpecification(joinMultilineFormulas(ivpRaw[DIF_EQUATIONS])) 


    # 7. Check results TODO: this section should be improved! 

    # 7.1) check differential variables
    if set(ivp[DIF_EQUATIONS].keys()) != set(ivp[INIT_VALS].keys()):
        raise Warning

    return ivp

def generateCppCode(ivp):
    """
    Generates C++-code that solves the problem.
    """

    def putExternBlock(puts, ivp):
        """
        Puts the 'extern "C"' block required for Emscripten compiling.
        """
        puts('\nextern "C" {\n')

        # name of solving function
        puts(f'{SPACE} int {SOLVER_PREFIX}{ivp[NAME]}(')

        # arguments
        puts(f'{DATA_TYPE} {ivp[ARGUMENT]["initial"]["name"]},')
        puts(f' {DATA_TYPE} {ivp[ARGUMENT]["final"]["name"]},')
        puts(f' {DATA_TYPE} {ivp[ARGUMENT]["step"]["name"]},\n')

        # initial values
        line = SUBSPACE

        for name in ivp[INIT_VALS].keys():
            line += f'{DATA_TYPE} _{name}Initial, '

        puts(line + '\n')

        # parameters
        line = SUBSPACE

        for name in ivp[PARAMS].keys():
            line += f'{DATA_TYPE} _{name}Val, '

        puts(line + '\n')

        puts(SUBSPACE + f'int _{ivp[ARGUMENT]["name"]}Count, ')
        puts('int _varsCount,\n')

        # dataframe
        puts(SUBSPACE + f'{DATA_TYPE} * {SOLVER_RESULT_NAME}, ')
        puts(f'int {SOLVER_RESULT_ROW_COUNT_NAME}, ')
        puts(f'int {SOLVER_RESULT_COL_COUNT_NAME}) noexcept;\n')

        puts('}\n')

    def putIVPcode(puts, ivp):
        """
        Puts C++-code specifying IVP.
        """

        # 1. PUT THE RIGHT PART OF ODEs

        puts(f'\n\nnamespace {ivp[NAME]}\n')
        puts('{\n')

        puts(SPACE + '// tollerance\n')
        puts(SPACE + f'{DATA_TYPE} {TOLLERANCE_NAME} = {TOLLERANCE_VALUE};\n\n')

        puts(SPACE + '// dimension of solution\n')
        puts(SPACE + f'const int {DIMENSION_NAME} = {len(ivp[INIT_VALS])};\n\n')

        puts(SPACE + '// parameters\n')

        # put parameters declaration
        for name in ivp[PARAMS]:
            puts(SPACE + f'{ARG_TYPE} {name} = {PARAMETERS_INIT_VALUE};\n')

        # operating name of the argument, for example, t
        argName = ivp['argument']['name']

        puts(f'\n{SPACE}//the right part of the ODEs\n')
        
        # put ODEs right part header, for example, f(...)
        puts(SPACE)
        puts(f'{VEC_TYPE} {ODES_RIGHT_PART_NAME}({ARG_TYPE} {argName}, ')
        puts(f'{VEC_TYPE} & {VEC_ARG_NAME}) noexcept\n')

        # put body of the function
        puts(SPACE + '{\n')

        puts(SUBSPACE + f'{VEC_TYPE} {VEC_OUTPUT_NAME}({DIMENSION_NAME});\n\n')

        puts(SUBSPACE + '// constants\n')

        # put constants
        for name in ivp[CONSTANTS]:
            puts(SUBSPACE + f'const {ARG_TYPE} {name} = {ivp[CONSTANTS][name]["value"]};\n')
      
        index = 0
        puts('\n' + SUBSPACE + '// extract variables\n')

        # transform input vector coordinates to the current names
        for name in ivp[DIF_EQUATIONS]:
            puts(SUBSPACE + f'{ARG_TYPE} {name} = {VEC_ARG_NAME}({index});\n')
            index += 1

        puts('\n' + SUBSPACE + '// expressions\n')

        # put expressions
        for name in ivp[EXPRESSIONS]:
            puts(SUBSPACE + f'{ARG_TYPE} {name} = {ivp[EXPRESSIONS][name]};\n')

        puts('\n' + SUBSPACE + '// output computation\n')
        index = 0

        # computation lines
        for name in ivp[DIF_EQUATIONS]:
            puts(SUBSPACE + f'{VEC_OUTPUT_NAME}({index}) = {ivp[DIF_EQUATIONS][name]};\n')
            index += 1        

        puts('\n' + SUBSPACE + f'return {VEC_OUTPUT_NAME};\n')
        
        puts(SPACE + '}' + f' // {ODES_RIGHT_PART_NAME}\n\n')

        # 1. PUT JACOBIAN

        puts(SPACE + '// Jacobian (it is required, when applying implicit method)\n')

        # header
        puts(SPACE + f'{MATRIX_TYPE} {JACOBIAN_NAME}(')                
        puts(f'{ARG_TYPE} {argName}, ')
        puts(f'{VEC_TYPE} & {VEC_ARG_NAME}, ')
        puts(f'{ARG_TYPE} {EPS_NAME}) noexcept\n')

        # body
        puts(SPACE + '{\n')

        puts(SUBSPACE + f'{MATRIX_TYPE} {MAT_OUTPUT_NAME}({DIMENSION_NAME}, {DIMENSION_NAME});\n')

        puts(SUBSPACE + f'{VEC_TYPE} {VAL_NAME} = {ODES_RIGHT_PART_NAME}')
        puts(f'({argName}, {VEC_ARG_NAME});\n')

        puts(SUBSPACE + f'{VEC_TYPE} {Y_DER_NAME} = {VEC_ARG_NAME};\n\n')

        # put loop-block
        puts(SUBSPACE + f'for (int i = 0; i < {DIMENSION_NAME}; i++)' + ' {\n')
        puts(SUBSUBSPACE + f'{Y_DER_NAME}(i) += {EPS_NAME};\n')
        puts(SUBSUBSPACE + f'{MAT_OUTPUT_NAME}.col(i) = (')
        puts(f'{ODES_RIGHT_PART_NAME}({argName}, {Y_DER_NAME}) - ')
        puts(f'{VAL_NAME}) / {EPS_NAME};\n')
        puts(SUBSUBSPACE + f'{Y_DER_NAME}(i) -= {EPS_NAME};\n')

        puts(SUBSPACE + '}\n\n')

        puts(SUBSPACE + f'return {MAT_OUTPUT_NAME};\n')

        puts(SPACE + '}' + f' // {JACOBIAN_NAME}\n\n')

        # 3. PUT T-DERIVATIVE BLOCK

        puts(SPACE + '// Derivative with respect to t (it is required, when applying implicit method)\n')

        # header
        puts(SPACE + f'{VEC_TYPE} {T_DERIVATIVE_NAME}(')                
        puts(f'{ARG_TYPE} {argName}, ')
        puts(f'{VEC_TYPE} & {VEC_ARG_NAME}, ')
        puts(f'{ARG_TYPE} {EPS_NAME}) noexcept\n')

        # body
        puts(SPACE + '{\n')

        puts(SUBSPACE + f'return ({ODES_RIGHT_PART_NAME}(')
        puts(f'{argName} + {EPS_NAME}, {VEC_ARG_NAME}) - ')
        puts(f'{ODES_RIGHT_PART_NAME}({argName}, {VEC_ARG_NAME}))')
        puts(f' / {EPS_NAME};\n')

        puts(SPACE + '}' + f' // {T_DERIVATIVE_NAME}\n\n')

        puts('}; ' + f'// {ivp[NAME]}\n')

    def putAnnotation(puts, ivp):
        """
        Puts ReachFunctionView annotation.
        """
        
        puts(f'\n//name: {SOLVER_PREFIX + ivp[NAME]}\n')

        # 1. PUT ANNOTATION CONCERNING ARGUMENT

        # initial value
        puts(f'//input: {ANNOT_ARG_TYPE} {ivp[ARGUMENT]["initial"]["name"]}')
        puts(f' = {ivp[ARGUMENT]["initial"]["value"]}')
        puts(' {caption: ' + ivp[ARGUMENT]["initial"]["name"] + ';')
        puts(f' category: {ivp[ARGUMENT]["title"]}' + '}\n')

        # final value
        puts(f'//input: {ANNOT_ARG_TYPE} {ivp[ARGUMENT]["final"]["name"]}')
        puts(f' = {ivp[ARGUMENT]["final"]["value"]}')
        puts(' {caption: ' + ivp[ARGUMENT]["final"]["name"] + ';')
        puts(f' category: {ivp[ARGUMENT]["title"]}' + '}\n')

         # step value
        puts(f'//input: {ANNOT_ARG_TYPE} {ivp[ARGUMENT]["step"]["name"]}')
        puts(f' = {ivp[ARGUMENT]["step"]["value"]}')
        puts(' {caption: ' + ivp[ARGUMENT]["step"]["name"] + ';')
        puts(f' category: {ivp[ARGUMENT]["title"]}' + '}\n')

        # 2. INITIAL VALUES

        # dict with initial values specification
        data = ivp[INIT_VALS]

        # put annotation
        for name in data.keys():
            puts(f'//input: {ANNOT_ARG_TYPE} _{name}Initial = {data[name]["value"]}')
            puts(' {' + f'units: {data[name]["units"]};')
            puts(f' caption: {name}; category: {INIT_VALS}' + '}\n')
        
        # 3. PARAMETERS

        # dict with parameters specification
        data = ivp[PARAMS]

        # put annotation
        for name in data.keys():
            puts(f'//input: {ANNOT_ARG_TYPE} _{name}Val = {data[name]["value"]}')
            puts(' {' + f'units: {data[name]["units"]};')
            puts(f' caption: {name}; category: {PARAMS}' + '}\n')

        # 4. SPECIAL INPUTS
        puts(f'//input: int _{ivp[ARGUMENT]["name"]}Count\n')
        puts(f'//input: int _varsCount\n')

        # 5. OUTPUT COLUMNS
        puts(f'//output: column_list {SOLVER_RESULT_NAME} ')
        puts(f'[new(_{ivp[ARGUMENT]["name"]}Count, _varsCount)')

        # create new columns names
        line = "; '" + ivp[ARGUMENT]["name"] + "'"       

        for name in ivp[DIF_EQUATIONS]:
            line += "; '" + name + f'({ivp[ARGUMENT]["name"]})' + "'"

        puts(line + ']\n')     

        # 6. OUTPUT DATAFRAME
        puts(f'//output: dataframe {DF_SOLUTION_NAME} [{SOLVER_RESULT_NAME}] {DF_OUTPUT_ANNOT}\n')

        # 7. EDITOR
        puts(EDITOR_LINE + '\n')
    
    def putSolverHeader(puts, ivp):
        """
        Puts IVP solving function header.
        """
        puts('EMSCRIPTEN_KEEPALIVE\n')

        # name of solving function
        puts(f'int {SOLVER_PREFIX}{ivp[NAME]}(')

        # arguments
        puts(f'{DATA_TYPE} {ivp[ARGUMENT]["initial"]["name"]},')
        puts(f' {DATA_TYPE} {ivp[ARGUMENT]["final"]["name"]},')
        puts(f' {DATA_TYPE} {ivp[ARGUMENT]["step"]["name"]},\n')

        # initial values
        line = SPACE

        for name in ivp[INIT_VALS].keys():
            line += f'{DATA_TYPE} _{name}Initial, '

        puts(line + '\n')

        # parameters
        line = SPACE

        for name in ivp[PARAMS].keys():
            line += f'{DATA_TYPE} _{name}Val, '

        puts(line + '\n')

        puts(SPACE + f'int _{ivp[ARGUMENT]["name"]}Count, ')
        puts('int _varsCount,\n')

        # dataframe
        puts(SPACE + f'{DATA_TYPE} * {SOLVER_RESULT_NAME}, ')
        puts(f'int {SOLVER_RESULT_ROW_COUNT_NAME}, ')
        puts(f'int {SOLVER_RESULT_COL_COUNT_NAME}) noexcept\n')
    
    def putSolverBody(puts, ivp):
        """
        Puts IVP solving function body.
        """
        puts('{\n')

        puts(SPACE + f'using namespace {ivp[NAME]};\n\n')

        puts(SPACE + f'{DATA_TYPE} {INIT_VALS_ARR_NAME}[DIM];\n\n')
        
        puts(SPACE + '// initial values\n')

        index = 0

        # put initial values 
        for name in ivp[INIT_VALS]:
            puts(SPACE + f'{INIT_VALS_ARR_NAME}[{index}] = _{name}Initial;\n')
            index += 1

        puts('\n' + SPACE + '// parameters\n')

        # put parameters
        for name in ivp[PARAMS]:
            puts(SPACE + name + f' = _{name}Val;\n')

        # return line
        line = '\n' + SPACE + f'return solveODE({ODES_RIGHT_PART_NAME}, ' 

        if ivp[METHOD] == 'implicit':
            line += f'{T_DERIVATIVE_NAME}, {JACOBIAN_NAME}, '

        line += f'{ivp[ARGUMENT]["initial"]["name"]}, '
        line += f'{ivp[ARGUMENT]["final"]["name"]}, '
        line += f'{ivp[ARGUMENT]["step"]["name"]}, '
        line += f'{INIT_VALS_ARR_NAME}, {TOLLERANCE_NAME},\n' + SUBSUBSPACE
        line += f'{SOLVER_RESULT_NAME}, _{ivp[ARGUMENT]["name"]}Count, _varsCount);\n'
        
        puts(line)
        
        puts('} //' + SOLVER_PREFIX + ivp[NAME])


    
    settings = {}

    # load settings
    with open(EXPORT_SETTINGS_FILE, 'r') as file:
         settings = json.load(file)
    
    # generate C++-code
    with open(settings['folder'] + '/' + ivp[NAME] + '.cpp', 'w') as file:
        puts = file.write

        puts(f'// {ivp[NAME]}.cpp\n')

        # STL includes TODO: consider other libs
        puts('\n#include <cmath>\n')
        puts('using namespace std;\n')

        # Emscripten tools
        puts(f'\n#include <{settings["emscripten"]}>\n')

        # the 'extern "C"' block - it's required for Emscripten
        putExternBlock(puts, ivp)

        # Eigen tools
        puts(f'''\n#include "{settings['eigen']}"\n''')
        puts('using namespace Eigen;\n')

        # ODEs solver tools
        puts(f'\n#include "{ODES_SOLVER_LIB}"\n')
        puts(f'using namespace {ODES_SOLVER_NAMESPACE};')

        # IVP code
        putIVPcode(puts, ivp)

        # RichFunctionView annotation
        putAnnotation(puts, ivp)

        # solving function's header
        putSolverHeader(puts, ivp)

        # solving function's body
        putSolverBody(puts, ivp)


def updateExportSettings(ivp):
    """
    Update C++-to-wasm export settings.
    """
    settings = {}

    with open(EXPORT_SETTINGS_FILE, 'r') as file:
         settings = json.load(file)
    
    settings['name'] = ivp[NAME]
    settings['source'] = [ivp[NAME] + '.cpp']

    with open(EXPORT_SETTINGS_FILE, 'w') as file:
            json.dump(settings, file)


def main(ivpFileName="IVP.txt"):
    """
    The main script: 
        1) get initial value problem from file;
        2) get specification of IVP from its raw description;
        3) generate C++ code;
        4) update settings for C++-to-wasm export;
        
    """
    try:
        # 1) get initial value problem from file
        ivpRaw = getRawIVPfromFile(ivpFileName)       

        with open('rawIVP.json', 'w') as file:
            json.dump(ivpRaw, file)
        
        # 2) get specification of IVP from its raw description
        ivp = getSpecificationFromRawIVP(ivpRaw)

        with open('IVP.json', 'w') as file:
            json.dump(ivp, file)

        # 3) generate C++ code
        generateCppCode(ivp)

        # 4) update settings for C++-to-wasm export
        updateExportSettings(ivp)



    except OSError:
        print("OS ERROR: check paths and file names!")
    
    except KeyError:
        print("PARSING ERROR: check specificaition of the initial value problem!")

    except IndexError:
        print("PARSING ERROR: check C/C++-code!")

    except Warning:
        print("ADD")


if __name__ == '__main__':
    main()
