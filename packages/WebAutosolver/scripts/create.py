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
            line += f'{DATA_TYPE} {name}, '

        puts(line + '\n')

        # parameters
        line = SUBSPACE

        for name in ivp[PARAMS].keys():
            line += f'{DATA_TYPE} {name}, '

        puts(line + '\n')


        puts(SUBSPACE + f'{DATA_TYPE} * result, int resultRowCount, int resultColCount) noexcept;\n')

        puts('}\n')

    def putIVPcode(puts, ivp):
        """
        Puts C++-code specifying IVP.
        """
        pass

    def putAnnotation(puts, ivp):
        """
        Puts ReachFunctionView annotation.
        """
        pass
    
    def putSolverHeader(puts, ivp):
        """
        Puts IVP solving function header.
        """
        pass
    
    def putSolverBody(puts, ivp):
        """
        Puts IVP solving function body.
        """
        pass
    
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
