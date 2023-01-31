""" create.py January 21, 2023.
    This script parses a file with initial value problem (IVP) 
    and creates webautosolver Datagrok function.
"""

import os
import json

# Names of special categories specifying IVP
NAME = 'name'
METHOD = 'method'
DIF_EQUATIONS = 'differential equations'
EXPRESSIONS = 'expressions'
ARGUMENT = 'argument'
INIT_VALS = 'initial values'
CONSTANTS = 'constants'
PARAMS = 'parameters'


def getIVPfromFile(ivpFileName):
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
            if '#' in line:                
                category, ignore, rest = line.partition(':') # extract category data
                category = category.lstrip('#').strip() # extract name of category

                # process specific categories
                if category in {NAME, METHOD}:                
                    ivpRaw[category] = rest.strip()
                elif category == ARGUMENT:
                    ivpRaw[ARGUMENT].append(rest.strip())                        

            else:
                line = line.strip()  
                print(category)              
                if(line != ''):
                    ivpRaw[category].append(line)        

    return ivpRaw


def main(ivpFileName="IVP.txt"):
    """
    The main script: 
        1) get initial value problem from file;
        
    """
    try:
        # 1) get initial value problem from file
        ivpRaw = getIVPfromFile(ivpFileName)       

        with open('rawIVP.json', 'w') as file:
            json.dump(ivpRaw, file)

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
