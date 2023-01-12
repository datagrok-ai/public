""" export.py  January 12, 2023.   
    This script exports C/C++-functions to Datagrok package. 
"""

import os
import json

# auxiliry map 
sizesMap = {'rowCount': 'numOfRows', 'columnCount': 'numOfColumns', 'data': 'data'}

def getSettings(nameOfFile):
    """
    Return dictionary with settings of C/C++-functions export from json-file with settings.
    """
    with open(nameOfFile, 'r') as file:
         return json.load(file)

class Function:
    """ Exported C/C++-function specification class.
    """
    def __init__(self, cFuncAnnotationLines, typeList):        
        self.annotation = []
        self.arguments = {}
        self.prototype = ''
        self.callArgs = '['
        self.typeList = typeList
        self.cFuncTypeIndex = 0
        self.output = {}
        self.paramsWithNewSpecifier = {}

        # process each line of annotation
        for line in cFuncAnnotationLines:
            self.processAnnotationLine(line)

        # complete output of specification
        self.finalizeOutput()      

        # final preparation of prototype
        self.prototype = self.prototype.rstrip(', ')
        self.prototype += ')'

        # final preparation of call args string
        self.callArgs = self.callArgs.rstrip(', ')
        self.callArgs += ']'

        # create specification dict
        self.specification = {'arguments': self.arguments, 
        'output': self.output, 
        'annotation': self.annotation,
        'prototype': self.prototype,
        'callArgs': self.callArgs}

    def processAnnotationLine(self, line):
        """ Process one line of annotation before exported C/C++-function.
        """
        if line.startswith('//input:'):
            self.processInputLine(line) 
        elif line.startswith('//output:'):
            self.processOutputLine(line)
        elif line.startswith('//name:'):
            self.processNameLine(line)        
        else:
            self.annotation.append(line)

    def processNameLine(self, line):
        """ Process line containing '//name:'.
        """
        ignore1, ignore2, name = line.partition(' ')        
        self.prototype += (name + '(')
        self.annotation.append(line)

    def processOutputLine(self, line):
        """ Process line containing '//output:'.
        """
        if 'new' in line:
            self.processOutputLineWithNew(line)
        else:
            self.processOutputLineWithoutNew(line)        

    def processInputLine(self, line):
        """ Process '//input'-line with arguments that are passed.
        """
        self.annotation.append(line)

        # get type and name part of line
        ignore1, ignore2, rest = line.partition(': ')        

        # get type and name of the argument
        argType, ignore, argName = rest.partition(' ')

        self.prototype += argName + ', '

        if argType == 'dataframe':
            return

        self.callArgs += argName + ', '

        currentType = self.typeList[self.cFuncTypeIndex]        

        if argType == 'column':
            self.arguments[argName] = {'type': currentType + 'Column'}
            self.cFuncTypeIndex += 2

        elif argType == 'column_list':
            self.arguments[argName] = {'type': currentType + 'Columns'}
            self.cFuncTypeIndex += 3

        else:
            self.arguments[argName] = {'type': 'num'}
            self.cFuncTypeIndex += 1

    def processOutputLineWithNew(self, line):
        """ Process '//output'-line with arguments that are newly created.
        """
         # get type and name part of line
        ignore1, ignore2, rest = line.partition(': ')        

        # get type and name of the argument
        argType, ignore, rest = rest.partition(' ')

        if argType == 'dataframe':
            return

        # extract name & specification
        argName, ignore, specification = rest.partition(' ')

        # store argument info in dictionary
        self.paramsWithNewSpecifier[argName] = argType

        # extract the main data
        specification = specification.lstrip('[new(').rstrip(']')        
    
        currentType = self.typeList[self.cFuncTypeIndex]
        
        if argType == 'column':
            specification = specification.strip()

            ref = None
            value = None

            if '.' not in specification:
                ref = specification
                value = 'data'
            else:
                ref, ignore, value = specification.partition('.')

            self.arguments[argName] = {'type': 'new' + currentType.capitalize() + 'Column'}
            self.arguments[argName]['numOfRows'] = {'ref': ref, 'value': sizesMap[value]} 
            self.cFuncTypeIndex += 2

        elif argType == 'column_list':
            
            first, ignore, rest = specification.partition(',')
            first = first.strip()
            rest = rest.strip()          

            self.arguments[argName] = {'type': 'new' + currentType.capitalize() + 'Columns'}

            ref = None
            value = None

            if '.' not in first:
                ref = first
                value = 'data'
            else:
                ref, ignore, value = first.partition('.')

            self.arguments[argName]['numOfRows'] = {'ref': ref, 'value': sizesMap[value]}

            colCountStr, ignore, namesStr = rest.partition(')')
            colCountStr = colCountStr.strip()

            if '.' not in colCountStr:
                ref = colCountStr
                value = 'data'
            else:
                ref, ignore, value = colCountStr.partition('.')                        

            self.arguments[argName]['numOfColumns'] = {'ref': ref, 'value': sizesMap[value]}

            namesStr = namesStr.strip(';').strip().replace(';', ',')           

            self.arguments[argName]['names'] = namesStr

            self.cFuncTypeIndex += 3

        else:
            self.arguments[argName] = {'type': 'num'}
            self.cFuncTypeIndex += 1

    def processOutputLineWithoutNew(self, line):
        """ Process '//output'-line.
        """
        ignore1, ignore2, specification = line.partition(':')       

        outputType, ignore, outputSource = specification.strip().partition(' ')

        outputType = outputType.strip()

        ignore1, ignore2, outputSource = outputSource.partition(' ')

        outputSource = outputSource.strip().strip('[]')
        
        if outputType == 'objects':
            self.output['type'] = outputType
            args = outputSource.split(',')
            sourceString = ''

            for arg in args:
                sourceString += f"'{arg.strip()}', "           

            self.output['source'] = f"[{sourceString.rstrip(', ')}]"

            return

        elif outputType == 'dataframe':
            self.output['type'] = 'tableFromColumns'
            self.output['source'] = f'{outputSource}'

        else:
            self.output['type'] = f'{outputType}'
            self.output['source'] = f'{outputSource}'

        stringToAnnotation, ignore1, ignore2 = line.partition('[')
        self.annotation.append(stringToAnnotation)

    def finalizeOutput(self):
        """ Complete output part of the function specification.
        """
        if len(self.output) == 0:
            if len(self.paramsWithNewSpecifier) == 0:
                raise Warning
            elif len(self.paramsWithNewSpecifier) == 1:
                argName = list(self.paramsWithNewSpecifier.keys())[0]
                self.output['type'] = self.paramsWithNewSpecifier[argName]
                self.output['source'] = argName
                self.annotation.append(f'//output: {self.output["type"]} {argName}')
            else:
                self.output['type'] = 'objects'
                sourceString = '['

                for name in self.paramsWithNewSpecifier.keys():
                    sourceString += f"'{name}'" + ', '

                sourceString = sourceString.rstrip(', ')
                sourceString += ']'
                self.output['source'] = sourceString

    def getSpecification(self):
        return self.specification

def getFunctionsFromFile(nameOfFile):
    """
    Parse C/C++-file and return specification of exported functions.
    """    
    with open(nameOfFile, 'r') as file:
        listOfLines = file.readlines()
        dictionaryOfFunctions = {}     

        for lineIndex in range(len(listOfLines)):
            if listOfLines[lineIndex] == "EMSCRIPTEN_KEEPALIVE\n":  # the next line is function's prototype

                # 1. Analyze function description
                j = lineIndex + 1
                funcLine = ''

                # 1.1) get one line prototype of exported function
                while True:
                    funcLine += listOfLines[j].strip().rstrip('\n')                    
                    if ')' in listOfLines[j]:
                        break
                    j += 1                

                # 1.2) get name of function & a string of its arguments
                ignore1, ignore2, nameAndArgs = funcLine.partition(' ')                
                functionName, ignore, argsLine = nameAndArgs.partition('(')
                
                # 1.3) get types of arguments
                cFuncArgs = argsLine.split(',')
                typeList = [] # list of argument's type 

                for s in cFuncArgs:
                    argType, ignore1, ignore2 = s.strip().partition(' ')
                    typeList.append(argType)

                # 1.4) get annotation lines
                j = lineIndex - 1
                annotationLines = []                

                while listOfLines[j].startswith('//'):
                    annotationLines.append(listOfLines[j].strip())
                    j -= 1

                annotationLines.reverse()                
                function = Function(annotationLines, typeList)
                dictionaryOfFunctions[functionName] = function.getSpecification()                       

        return dictionaryOfFunctions                      
      
def getFunctionsFromListOfFiles(settings):
    """ Loads a list of C-file names and gets functions from each of these files
    """
    functions = {} # dictionary: key - name of file; value - functions data  
    
    for name in settings["source"]:
        functions[name] = getFunctionsFromFile(settings["folder"] + '/' + name)

    return functions

def saveFunctionsDataToFile(functionsData, nameOfFile):
    """
    Save C-functions descriptors to json-file.
    """
    with open(nameOfFile, 'w') as file:
            json.dump(functionsData, file)

def getCommand(settings, functionsData):
    """
    Return Emscripten command based on settings and functions data.
    """
    command = 'em++' # Emscripten command, starts with "emcc"

    # set optimization mode
    command += ' ' + settings['optimizationMode']

    # add names of C-files
    for nameOfFile in settings['source']:
        command += f' {settings["folder"]}/{nameOfFile}'        
    
    # add name of JS-library that is created by Emscripten
    command += f' -o {settings["folder"]}/{settings["name"]}.js'     

    # specify memory restrictions
    command += ' -s TOTAL_MEMORY=' + settings['totalMemory']

    # set that wasm-file must be created
    command += ' -s WASM=1'

    # allow memory growth
    command += ' -s ALLOW_MEMORY_GROWTH=1'

    # create module (? should be checked!)
    command += ' -s MODULARIZE=1'

    # set export name
    command += ' -s EXPORT_NAME="export' + settings['name'] + '"'

    # set a list of exported functions
    command += ' -s EXPORTED_FUNCTIONS=['

    # add names of functions
    for nameOfFile in functionsData.keys():
        for nameOfFunction in functionsData[nameOfFile].keys():
            command += '"_' + nameOfFunction + '",'
    
    # add malloc and free - they are used for memory operating
    command += '"_malloc","_free"]' 

    # add exported runtime methods
    command += ' -s EXPORTED_RUNTIME_METHODS=["cwrap","ccall"]'  # also, "ccall" can be added 

    return command
   
def saveCommand(command, nameOfFile):
    """
    Save Emscripten command to txt-file.
    """
    with open(nameOfFile, "w") as file:
        file.write(command)

def updatePackageJsonFile(settings):
    """ Add JS-file with exported C/C++-functions to "sources" of package.json.
    """
    # data from package.json
    packageData = None 

    # get package data from package.json 
    with open(settings['packageJsonFile'], 'r') as file:
        packageData = json.load(file)

    # create full name of JS-lib file that is added to package.json
    fullNameOfLibFile = f"{settings['folder'].lstrip('../')}/{settings['name']}.js"

    # add dependence to package data
    if "sources" in packageData.keys():

        if fullNameOfLibFile not in packageData["sources"]: # add JS-file if "source" does not contain it yet
            packageData["sources"].append(fullNameOfLibFile)

    else:
        packageData["sources"] = [ fullNameOfLibFile ]

    # update package.json
    with open(settings["packageJsonFile"], 'w') as file:
        json.dump(packageData, file)

def completeJsWasmfile(settings, functionsData):
    """
    Add exported C/C++-function specifications and the corresponding init-function
    to JS-file created by Emscripten.    
    """
    with open(f'{settings["folder"]}/{settings["name"]}.js', 'a') as file:
        put = file.write
        put('\n\n')

        exportedFuncsNames = []

        for cppFile in functionsData.keys():
            for funcName in functionsData[cppFile].keys():
                exportedFuncsNames.append(funcName)
                put(f'var {funcName} = ' + '{\n')
                arguments = functionsData[cppFile][funcName]['arguments']
                put('  arguments: {\n')
                argCommasCount = len(arguments.keys()) - 1

                for arg in arguments.keys():
                    put(f'    {arg}: ' + '{\n')
                    attrCommasCount = len(arguments[arg].keys()) - 1

                    for attribute in arguments[arg].keys():                        

                        if attribute == 'type':
                            put(f"      {attribute}: '{arguments[arg][attribute]}'")

                        elif attribute == 'names':
                            put(f'      {attribute}: [{arguments[arg][attribute]}]')                           

                        else:
                            put(f'      {attribute}: ' + '{\n')
                            ref = arguments[arg][attribute]['ref']
                            put('        ref: ' + f"'{ref}',\n")
                            value = arguments[arg][attribute]['value']
                            put('        value: ' + f"'{value}'\n")
                            put('      }')

                        if attrCommasCount > 0:
                            put(',')
                            attrCommasCount -= 1
                                
                        put('\n')

                    put('    }')

                    if argCommasCount > 0:
                        put(',')
                        argCommasCount -= 1

                    put('\n')

                put('  },\n')
                put('  output: {\n')
                outputType = functionsData[cppFile][funcName]['output']['type']
                put(f"    type: '{outputType}',\n")
                outputSource = functionsData[cppFile][funcName]['output']['source']
                put('    source: ')

                if outputType != 'objects':
                    put(f"'{outputSource}'")
                else:
                    put(f'{outputSource}')
                
                put('\n  }\n')
                put('}; ' + f'// {funcName}\n\n')

        # add init-function
        put(f'var {settings["name"]} = undefined;\n\n')
        put(f'async function init{settings["name"]}() ' + '{\n')
        put(f'  if ({settings["name"]} === undefined) ' + '{\n')
        put(f'    console.log("Wasm not Loaded, Loading");\n')
        put(f'    {settings["name"]} = await export{settings["name"]}();\n')

        for name in exportedFuncsNames:
            put(f'    {settings["name"]}.{name} = {name};\n')
        
        put('  } else {\n')
        put('    console.log("Wasm Loaded, Passing");\n')
        put('  }\n')
        put('}') 
        
def addExportedFunctionsToPackageFile(settings, functionsData):
    """ Generate code in the package-file.
    """
    with open(settings['packageFile'], 'a') as file:
        put = file.write

        # put 'import...' line
        put('\nimport { callWasm } from ')
        callWasm = settings['runtimeSystemFile'].rstrip('.js')
        put(f"'{callWasm}';\n\n")

        # put init-function
        put('//tags: init\n')
        put('export async function init() {\n')
        put(f'  await init{settings["name"]}();\n')
        put('}\n\n')

        # put call of each exported C/C++-function
        for cppFileName in functionsData.keys():
            for funcName in functionsData[cppFileName].keys():
                annotation = functionsData[cppFileName][funcName]['annotation']
                prototype = functionsData[cppFileName][funcName]['prototype']
                callArgs = functionsData[cppFileName][funcName]['callArgs']

                # put annotation lines
                for line in annotation:
                    put(line + '\n')

                # put package function code
                put(f'export function {prototype} ' + '{\n')
                put(f"  return callWasm({settings['name']}, '{funcName}', {callArgs});\n")
                put('}\n\n')


def main(nameOfSettingsFile="module.json"):
    """
    The main script: 
        1) load export settings from json-file;
        2) parse C/C++-files and get C/C++-functions to be exported;        
        3) create Emscripten command;        
        4) execute Emscripten command;
        5) complete JS-file created by Emscripten;
        6) append Datagrok package file with exported functions;
        7) update the file package.json: add JS-file name created by Emscripten to 'sources'.
    """
    try:
        # 1) load settings
        settings = getSettings(nameOfSettingsFile)        

        # 2) parse C-files and get exported functions data
        functionsData = getFunctionsFromListOfFiles(settings)
    
        # write functions descriptors to json-file
        saveFunctionsDataToFile(functionsData, 'func.json')

        # 3) get command for Emscripten tool
        command = getCommand(settings, functionsData)  

        # save command to file
        #saveCommand(command, 'command.txt')  

        # 4) execute command by Emscripten
        os.system(command)

        # 5) complete JS-file created by Emscripten
        completeJsWasmfile(settings, functionsData) 

        # 6) append Datagrok package file with exported functions
        addExportedFunctionsToPackageFile(settings, functionsData) 

        # 7) update the file package.json
        updatePackageJsonFile(settings)

    except OSError:
        print("OS ERROR: check paths and file names!")
    
    except KeyError:
        print("PARSING ERROR: check C/C++-code!")

    except IndexError:
        print("PARSING ERROR: check C/C++-code!")

    except Warning:
        print("Check annotation of exported C/C++-function!")


if __name__ == '__main__':
    main()