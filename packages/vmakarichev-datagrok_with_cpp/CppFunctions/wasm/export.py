""" export.py  Release 3.0.0, Ocotber 21, 2022

    SMART EXPORT SCRIPT
    
    This script exports C-functions to Datagrok package

    The following steps are performed: 
        1) load export settings from json-file;
        2) parse C-files and get C-functions to be exported;
        3) save C-functions descriptors (name, type of return value, arguments descriptors: type, name);
        4) create Emscripten command;
        5) save Emscripten command to text file;
        6) execute Emscripten command;
        7) append Datagrok package file with exported functions;
        8) add dependecnies to the file package.json;
        9) add module init-function to JS-file with exported C-functions. 
    """

import os
import json

# map between C- and JavaScript-types
typesMap = {"int *": "Int32Array", "float *": "Float32Array"}

# sets of types that represents data structures
structureTypes={"column", "array", "columns", "arrays"}
simpleStructureTypes = {"column", "array"}
complexStructureTypes = {"columns", "arrays"}

# name of C-functions wrapper
C_FUNC_WRAPPER_NAME = "cFuncWrapper"

def getFunctionsFromFile(nameOfFile):
    """Process the C-file and returns dictionary of functions, type of return, types of arguments.
       nameOfFile - name of C-file.
    """
    
    with open(nameOfFile, 'r') as file:

        listOfLines = file.readlines()

        numOfLines = len(listOfLines)

        i = 0 # index of the current line

        dictionaryOfFunctions = {}

        # analyse all lines of file: look for functions marked with EMSCRIPTEN_KEEPALIVE
        while i < numOfLines:

            if listOfLines[i] == "EMSCRIPTEN_KEEPALIVE\n":  # the next line is function's prototype

                # 1. Analyze function description
                
                functionPrototype = listOfLines[i + 1]
                #print("Prototype: ", functionPrototype, end='')

                # detection of function descriptor and a list of arguments
                functionDescriptor, ignore, rest = functionPrototype.partition("(")
                arguments, ignore, trash = rest.partition(")")

                typeOfReturnValue, ignore, nameOfFunction = functionDescriptor.rpartition(" ")
                #print(f"returnType: {typeOfReturnValue}, name of function: {nameOfFunction}")

                listOfArgumentDescriptors = arguments.split(',') 

                """ Dictionary of the current function args that are numbers.
                    key is a name of argument, 
                    value is a dictionary: {"type": typeName, 
                                            "map": nameOfArgumentFromJS}
                """
                variablesData = {} 

                """ Dictionary of the current function args that are arrays.
                    key is a name of argument, 
                    value is a dictionary: {"type": typeName, 
                                            "map": nameOfArgumentFromJS, 
                                            "length": lengthOfArray}
                """
                arraysData = {}

                cFuncArguments = [] # list of arguments of the function.

                # process each argument descriptor: get name, type, ...
                for argument in listOfArgumentDescriptors:

                    typeDescriptor, ignore, nameOfArgument = argument.rpartition(" ")
                    typeDescriptor = typeDescriptor.strip(" ")

                    # store name of the current argument
                    cFuncArguments.append(nameOfArgument)

                    # check category of the argument: array or number 
                    if '*' in typeDescriptor:
                        # argument is an array
                        arraysData[nameOfArgument] = {"type": typeDescriptor, 
                        "map": None, "length": None}
                    else:
                        # argument is a number, not an array
                        variablesData[nameOfArgument] = {"type": typeDescriptor,
                        "map": None}

                # 2. Analyze preambula

                j = i - 1                

                inputs = {} # dictionary of inputs data

                outputs = {} # dictionary of outputs data

                packageFunctionPreambula = [] # list of lines that is before each Datagrok function

                # dictionary that contains package function data
                packageFunctionData = {"name": None, "export": True,
                "input": inputs, "output": outputs, "preambula": None }
 
                j = i - 1 # number of the line, previous to EMSCRIPTEN_KEEPALIVE

                # consider each line that starts with '//'
                while listOfLines[j].startswith("//"): 

                    line = listOfLines[j].rstrip('\n') # current string without '\n'

                    packageFunctionPreambula.append(line)

                    j -= 1

                # reverse operation is used, since its was obtained in a reverse manner
                packageFunctionPreambula = packageFunctionPreambula[::-1]

                packageFuncArguments = []

                k = 0 # index of the current C argument

                # process package function's preambula
                for line in packageFunctionPreambula:

                    if line.startswith("//name:"):
                        ignore1, ignore2, rest = line.partition(": ")

                        # check function status: export or not
                        if "non-export" in rest: 
                            packageFunctionData["export"] = False

                            nameOfPackageFunc, ignore1, ignore2 = rest.partition(' ')
                            packageFunctionData["name"] = nameOfPackageFunc

                        else:
                            packageFunctionData["name"] = rest

                    elif line.startswith("//input:"):
                        ignore1, ignore2, argsData = line.partition(": ")
                        argType, ignore1, argName = argsData.partition(" ")
                        inputs[argName] = {"type": argType, "update": None}
                        packageFuncArguments.append(argName)

                        if argType == "dataframe": # i.e. this is dataframe
                            pass
                        elif argType in structureTypes: # i.e. this is a column or array
                            arraysData[ cFuncArguments[k] ]["map"] = argName
                            k +=1
                            variablesData[ cFuncArguments[k] ]["map"] = argName + ".length"
                            k += 1
                        else: # i.e. this is a number
                            variablesData[ cFuncArguments[k] ]["map"] = argName
                            k += 1                    

                    elif line.startswith("//output:"):
                        ignore1, ignore2, argsData = line.partition(": ")
                        argType, ignore1, rest = argsData.partition(" ")

                        if '[' not in rest: # rest doesn't contain an expression like [...], i.e. simple name             
                            outputs[rest] = {"type": argType}

                        else:
                            # get expression in [...]
                            line, ignore1, ignore2 = line.partition(" [")
                            argName, ignore, rest = rest.partition(' [')
                            rest,ignore1,ignore2 = rest.partition(']')                            

                            cArgNames = []  # list of C-arguments that defines output

                            prevLength = None  # length of the previous array <- it is used to make preambula shorter

                            for item in rest.split(','):
                                item = item.strip(" ")

                                if 'new' in item:
                                    name, ignore, initCommand = item.partition(" ")
                                    cArgNames.append(name)
                                    ignore1, ignore2, currentLength = initCommand.partition('(')
                                    currentLength, ignore1, ignore2 = currentLength.partition(')')

                                    if currentLength == '...':  # i.e. it is the same as th previous one
                                        currentLength = prevLength
                                    else:
                                        prevLength = currentLength
                                    
                                    variablesData[name + "Length"]["map"] = currentLength

                                    arraysData[name]["length"] = currentLength

                                else:                                    
                                    cArgNames.append(item)
                            
                            # check number of C-arguments
                            if len(cArgNames) > 1:
                                outputs[argName] = {"type": argType, "update": cArgNames}
                            else:
                                outputs[argName] = {"type": argType, "update": cArgNames[0]}
 
                
                # check the last line of preambule: it may containe smth like [...] - it should be removed

                m = len(packageFunctionPreambula) - 1  # index of last line of preambula

                # remove '[...]' from preambula last line, if exists
                if '[' in packageFunctionPreambula[m]:
                    packageFunctionPreambula[m], ignore1, ignore2 = packageFunctionPreambula[m].partition(" [")


                # add package function preambula: //name: ..., etc.
                packageFunctionData["preambula"] = packageFunctionPreambula

                dictionaryOfFunctions[nameOfFunction] = {"type": typeOfReturnValue, 
                "arguments": {"variables": variablesData, "arrays": arraysData},
                "orderOfArguments": cFuncArguments,
                "packageFunction": packageFunctionData}                

            i += 1
        
        #print(dictionaryOfFunctions)

        return dictionaryOfFunctions

def getFunctionsFromListOfFiles(listOfFileNames):
    """ Loads a list of C-file names and gets functions from each of these files
    """
    functions = {} # dictionary: key - name of file; value - functions data

    """with open(nameOfFile, 'r') as file:        
        for name in file.readlines():
            name = name.rstrip('\n') # new line symbol at the end is deleted

            functions[name] = getFunctionsFromFile(name)"""
    
    for name in listOfFileNames:
        functions[name] = getFunctionsFromFile(name)

    return functions

def appendPackageFileWithCfunctions(nameOfFile, moduleName, functions):
    """ Append Datagrok package file with C-functions.
        nameOfFile - name of Datagrok package file to be appended with C-functions;
        moduleName - name of JS-module that will be used in the package file;
        functions - functions data.

        For each exported C-function the following steps are performed:
          1) writing preambula of package function: //name: ...
          2) creating a string of arguments, for example, "arg1, arg2, ..."
          3) writing the first line of function, for example, "function minOfColumn(col) {"
          4) generating and writing a body of the function:
             4.1) preparing variables (non-arrays);
             4.2) preparing arrays that are passed to C-function: memory allocation, etc.;
             4.3) executing C-function;
             4.4) updating variables & arrays;
             4.5) creating output (if it exists)
             4.6) cleaning allocated memory;
             4.7) returning output (if it exists)
    """

    # open package file and append it with C-functions
    with open(nameOfFile, 'a') as file:

        put = file.write # a short form of the function that writes data to file

        put('\n' * 3 + '// EXPORTED C-FUNCTIONS')

        # consider each exported file
        for cFileName in functions.keys():

            # write comment-line with a name of C-file
            put('\n' * 2 + "// Functions from " + cFileName + '\n\n')

            functionsFromCurrentFile = functions[cFileName]

            # add each exported function
            for nameOfFunction in functionsFromCurrentFile.keys():

                # 0. Extract the current function data
                functionData = functionsFromCurrentFile[nameOfFunction]
                arguments = functionData["arguments"] 
                typeOfOutput = functionData["type"]
                packageFunctionData = functionData["packageFunction"]
                isExportFunc = packageFunctionData["export"]

                # 1. Write preambula of the current package function
                if isExportFunc:
                    for line in packageFunctionData["preambula"]:
                        put(line + '\n')
                else:
                    put('/*\n')  # in this case, annotation is written in /*...*/
                    for line in packageFunctionData["preambula"]:
                        put('   ' + line.strip("//") + '\n')
                    put("*/\n")

                # 2. Prepare string of arguments of the package function
                stringOfArguments = ""

                for argument in packageFunctionData["input"].keys():
                    stringOfArguments += argument + ", "

                # remove redundant ',' at the end
                stringOfArguments = stringOfArguments.rstrip(", ")


                # 3. Write the first line of package function
                if isExportFunc:
                    put("export function " + packageFunctionData["name"] + "(" + stringOfArguments + ") {\n")
                else:
                    put("function " + packageFunctionData["name"] + "(" + stringOfArguments + ") {\n")


                # 4. GENERATION OF BODY OF PACKAGE FUNCTION

                """ Dictionary that containes correspondance of C- 
                    and JS/TS- names of variables:
                        key - name of C-variable,
                        value - the corresponding JS/TS- name.
                    For the case of array, value is a name of pointer.
                    Otherwise, key = value
                """
                namesCorrespondence = {}

                # 4.1) preparing variables, non-arrays
                for varName, specification in arguments["variables"].items():
                    namesCorrespondence[varName] = varName
                    put(f"  let {varName} = {specification['map']};\n")

                # 4.2) preparing arrays that are passed to C-function

                # Dictionary: key - name of argument, value - (type, length, heap, ptr) for JS/TS code 
                arraysRoutine = {}

                for arrName, specification in arguments["arrays"].items():
                    put(f"\n  // Buffer routine: {arrName}\n")

                    if specification["map"] is None: # C-array is not initialized

                        # set length
                        length = "length_" + arrName
                        put(f"  let {length} = {specification['length']};\n")
                        
                        # set number of bytes
                        numOfBytes = "numOfBytes_" + arrName
                        put(f"  let {numOfBytes} = {typesMap[specification['type']]}.BYTES_PER_ELEMENT * {length};\n")

                        # allocate memory
                        dataPtr = "dataPtr_" + arrName
                        put(f"  let {dataPtr} = {moduleName}._malloc({numOfBytes});\n")

                        # dataHeap 
                        dataHeap = "dataHeap_" + arrName
                        put(f"  let {dataHeap} = new Uint8Array({moduleName}.HEAPU8.buffer, {dataPtr}, {numOfBytes});\n")

                        # store type, length, heap and ptr 
                        arraysRoutine[arrName] = (typesMap[specification['type']], length, dataHeap, dataPtr)

                        # set name correspondence
                        namesCorrespondence[arrName] = dataPtr


                    else: # C-array is initialized using some argument of package function

                        packArg = specification['map'] # package function argument

                        # raw data
                        rawData = "rawData_" + arrName

                        if packageFunctionData["input"][packArg]["type"] == "column":
                            put(f"  let {rawData} = {packArg}.getRawData();\n")
                        else:
                            put(f"  let {rawData} = {packArg};\n")  # <- TO DISCUSS OTHER STRUCTURES

                        # set length
                        length = "length_" + arrName
                        put(f"  let {length} = {rawData}.length;\n")                        

                        # set number of bytes
                        numOfBytes = "numOfBytes_" + arrName
                        put(f"  let {numOfBytes} = {rawData}.BYTES_PER_ELEMENT * {length};\n")

                        # allocate memory
                        dataPtr = "dataPtr_" + arrName
                        put(f"  let {dataPtr} = {moduleName}._malloc({numOfBytes});\n")

                        # dataHeap 
                        dataHeap = "dataHeap_" + arrName
                        put(f"  let {dataHeap} = new Uint8Array({moduleName}.HEAPU8.buffer, {dataPtr}, {numOfBytes});\n")

                        # set heap
                        put(f"  {dataHeap}.set(new Uint8Array({rawData}.buffer));\n")

                        # store type, length, heap and ptr 
                        arraysRoutine[arrName] = (typesMap[specification['type']], length, dataHeap, dataPtr)

                        # set name correspondence
                        namesCorrespondence[arrName] = dataPtr


                # 4.3) executing C-function

                # 4.3.1) create a string of arguments, separeted with comma

                stringOfArguments = ""

                for name in functionData["orderOfArguments"]:
                    stringOfArguments += namesCorrespondence[name] + ", "

                # remove redundant ',' at the end
                stringOfArguments = stringOfArguments.rstrip(", ")  

                # 4.3.2) execute C-function 
                put("\n  // Call exported C-function\n")                
                if typeOfOutput != "void":
                    put(f"  let exportedFunctionOutput = {moduleName}._{nameOfFunction}({stringOfArguments});\n")
                else:
                    put(f"  {moduleName}._{nameOfFunction}({stringOfArguments});\n")      


                # 4.4) updating variables & arrays - TO BE DISCUSSED
                for name in packageFunctionData["input"].keys():
                    if packageFunctionData["input"][name]["update"] is not None:
                        # TODO: add modification
                        pass

                # 4.5) creating output structure: column, array ... (if it exists) - TO BE DISCUSSED
                
                # process each output that is a structure
                for name, specification in packageFunctionData["output"].items():                    

                    if specification["type"] == "column":

                        put(f'\n  // Creating output: the {specification["type"]} "{name}"\n')

                        update = specification["update"] # name of C-function argument for updating

                        if update not in arguments["arrays"]:
                            print(f"Error: {update} must be an array!") # <- TODO: handle exception

                        type, length, heap, ignore = arraysRoutine[update]

                        dataArray = "dataArray_" + name
                        put(f"  let {dataArray} = new {type}({heap}.buffer, {heap}.byteOffset, {length});\n")

                        put(f'  let {name} = DG.Column.from{type}("{name}", {dataArray});\n')

                    elif specification["type"] == "array":

                        put(f'\n  // Creating output: the {specification["type"]} "{name}"\n')
                        
                        update = specification["update"] # name of C-function argument for updating

                        if update not in arguments["arrays"]:
                            print(f"Error: {update} must be an array!") # <- TODO: handle exception

                        type, length, heap, ignore = arraysRoutine[update]

                        put(f"  let {name} = new {type}({heap}.buffer, {heap}.byteOffset, {length});\n")

                    elif specification["type"] == "columns":

                        put(f'\n  // Creating output: the {specification["type"]} "{name}"\n')

                        put(f"  let {name} = [];\n")

                        for arrName in specification["update"]:
                            type, length, heap, ignore = arraysRoutine[arrName]
                            put(f"  {name}.push(DG.Column.from{type}('{arrName}', new {type}({heap}.buffer, {heap}.byteOffset, {length})));\n")

                    elif specification["type"] == "arrays":

                        put(f'\n  // Creating output: the {specification["type"]} "{name}"\n')

                        put(f"  let {name} = [];\n")

                        for arrName in specification["update"]:
                            type, length, heap, ignore = arraysRoutine[arrName]
                            put(f"  {name}.push(new {type}({heap}.buffer, {heap}.byteOffset, {length}));\n")


                 

                # 4.6) cleaning allocated memory
                if len(arguments["arrays"]) > 0:

                    put("\n  // Cleaning memory allocated above\n")

                    for name in arguments["arrays"].keys():

                        # 0-th is type, 1-st is length, 2-nd is heap, 3rd is pointer
                        ptr = arraysRoutine[name][3]

                        put(f"  {moduleName}._free({ptr});\n")
                

                # 4.7) returning output (if it exists)
                for name, specification in packageFunctionData["output"].items():

                    # TODO: here, creating of more complex output can be implemented

                    if specification["type"] in structureTypes:  # i.e. this is column or array
                        put(f"\n  return {name};\n")

                    else: # this is a number; so, return is equal to return of C-function
                        put("\n  return exportedFunctionOutput;\n")                                   
                      
                # Finish!
                put("}\n\n")

def appendPackage(settings, functions):
    """ Append Datagrok package file with C-functions.
        settings - dictionary of export settings
        functions - functions data.

        This function opens Datagrok package file and performes the following steps:
          1) add the wrapper for exported C-functions - runtime system;
          2) generate code for exported function call; 
             for this purpose, for each C-file and each exported function the next actions are done:
                2.1) write preambula or annotation of the function;
                2.2) prepare a string of arguments of the function;
                2.3) write the first line of function: <name>(args);
                2.4) write declaration of new arrays (if they are used);
                2.5) write C-function call lines & output lines.

    """
    
    # open package file and append it with C-functions
    with open(settings["packageFile"], 'a') as packageFile:

        put = packageFile.write # a short form of the function that writes data to file
        
        # 1. ADD RUNTIME SYSTEM
        with open(settings["runtimeSystemFile"], 'r') as runtimeSystemFile:
            for line in runtimeSystemFile.readlines():
                put(line)

        # 2. GENERATE CODE FOR CALLING EXPORTED C-FUNCTIONS

        put('\n' * 3 + '// EXPORTED C-FUNCTIONS')

        # consider each exported file
        for cFileName in functions.keys():

            # write comment-line with a name of C-file
            put('\n' * 2 + "// Functions from " + cFileName + '\n\n')

            functionsFromCurrentFile = functions[cFileName]

            # add each exported function
            for nameOfFunction in functionsFromCurrentFile.keys():

                # 2.0. Extract the current function data
                functionData = functionsFromCurrentFile[nameOfFunction]
                arguments = functionData["arguments"] 
                typeOfOutput = functionData["type"]
                packageFunctionData = functionData["packageFunction"]
                cfuncArgs = functionData["orderOfArguments"]
                isExportFunc = packageFunctionData["export"]                

                # 2.1. Write preambula of the current package function
                if isExportFunc:
                    for line in packageFunctionData["preambula"]:
                        put(line + '\n')
                else:
                    put('/*\n')  # in this case, annotation is written in /*...*/
                    for line in packageFunctionData["preambula"]:
                        put('   ' + line.strip("//") + '\n')
                    put("*/\n")

                # 2.2. Prepare string of arguments of the package function
                stringOfArguments = ""

                for argument in packageFunctionData["input"].keys():
                    stringOfArguments += argument + ", "

                # remove redundant ',' at the end
                stringOfArguments = stringOfArguments.rstrip(", ")

                # 2.3. Write the first line of package function
                firstLine = "function " + packageFunctionData["name"] + "(" + stringOfArguments + ") {\n"
                if isExportFunc:
                    put("export " + firstLine)
                else:
                    put(firstLine)

                # 2.4. Write declaration of new arrays (if they are used)
                for nameOfArray, specification in arguments["arrays"].items():
                    if specification["map"] is None:
                        jsType = typesMap[specification["type"]]
                        length =specification["length"]
                        put(f"  let {nameOfArray} = new {jsType}({length});\n")
                
                # 2.5. Write C-function call lines & output lines

                # return types flags
                doExportFuncReturnNumber = (typeOfOutput != 'void')
                doPackageFuncReturnNumber = False 
                doPackageFuncReturnStructure = False      
                
                # 2.5.1. Prepare cFuncWrapper lines
                firstLine = f"{C_FUNC_WRAPPER_NAME}({settings['moduleName']}, '{nameOfFunction}', "
                
                if doExportFuncReturnNumber:
                    firstLine += " 'number',\n"
                else:
                    firstLine += " null,\n"

                # line of argument type in the form like: ['Int32Array', 'number', ...]
                argsTypesLine = '['

                # line of arguments in the form: [arg1, arg2, ...]
                argsLine = '['

                # line of indeces of args that are to update
                argsToUpdateLine = '['

                # the last line(s) of the package function
                returnLine = "  return "

                # line with expression to be returned in the case of column(s) or array(s) return
                expressionToBeReturned = ''

                # complement argsTypesLine & argsLine with the required data defined by arguments
                for arg in cfuncArgs:
                    if arg in arguments["variables"].keys(): # arg is a variable
                        argsTypesLine += "'number', "
                        argsLine += arguments["variables"][arg]["map"] + ", "

                    else: # arg is an array
                        argsTypesLine += f"'{typesMap[arguments['arrays'][arg]['type']]}', "

                        mapStructure = arguments['arrays'][arg]['map']

                        if mapStructure is None: # it's newly created array
                            argsLine += arg + ", "
                        else:
                            typeOfStructure = packageFunctionData["input"][mapStructure]["type"]

                            if typeOfStructure == "column":
                                argsLine += mapStructure + ".getRawData(), "
                            elif typeOfStructure == "array":
                                argsLine += mapStructure + ", "
                          
                # check output & complement argsToUpdateLine with the required indeces
                                
                for name, specification in packageFunctionData["output"].items():
                    outputType = specification["type"]
                    
                    if outputType in simpleStructureTypes: # just one array or column 
                        cArgName = specification["update"]                        
                        argsToUpdateLine += str(cfuncArgs.index(cArgName))
                        doPackageFuncReturnStructure = True

                        if outputType == "array":
                            expressionToBeReturned = cArgName + ";\n"
                        elif outputType == "column":
                            expressionToBeReturned = "DG.Column.from"
                            expressionToBeReturned += typesMap[arguments['arrays'][cArgName]['type']]
                            expressionToBeReturned += f"('{cArgName}', {cArgName});\n"

                    elif outputType in complexStructureTypes: # arrays or columns 
                        expressionToBeReturned = "[" 
                        extraSpace = ' ' * (len(returnLine) + 1)    

                        for cArgName in specification["update"]:
                            argsToUpdateLine += str(cfuncArgs.index(cArgName)) + ', '

                            if outputType == "arrays":
                                expressionToBeReturned += cArgName + ", "
                            elif outputType == "columns":
                                expressionToBeReturned += "DG.Column.from"
                                expressionToBeReturned += typesMap[arguments['arrays'][cArgName]['type']]
                                expressionToBeReturned += f"('{cArgName}', {cArgName}),\n"
                                expressionToBeReturned += extraSpace
                        
                        if outputType == "arrays":
                            expressionToBeReturned = expressionToBeReturned.rstrip(", ") + "];\n"
                        elif outputType == "columns":
                            expressionToBeReturned = expressionToBeReturned.rstrip(",\n" + extraSpace)
                            expressionToBeReturned += "];\n"
                       
                        argsToUpdateLine = argsToUpdateLine.rstrip(', ')  # remove extra ', '
                        doPackageFuncReturnStructure = True                       

                    else: # number
                        doPackageFuncReturnNumber = True             
                              
                # remove extra ', ' (if required) and add the required ending
                argsTypesLine = argsTypesLine.rstrip(', ') + '],\n'
                argsLine = argsLine.rstrip(', ') + '],\n'
                argsToUpdateLine += ']);\n'

                # 2.5.2. Write cFuncWrapper lines
                extraSpace = ' ' # space that makes generated code nicer

                if doPackageFuncReturnNumber and doExportFuncReturnNumber: 
                    extraSpace *= len(returnLine) + len(C_FUNC_WRAPPER_NAME) + 1 
                    put(returnLine + firstLine)  # put "return cFuncWrapper(... "                  
                else: 
                    extraSpace *= len(C_FUNC_WRAPPER_NAME) + 3
                    put(f"  {firstLine}")  # put "cFuncWrapper(... "
                
                put(extraSpace + argsTypesLine)
                put(extraSpace + argsLine)
                put(extraSpace + argsToUpdateLine)
                
                # 2.5.3. Return lines <- in the case of column(s) or array(s) return
                if doPackageFuncReturnStructure:
                    put(returnLine + expressionToBeReturned)
                
                # Finish!
                put("}\n\n")

def getSettings(nameOfFile):
    """
    Return dictionary with settings of C-functions export from json-file with settings.
    """

    with open(nameOfFile, 'r') as file:
         return json.load(file)

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
    command = "emcc" # Emscripten command, starts with "emcc"

    # set optimization mode
    command += " " + settings["optimizationMode"]

    # add names of C-files
    for nameOfFile in settings["source"]:
        command += " " + nameOfFile
    
    # add name of JS-library that is created by Emscripten
    command += " -o " + settings["nameOfLibFile"]

    # set that wasm-file must be created
    command += " -s WASM=1"

    # allow memory growth
    command += " -s ALLOW_MEMORY_GROWTH=1"

    # create module (? should be checked!)
    command += " -s MODULARIZE=1"

    # set export name
    command += ' -s EXPORT_NAME="' + settings["exportName"] + '"'

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

def addDependenciesToPackageJsonFile(settings):
    """ Add JS-file with exported C-functions to "sources" of package.json
    """

    # data from package.json
    packageData = None 

    # get package data from package.json 
    with open(settings["packageJsonFile"], 'r') as file:
        packageData = json.load(file)

    # create full name of JS-lib file that is added to package.json
    fullNameOfLibFile = settings["folder"] + "/" + settings["nameOfLibFile"]

    # add dependence to package data
    if "sources" in packageData.keys():

        if fullNameOfLibFile not in packageData["sources"]: # add JS-file if "source" does not contain it yet
            packageData["sources"].append(fullNameOfLibFile)

    else:
        packageData["sources"] = [ fullNameOfLibFile ]

    # update package.json
    with open(settings["packageJsonFile"], 'w') as file:
        json.dump(packageData, file)

def addModuleInitFunctionToJavaScriptLibFile(settings):
    """
    Add init-function to JS-file of the library.

    The following code is added:

    var [module name] = undefined;

    async function init[module name]() {
      if ([module name] === undefined) {
        console.log("Wasm not Loaded, Loading");
        [module name]  = await [export name]();
      } else {
          console.log("Wasm Loaded, Passing");
        }
    }

    """

    with open(settings["nameOfLibFile"], 'a') as file:

        moduleName = settings["moduleName"]
        exportName = settings["exportName"]

        # add declaration of module
        file.write("\n\nvar " + moduleName + " = undefined;\n\n")

        # add init-function - this approach is copied from the package DSP
        file.write("async function init" + moduleName + "() {\n")
        file.write("  if (" + moduleName + " === undefined) {\n")
        file.write('    console.log("Wasm not Loaded, Loading");\n')
        file.write("    " + moduleName + "  = await " + exportName + "();\n")
        file.write("  } else {\n")
        file.write('    console.log("Wasm Loaded, Passing");\n')
        file.write('  }\n')
        file.write('}')


def main(nameOfSettingsFile="module.json"):
    """
    The main script: 
        1) load export settings from json-file;
        2) parse C-files and get C-functions to be exported;
        3) save C-functions descriptors (name, type of return value, arguments descriptors: type, name);
        4) create Emscripten command;
        5) save Emscripten command to text file;
        6) execute Emscripten command;
        7) append Datagrok package file with exported functions;
        8) add dependecnies to the file package.json;
        9) add module init-function to JS-file with exported C-functions. 
    """
    try:
        # load settings
        settings = getSettings(nameOfSettingsFile)

        # parse C-files and get exported functions data
        functionsData = getFunctionsFromListOfFiles(settings["source"])
    
        # write functions descriptors to json-file
        saveFunctionsDataToFile(functionsData, settings["nameOfFileForFunctionsData"])

        # get command for Emscripten tool
        command = getCommand(settings, functionsData)  

        # save command to file
        saveCommand(command, settings["fileWithEmscriptenCommand"])  

        # execute command by Emscripten (IT MUST BE INSTALLED!)   
        os.system(command)
      
        # append Datagrok package file with exported functions: old style
        #appendPackageFileWithCfunctions(settings["packageFile"], settings["moduleName"], functionsData)

        # append Datagrok package file with exported functions: runtime system use
        appendPackage(settings, functionsData)

        # add dependencies to the file package.json
        addDependenciesToPackageJsonFile(settings)

        # add module init-function to JS-file of the library
        addModuleInitFunctionToJavaScriptLibFile(settings) 

    except OSError:
        print("OS ERROR: check paths and file names!")
    
    except KeyError:
        print("PARSING ERROR: check preambula!")

    except IndexError:
        print("PARSING ERROR: check arguments")


if __name__ == '__main__':
    main()