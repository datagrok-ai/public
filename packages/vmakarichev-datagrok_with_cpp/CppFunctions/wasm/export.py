""" export.py  Ocotber 12, 2022

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

# a set of types that represents data structures
structureTypes={"column", "array"}

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

                """ List of arguments of the function.
                    This list is further used, when executing C-function in package 
                """
                orderOfArguments = []

                # process each argument descriptor: get name, type, ...
                for argument in listOfArgumentDescriptors:

                    typeDescriptor, ignore, nameOfArgument = argument.rpartition(" ")
                    typeDescriptor = typeDescriptor.strip(" ")

                    # store name of the current argument
                    orderOfArguments.append(nameOfArgument)

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
                packageFunctionData = {"name": None, 
                "input": inputs, "output": outputs, "preambula": None }
 
                mapInstructions = []  # lines that contains 'map'

                updateInstructions = []  # lines that contains 'update'

                # consider each line that starts with '//'
                while listOfLines[j].startswith("//"): 

                    line = listOfLines[j].rstrip('\n') # current string without '\n'

                    if "map" in line:
                        mapInstructions.append(line)
                    elif "update" in line:
                        updateInstructions.append(line)
                    else:
                        packageFunctionPreambula.append(line)

                    j -= 1


                # process package function's preambula: reversed pass is used
                for line in packageFunctionPreambula[::-1]:

                    if line.startswith("//name:"):
                        ignore1, ignore2, nameOfFunc = line.partition(": ")
                        packageFunctionData["name"] = nameOfFunc

                    elif line.startswith("//input:"):
                        ignore1, ignore2, argsData = line.partition(": ")
                        argType, ignore1, argName = argsData.partition(" ")
                        inputs[argName] = {"type": argType, "update": None}

                    elif line.startswith("//output:"):
                        ignore1, ignore2, argsData = line.partition(": ")
                        argType, ignore1, argName = argsData.partition(" ")
                        outputs[argName] = {"type": argType}


                # process map-instructions
                for line in mapInstructions:
                    ignore1, ignore2, argMap = line.partition(": ")                        
                    cFuncArg, ignore1, packageFuncArg = argMap.partition(" is ")

                    if cFuncArg in variablesData.keys(): # argument is a variable, NOT an array
                        variablesData[cFuncArg]["map"] = packageFuncArg

                    else: # argument is an array
                        if "new" in packageFuncArg:

                            # Extract an expression in (...) after the word 'new'
                            ignore1, ignore2, important = packageFuncArg.partition('(')
                            lengthValue, ignore1, ignore2 = important.partition(')')
                            arraysData[cFuncArg]["length"] = lengthValue

                        else:
                            arraysData[cFuncArg]["map"] = packageFuncArg             
                

                # process 'update'-instructions
                for line in updateInstructions:
                    ignore1, ignore2, argMap = line.partition(": ")
                    packageFuncArg, ignore1, cFuncArg = argMap.partition(" is ")

                    if packageFuncArg in inputs.keys():
                        inputs[packageFuncArg]["update"] = cFuncArg
                    else:
                        outputs[packageFuncArg]["update"] = cFuncArg
     

                # add package function preambula: //name: ..., etc.
                packageFunctionData["preambula"] = packageFunctionPreambula[::-1]

                dictionaryOfFunctions[nameOfFunction] = {"type": typeOfReturnValue, 
                "arguments": {"variables": variablesData, "arrays": arraysData},
                "orderOfArguments": orderOfArguments,
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

                # 1. Write preambula of the current package function
                for line in packageFunctionData["preambula"]:
                    put(line + '\n')

                # 2. Prepare string of arguments of the package function
                stringOfArguments = ""

                for argument in packageFunctionData["input"].keys():
                    stringOfArguments += argument + ", "

                # remove redundant ',' at the end
                stringOfArguments = stringOfArguments.rstrip(", ")


                # 3. Write the first line of package function
                put("export function " + packageFunctionData["name"] + "(" + stringOfArguments + ") {\n")


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

                    elif specification["type"] == "array": # <- currently, just array is considered; TODO

                        put(f'\n  // Creating output: the {specification["type"]} "{name}"\n')
                        
                        update = specification["update"] # name of C-function argument for updating

                        if update not in arguments["arrays"]:
                            print(f"Error: {update} must be an array!") # <- TODO: handle exception

                        type, length, heap, ignore = arraysRoutine[update]

                        put(f"  let {name} = new {type}({heap}.buffer, {heap}.byteOffset, {length});\n")
                 

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
    command += ' -s EXPORTED_RUNTIME_METHODS=["cwrap"]'  # also, "ccall" can be added 

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
      
        # append Datagrok package file with exported functions
        appendPackageFileWithCfunctions(settings["packageFile"], settings["moduleName"], functionsData)

        # add dependencies to the file package.json
        addDependenciesToPackageJsonFile(settings)

        # add module init-function to JS-file of the library
        addModuleInitFunctionToJavaScriptLibFile(settings) 

    except OSError:
        print("OS ERROR: check paths and file names!")
    
    except KeyError:
        print("PARSING ERROR: check preambula!")


if __name__ == '__main__':
    main()