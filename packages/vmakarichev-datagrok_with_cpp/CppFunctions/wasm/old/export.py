""" export.py
    
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

def getFunctionsFromFile(nameOfFile):
    """Process the C-file and returns dictionary of functions, type of return, types of arguments.
       nameOfFile - name of C-file.
    """
    
    with open(nameOfFile, 'r') as file:

        listOfLines = file.readlines()

        numOfLines = len(listOfLines)

        i = 0 # index of the current line

        dictionaryOfFunctions = {}

        # analyse all lines of file
        while i < numOfLines:

            if listOfLines[i] == "EMSCRIPTEN_KEEPALIVE\n":  # the next line is function's prototype
                functionPrototype = listOfLines[i + 1]
                #print("Prototype: ", functionPrototype, end='')

                # detection of function descriptor and a list of arguments
                functionDescriptor, ignore, rest = functionPrototype.partition("(")
                arguments, ignore, trash = rest.partition(")")

                typeOfReturnValue, ignore, nameOfFunction = functionDescriptor.rpartition(" ")
                #print(f"returnType: {typeOfReturnValue}, name of function: {nameOfFunction}")

                listOfArgumnetDescriptors = arguments.split(',') 

                """print("  descriptor: ", functionDescriptor)
                print(" arguments: ", listOfArgumnetDescriptors)"""

                listOfArgumentsData = []

                for argument in listOfArgumnetDescriptors:
                    #print("    ", argument.rpartition(" "))
                    typeDescriptor, ignore, nameOfArgument = argument.rpartition(" ")
                    typeDescriptor = typeDescriptor.strip(" ")
                    #print(f"type: {typeDescriptor}, name: {nameOfArgument}")
                    listOfArgumentsData.append([typeDescriptor, nameOfArgument])

                #print()

                dictionaryOfFunctions[nameOfFunction] = {"type": typeOfReturnValue, "arguments": listOfArgumentsData}

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
    """

    # open file and append it with C-functions
    with open(nameOfFile, 'a') as file:

        file.write('\n' * 3 + '// C-FUNCTIONS')

        # consider each exported file
        for cFileName in functions.keys():

            # write comment-line with a name of C-file
            file.write('\n' * 2 + "// Functions from " + cFileName + '\n')

            functionsFromCurrentFile = functions[cFileName]

            # add each exported function
            for nameOfFunction in functionsFromCurrentFile.keys():

                # 1. WRITE ANNOTATION OF THE CURRENT FUNCTION

                file.write("\n/* Exported C-function\n")

                # write name of function
                file.write(f"name: {nameOfFunction}\n")

                functionData = functionsFromCurrentFile[nameOfFunction]

                listOfArguments = functionData["arguments"]

                typeOfOutput = functionData["type"]

                stringOfArguments = ""

                # write each input value data
                for argument in listOfArguments:
                    file.write(f"input: {argument[0]} {argument[1]}\n")
                    stringOfArguments += argument[1] + ", " # complete a string with arguments

                # remove redundant ',' at the end
                stringOfArguments = stringOfArguments.rstrip(", ")
                
                # write output type specification
                file.write(f"output: {typeOfOutput} result  */\n")       


                # 2. WRITE NAME OF THE CURRENT FUNCTION         

                # prepare the first line of the function: async case - old approach
                #firstLineOfFunction = "async function " + nameOfFunction + "(" + stringOfArguments + ") { // <--- CHECK THIS!\n"

                # prepare the first line of the function
                firstLineOfFunction = "function " + nameOfFunction + "(" + stringOfArguments + ") { // <--- CHECK THIS!\n"


                # 3. BODY OF THE FUNCTION

                # write the first line
                file.write(firstLineOfFunction)

                # write module initialization line: async case - old approach
                #file.write("  await init" + moduleName + "();\n")


                # 4. Extract arrays from all arguments & create execution string & create execution string

                stringOfArguments = '' # this is a string of arguments that will be put into exported functions

                arraysAttributes = [] # attributes of arrays: type, ...

                # analyze each argument
                for argument in listOfArguments:
                    if "*" in argument[0]:
                        type = argument[0]
                        name = argument[1]
                        numOfBytes = "numOfBytes_" + name
                        length = name + ".length"
                        bytesPerEl = name + ".BYTES_PER_ELEMENT"
                        dataPtr = "dataPtr_" + name
                        dataHeap = "dataHeap_" + name

                        stringOfArguments += dataPtr + ", "

                        file.write("\n  // Buffer routine: " + name + ". Check it!\n")

                        currentString = f"  let {numOfBytes} = {length} * {bytesPerEl}; // <--- CHECK THIS!\n"
                        file.write(currentString)

                        currentString = f"  let {dataPtr} = {moduleName}._malloc({numOfBytes});\n"
                        file.write(currentString)

                        currentString = f"  let {dataHeap} = new Uint8Array({moduleName}.HEAPU8.buffer, {dataPtr}, {numOfBytes});\n"
                        file.write(currentString)

                        currentString = f"  {dataHeap}.set(new Uint8Array({name}.buffer)); // <--- CHECK THIS!\n"
                        file.write(currentString)

                        arraysAttributes.append((type, name, numOfBytes, length, bytesPerEl, dataPtr, dataHeap))
                    else:
                        stringOfArguments += argument[1] + ', '

                # remove redundant ',' at the end
                stringOfArguments = stringOfArguments.rstrip(", ")


                # 5. Write line with execution of exported function
                if typeOfOutput != "void":
                    file.write("\n  let result = " + moduleName + "._" + nameOfFunction + "(" + stringOfArguments + ");\n")
                else:
                    file.write("\n  " + moduleName + "._" + nameOfFunction + "(" + stringOfArguments + ");\n")


                # 6. Write lines: copying data from heap to JS-arrays - optional

                if len(arraysAttributes) > 0:
                    file.write("\n  /* OPTIONAL: each array could be modified when executing exported C-function.")
                    file.write("\n     Remove '//' in order to get the modified array.")
                    file.write("\n     The modified array further can be returned.")
                    file.write("\n     MODIFY if it is required! */\n\n")

                for arrayData in arraysAttributes:
                    type, name, numOfBytes, length, bytesPerEl, dataPtr, dataHeap = arrayData

                    jsTypeName = "?????????"
                    if type in typesMap.keys():
                        jsTypeName = typesMap[type]


                    file.write(f"  //let {name}Modified = new {jsTypeName}({dataHeap}.buffer, ")
                    file.write(f"{dataHeap}.byteOffset, {length});\n")


                # 7. Write lines with memory cleaning

                if len(arraysAttributes) > 0:
                    file.write("\n  // Cleaning allocated memory.\n")

                for arrayData in arraysAttributes:
                    type, name, numOfBytes, length, bytesPerEl, dataPtr, dataHeap = arrayData

                    currentString = f"  {moduleName}._free({dataPtr});\n"
                    file.write(currentString)


                # 8. Finally, the last row of the function: "return ... "
                if typeOfOutput == "void":
                    file.write("\n//  return ...;   //  <--- here, smth can be returned\n}\n")
                else:
                    file.write("\n  return result;\n}\n")

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
            print("OS error!")

if __name__ == '__main__':
    main()