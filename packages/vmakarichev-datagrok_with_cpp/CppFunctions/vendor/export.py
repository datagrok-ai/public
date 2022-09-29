# parser.py

import os
import json

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

                dictionaryOfFunctions[nameOfFunction] = {"type": typeDescriptor, "arguments": listOfArgumentsData}

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

        for cFileName in functions.keys():
            file.write('\n' * 2 + "// Functions from " + cFileName + '\n')

            functionsFromCurrentFile = functions[cFileName]

            for nameOfFunction in functionsFromCurrentFile.keys():
                file.write(f"\n//name: {nameOfFunction}\n")

                functionData = functionsFromCurrentFile[nameOfFunction]

                listOfArguments = functionData["arguments"]

                typeOfOutput = functionData["type"]

                stringOfArguments = ""

                for argument in listOfArguments:
                    file.write(f"//input: {argument[0]} {argument[1]}\n")
                    stringOfArguments += argument[1] + ", "

                stringOfArguments = stringOfArguments.rstrip(", ")
                
                file.write(f"//output: {typeOfOutput} result\n")                

                firstLineOfFunctions = "export async function " + nameOfFunction + "(" + stringOfArguments + ") {\n"
                file.write(firstLineOfFunctions)

                file.write("  await init" + moduleName + "();\n")

                file.write("  return " + moduleName + "._" + nameOfFunction + "(" + stringOfArguments + ");\n}\n")

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
    
     # !!!! May be some other functions should also be added
    
    # delete ',' at the end of the command
    command = command.rstrip(',') 

    # add closing bracket
    command += "]" 

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