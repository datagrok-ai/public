lib1.c, lib2.c - C-files with functions to be exported to Datagrok package

test.cpp, test.exe - C++ tools for testing lib1.c, lib2.c

listOfFilesToAdd.txt - text file with a list of C-libs to be exported

export.py - the main script (basic version)
 the following actions are performed:
  1) names of C-libs are read;
  2) C-files are parsed and the desired functions data is extracted;
  3) the file package.js is appended with JS-versions of C-functions;
  4) the file functionsData.json, which contains parsed functions data, is created.

