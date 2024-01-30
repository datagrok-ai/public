# Web Assembly

Package contains webAssembly code (compiled from cpp) of hierarchical clustering using reduced distance matrix.
the files are contained in ```wasm``` folder and corresponding cpp files are in ```wasm/cpp```. to modify the cpp code and compile it, following steps are to be done. 

1. Emscripten has to be installed on system and found in global PATH

2. In command line prompt, cd into the ```wasm/cpp``` directory and execute command ```emsdk activate latest``` .  This will init the working directory.

3. Then execute command found in same folder under the name ```emsc command.txt```.

You will generate two files, ```wasmCluster.js``` and ```wasmCluster.wasm```. In order to allow the wasm to work with webworkers, following has to be done:

1. Move the ```wasmCluster.js``` and ```wasmCluster.wasm``` to ```wasm``` directory.

2. Open ```wasmCluster.js``` and before the second line ```var exportCppLib = (() => {``` add ```export```

Done.
