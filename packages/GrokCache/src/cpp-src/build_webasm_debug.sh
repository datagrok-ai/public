#!/bin/sh
#
# Utility script to rebuild WASM module
#
emcc -std=c++11 -O0 --bind -s "MODULARIZE=1" -s "EXPORT_NAME=createGrokCache" -s EXTRA_EXPORTED_RUNTIME_METHODS=['ccall'] -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' -s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s DISABLE_EXCEPTION_CATCHING=0 -s NO_EXIT_RUNTIME=1 -s ENVIRONMENT="web" grokcache.cpp -o grokcache.js
