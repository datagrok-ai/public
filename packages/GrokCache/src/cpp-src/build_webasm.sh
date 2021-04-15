#!/bin/sh
#
# Utility script to rebuild WASM module
#
# ------------------------------------
# -s ENVIRONMENT="web"
#  - needed NOT to generate node.js elemenst in .js file, maing WASM portable
#    between browser and Node.js
#  - if not used webpack may need package 'fs'
#
# package.json:
#  "dependencies": {
#    "fs": "^0.0.1-security"    // <--- !!
#  }
# ------------------------------------
#
# -s EXTRA_EXPORTED_RUNTIME_METHODS=['ccall']
#   not used now, but in some sources considered necessary
#

emcc -std=c++11 -O2 --bind -s "MODULARIZE=1" -s "EXPORT_NAME=createGrokCache"  -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' -s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s DISABLE_EXCEPTION_CATCHING=0 -s NO_EXIT_RUNTIME=1 -s ENVIRONMENT="web"  grokcache.cpp -o grokcache.js
