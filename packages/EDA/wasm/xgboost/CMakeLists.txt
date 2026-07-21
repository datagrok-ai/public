# Emscripten build of the XGBoost wasm module for the Datagrok EDA package.
#
# XGBoost sources are NOT stored in this repo: pass the path to a clone
# (with the minimal-build patch applied, see patches/) via XGBOOST_SOURCE_DIR;
# the built artifacts (XGBoostAPI.js + .wasm) are copied into wasm/ by build.ps1.
# Manual use:
#   emcmake cmake -B build -DXGBOOST_SOURCE_DIR=<path-to-xgboost>
#   cmake --build build --parallel
cmake_minimum_required(VERSION 3.18)
project(XGBoostAPI CXX)

set(CMAKE_CXX_STANDARD 17)
if(NOT DEFINED EMSCRIPTEN)
  message(FATAL_ERROR "This project is wasm-only; configure with emcmake.")
endif()

set(XGBOOST_SOURCE_DIR "" CACHE PATH "Path to the xgboost source tree (patched)")
if(NOT XGBOOST_SOURCE_DIR)
  message(FATAL_ERROR "Pass -DXGBOOST_SOURCE_DIR=<path to xgboost clone>")
endif()

# Speed-first flags decided by benchmark (see README.md here):
# -O3 + wasm SIMD + LTO is the fastest fit AND predict variant; SIMD without
# LTO regresses fit on wide data, -Os makes predict ~5x slower.
# C++ exceptions stay disabled (Emscripten default) - all validation happens
# on the TypeScript side. Overridable for flag benchmarking (smoke/bench.mjs).
set(WASM_OPT_FLAGS "-O3 -msimd128 -flto" CACHE STRING "Optimization flags for compile and link")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WASM_OPT_FLAGS}")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Minimal xgboost: gbtree + hist, reg/binary/multi objectives only.
set(XGBOOST_WASM_MINIMAL ON CACHE BOOL "" FORCE)
set(BUILD_STATIC_LIB ON CACHE BOOL "" FORCE)
set(USE_OPENMP OFF CACHE BOOL "" FORCE)
set(KEEP_BUILD_ARTIFACTS_IN_BINARY_DIR ON CACHE BOOL "" FORCE)
add_subdirectory(${XGBOOST_SOURCE_DIR} xgboost-build)

# Classic MODULARIZE (no EXPORT_ES6): the generated UMD tail makes the same
# file importable by webpack (CommonJS) in both the main bundle and the
# worker bundle, and require()-able under node for tests.
#
# Two link targets from the same static library:
#  - XGBoostAPI      (web,worker)      - the shipped artifact. Must NOT
#    include the node environment: its loader would contain
#    require("node:fs") calls that break webpack bundling.
#  - XGBoostAPI-node (web,worker,node) - used only by smoke/api-test.mjs
#    in build.ps1; never copied to wasm/.
#
# GROWABLE_ARRAYBUFFERS defaults to 1 in Emscripten 6.x: the loader then backs
# the heap with a resizable ArrayBuffer (wasmMemory.toResizableBuffer) in new
# browsers, and TextDecoder.decode over such a buffer throws per spec
# ("The provided ArrayBuffer value must not be resizable") on the first string
# XGBoost logs. Node lacks toResizableBuffer, so tests pass while real
# browsers crash. Force the classic grow-and-replace behavior.
set(XGB_COMMON_LINK_FLAGS "${WASM_OPT_FLAGS} \
-sWASM=1 \
-sALLOW_MEMORY_GROWTH=1 \
-sGROWABLE_ARRAYBUFFERS=0 \
-sINITIAL_MEMORY=33554432 \
-sMODULARIZE=1 \
-sEXPORT_NAME=XGBoost \
-sEXPORTED_FUNCTIONS=_xgbTrain,_xgbModelSize,_xgbModelCopy,_xgbLoadModel,_xgbPredict,_xgbFreeModel,_malloc,_free \
-sEXPORTED_RUNTIME_METHODS=HEAPF32,HEAPU8")

foreach(variant IN ITEMS XGBoostAPI XGBoostAPI-node)
  add_executable(${variant} xgboost-api.cpp)
  target_include_directories(${variant} PRIVATE
    ${XGBOOST_SOURCE_DIR}/include
    ${XGBOOST_SOURCE_DIR}/dmlc-core/include)
  target_link_libraries(${variant} PRIVATE xgboost)
endforeach()

set_target_properties(XGBoostAPI PROPERTIES
  LINK_FLAGS "${XGB_COMMON_LINK_FLAGS} -sENVIRONMENT=web,worker")
set_target_properties(XGBoostAPI-node PROPERTIES
  LINK_FLAGS "${XGB_COMMON_LINK_FLAGS} -sENVIRONMENT=web,worker,node")
