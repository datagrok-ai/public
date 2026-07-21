# Native smoke test for the minimal (XGBOOST_WASM_MINIMAL) build of XGBoost.
#
# Prerequisite: a native build of the patched xgboost tree, e.g. (in the clone)
#   cmake -B build-minimal -G "Visual Studio 17 2022" -A x64 \
#         -DXGBOOST_WASM_MINIMAL=ON -DBUILD_STATIC_LIB=ON
#   cmake --build build-minimal --config Release --parallel
#
# Then:
#   cmake -B build -G "Visual Studio 17 2022" -A x64 -DXGB=<path-to-xgboost-clone>
#   cmake --build build --config Release
#   build/Release/smoke.exe        # expects: SMOKE TEST PASSED
cmake_minimum_required(VERSION 3.18)
project(xgb_smoke CXX)
set(CMAKE_CXX_STANDARD 17)

set(XGB "" CACHE PATH "Path to the xgboost clone with a build-minimal/ native build")
if(NOT XGB)
  message(FATAL_ERROR "Pass -DXGB=<path to xgboost clone>")
endif()

# xgboost is built with static CRT (/MT) by default; match it (set before the target)
cmake_policy(SET CMP0091 NEW)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

add_executable(smoke smoke.cpp)
target_include_directories(smoke PRIVATE ${XGB}/include ${XGB}/dmlc-core/include)
target_link_libraries(smoke PRIVATE
  ${XGB}/lib/xgboost.lib
  ${XGB}/build-minimal/src/objxgboost.dir/Release/objxgboost.lib
  ${XGB}/build-minimal/dmlc-core/Release/dmlc.lib)

# The native verification build keeps OpenMP on (the wasm build uses USE_OPENMP=OFF)
find_package(OpenMP REQUIRED CXX)
target_link_libraries(smoke PRIVATE OpenMP::OpenMP_CXX)
