@echo off
REM Builds sampler.js + sampler.wasm: the OptGP flux sampler (sampler.cpp) with GLPK 5.0 compiled in.
REM Run it from an emsdk shell (em++ on PATH, e.g. after emsdk_env.bat).
REM GLPK 5.0 sources: https://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz, extracted to GLPK_DIR (default below).
REM GLPK is compiled once into %GLPK_DIR%/build-wasm; delete that folder to rebuild it.
REM -ffp-contract=off keeps the arithmetic identical to cobra's (no fused multiply-adds).
setlocal EnableDelayedExpansion
cd /d "%~dp0"
if "%GLPK_DIR%"=="" set "GLPK_DIR=C:/Users/rizhi/Desktop/experiments/cpp/glpk-5.0"
set "GLPK_DIR=%GLPK_DIR:\=/%"
where em++ >nul 2>nul || (echo em++ not found: run emsdk_env.bat first & exit /b 1)
if not exist "%GLPK_DIR%/src/glpk.h" (echo GLPK 5.0 sources not found in %GLPK_DIR% & exit /b 1)

set "OBJ=%GLPK_DIR%/build-wasm"
REM cmd built-ins (del, rmdir, mkdir, pushd) need backslash paths
set "OBJW=%OBJ:/=\%"
if not exist "%OBJ%/glpk.done" (
    echo Compiling GLPK...
    if exist "%OBJW%" rmdir /s /q "%OBJW%"
    mkdir "%OBJW%"
    echo /* GLPK configuration for a single-threaded WebAssembly build */> "%OBJ%/config.h"
    for /r "%GLPK_DIR%/src" %%f in (*.c) do (
        set "p=%%~ff"
        REM src/proxy/main.c is a standalone program, not part of the library
        if /i not "%%~nxf"=="main.c" echo "!p:\=/!">> "%OBJ%/sources.rsp"
    )
    set "INC=-I%OBJ% -I%GLPK_DIR%/src"
    for %%d in (amd api bflib colamd draft env intopt minisat misc mpl npp proxy simplex zlib) do set "INC=!INC! -I%GLPK_DIR%/src/%%d"
    pushd "%OBJW%"
    call emcc -std=gnu11 -O2 -ffp-contract=off -DHAVE_CONFIG_H !INC! -c @sources.rsp -Wno-everything || (popd & exit /b 1)
    popd
    echo done> "%OBJ%/glpk.done"
)
if exist "%OBJW%\objects.rsp" del "%OBJW%\objects.rsp"
for %%f in ("%OBJ%/*.o") do echo "%OBJ%/%%~nxf">> "%OBJ%/objects.rsp"

echo Linking sampler...
call em++ -O3 -std=c++17 -ffp-contract=off -I"%GLPK_DIR%/src" sampler.cpp @"%OBJ%/objects.rsp" -o sampler.js ^
    -sMODULARIZE=1 -sEXPORT_NAME=exportCppSampler -sENVIRONMENT=web,worker ^
    -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=67108864 -sMAXIMUM_MEMORY=2147483648 -sSTACK_SIZE=1048576 ^
    -sEXPORTED_FUNCTIONS=_sample,_sample_with_warmup,_malloc,_free ^
    -sEXPORTED_RUNTIME_METHODS=cwrap,HEAPF64,HEAPF32,HEAP32 || exit /b 1

call "%~dp0prepend.bat"
echo Built sampler.js and sampler.wasm
