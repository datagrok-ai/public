@echo off
setlocal

REM Save the current and escher directories
set PACKAGE_DIR=%cd%
set "ESCHER_DIR=..\Escher"
set "PACKAGE_ESCHER_DIST_DIR=%PACKAGE_DIR%\escher_dist"

REM Navigate to the escher directory
cd /d "%ESCHER_DIR%"

REM Run npm build
echo Building escher...
set NODE_OPTIONS=--openssl-legacy-provider
call npm run build

REM Check if the build was successful by looking for escher.js in the dist directory
if exist dist\escher.js (
    echo Build successful. Copying files...

    REM Copy files back to the starting location
    copy /y dist\escher.js "%PACKAGE_ESCHER_DIST_DIR%"
    copy /y dist\escher.js.map "%PACKAGE_ESCHER_DIST_DIR%"

    echo Files copied successfully to %PACKAGE_ESCHER_DIST_DIR%.
) else (
    echo Build failed or dist\escher.js not found.
)

REM Return to the starting directory
cd /d "%PACKAGE_DIR%"

REM Build the package
echo Building package...
call npm run build

REM Publish the package
echo Publishing package...
call grok publish

endlocal