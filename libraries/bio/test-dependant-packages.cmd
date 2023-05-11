@echo on
call ./setup-unlink-clean.sh

set base_dir=%cd%

set libDirs=^
\..\..\js-api\ ^
\

set packageDirs=^
\..\..\packages\Bio\ ^
\..\..\packages\BiostructureViewer\ ^
\..\..\packages\Dendrogram\ ^
\..\..\packages\Peptides\ ^
\..\..\packages\PhyloTreeViewer\
rem \..\..\packages\SequenceTranslator\
rem \..\..\packages\Helm

for %%p in (%libDirs%) do cd %base_dir%\%%p & call npm install & if %errorlevel% neq 0 cmd /k
for %%p in (%libDirs%) do cd %base_dir%\%%p & call npm link datagrok-api & if %errorlevel% neq 0 cmd /k
for %%p in (%libDirs%) do cd %base_dir%\%%p & call npm run build & if %errorlevel% neq 0 cmd /k
for %%p in (%libDirs%) do cd %base_dir%\%%p & call npm link & if %errorlevel% neq 0 cmd /k

for %%p in (%packageDirs%) do cd %base_dir%\%%p & call npm install & if %errorlevel% neq 0 cmd /k
for %%p in (%packageDirs%) do cd %base_dir%\%%p & call npm link datagrok-api @datagrok-libraries/bio & if %errorlevel% neq 0 cmd /k
for %%p in (%packageDirs%) do cd %base_dir%\%%p & call npm run build & if %errorlevel% neq 0 cmd /k
for %%p in (%packageDirs%) do cd %base_dir%\%%p & call grok test --skip-build --host local & if %errorlevel% neq 0 cmd /k
