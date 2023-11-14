set package_dir=%cd%

set dirs=^
\..\..\js-api\ ^
\..\..\libraries\utils\ ^
\..\..\libraries\statistics\ ^
\..\..\libraries\tutorials\ ^
\..\..\libraries\ml\ ^
\..\..\libraries\bio\ ^
\..\..\packages\Bio\ ^
\..\..\packages\Helm\ ^
\..\..\packages\Dendrogram\ ^
\..\..\packages\Chem\ ^
\

call npm uninstall -g ^
    datagrok-api ^
    @datagrok-libraries/utils ^
    @datagrok-libraries/statistics ^
    @datagrok-libraries/tutorials ^
    @datagrok-libraries/ml ^
    @datagrok-libraries/bio ^

for %%p in (%dirs%) do cd %package_dir%\%%p & rmdir /s /q node_modules
for %%p in (%dirs%) do cd %package_dir%\%%p & rmdir /s /q dist

for %%p in (%dirs%) do cd %package_dir%\%%p & del "package-lock.json"
