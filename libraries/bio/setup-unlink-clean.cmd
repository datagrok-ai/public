set package_dir=%cd%

set dirs=^
\..\..\js-api\ ^
\..\..\libraries\utils\ ^
\..\..\libraries\gridext\ ^
\

call npm uninstall -g datagrok-api @datagrok-libraries/utils @datagrok-libraries/gridext

for %%p in (%dirs%) do cd %package_dir%\%%p & rmdir /s /q node_modules
for %%p in (%dirs%) do cd %package_dir%\%%p & rmdir /s /q dist

rem for %%p in (%dirs%) do cd %package_dir%\%%p & del "package-lock.json"
