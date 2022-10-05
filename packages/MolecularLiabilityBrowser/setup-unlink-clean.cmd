set package_dir=%cd%

set dirs=^
\..\..\js-api\ ^
\..\..\libraries\ml\ ^
\..\..\libraries\utils\ ^
\..\..\libraries\bio\ ^
\..\..\packages\Bio\ ^
\..\..\packages\Helm\ ^
\..\..\packages\MolecularLiabilityBrowserData ^
\

call npm uninstall -g ^
    datagrok-api ^
    @datagrok-libraries/ml ^
    @datagrok-libraries/utils ^
    @datagrok-libraries/bio ^
    @datagrok/bio ^
    @datagrok/helm ^
    @datagrok/molecular-liability-browser-data

for %%p in (%dirs%) do cd %package_dir%\%%p & rmdir /s /q node_modules
for %%p in (%dirs%) do cd %package_dir%\%%p & rmdir /s /q dist

for %%p in (%dirs%) do cd %package_dir%\%%p & del "package-lock.json"
