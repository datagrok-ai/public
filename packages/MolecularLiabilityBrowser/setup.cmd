call setup-unlink-clean.cmd

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

for %%p in (%dirs%) do cd %package_dir%\%%p & call npm install
for %%p in (%dirs%) do cd %package_dir%\%%p & call npm link
for %%p in (%dirs%) do cd %package_dir%\%%p & call npm run link-all
for %%p in (%dirs%) do cd %package_dir%\%%p & call npm run build
