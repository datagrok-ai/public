call setup-unlink-clean.cmd

set package_dir=%cd%

set dirs=^
\..\..\js-api\ ^
\..\..\libraries\utils\ ^
\..\..\libraries\bio\ ^
\

for %%p in (%dirs%) do cd %package_dir%\%%p & call npm install
for %%p in (%dirs%) do cd %package_dir%\%%p & call npm link
for %%p in (%dirs%) do cd %package_dir%\%%p & call npm run link-all
for %%p in (%dirs%) do cd %package_dir%\%%p & call npm run build
