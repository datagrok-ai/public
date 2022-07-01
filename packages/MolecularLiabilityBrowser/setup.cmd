set package_dir=%cd%

call unlink.cmd
pause

cd %package_dir%/../../js-api
call npm install
call npm link

cd %package_dir%/../../libraries/utils
call npm install
call npm run link-all
call npm link

cd %package_dir%/../../libraries/bio
call npm install
call npm run link-all
call npm link

cd %package_dir%/../../packages/Bio
call npm install
call npm run link-all
call npm link

cd %package_dir%/../../packages/MolecularLiabilityBrowserData
call npm install
call npm run link-all
call npm link

cd %package_dir%
call npm install
call npm run link-all

pause

echo js-api ...
cd %package_dir%/../../js-api
call npm run build

echo Utils ...
cd %package_dir%/../../libraries/utils
call npm run build

echo bio ...
cd %package_dir%/../../libraries/bio
call npm run build

echo Bio ...
cd %package_dir%/../../packages/Bio
call npm run build

echo MLB-Data ...
cd %package_dir%/../../packages/MolecularLiabilityBrowserData
call npm run build

echo MLB ...
cd %package_dir%
call npm run build
