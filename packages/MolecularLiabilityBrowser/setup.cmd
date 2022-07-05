set package_dir=%cd%

call setup-unlink-clean.cmd

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


echo Build js-api ...
cd %package_dir%/../../js-api
call npm run build

echo Build  utils ...
cd %package_dir%/../../libraries/utils
call npm run build

echo Build bio ...
cd %package_dir%/../../libraries/bio
call npm run build

echo Build Bio ...
cd %package_dir%/../../packages/Bio
call npm run build

echo Build MLB-Data ...
cd %package_dir%/../../packages/MolecularLiabilityBrowserData
call npm run build

echo Build MLB ...
cd %package_dir%
call npm run build
