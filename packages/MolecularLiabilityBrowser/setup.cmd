set package_dir=%cd%

call npm uninstall -g datagrok-api @datagrok-libraries/utils @datagrok-libraries/bio @datagrok/bio @datagrok/molecular-liability-browser-data

cd "%package_dir%/../../js-api/"                                 & rmdir /s /q "node_modules"
cd "%package_dir%/../../libraries/utils/"                        & rmdir /s /q "node_modules"
cd "%package_dir%/../../libraries/bio/"                          & rmdir /s /q "node_modules"
cd "%package_dir%/../../packages/Bio/"                           & rmdir /s /q "node_modules"
cd "%package_dir%/../../packages/MolecularLiabilityBrowserData"  & rmdir /s /q "node_modules"
cd "%package_dir%/"                                              & rmdir /s /q "node_modules"

pause

cd "%package_dir%/../../js-api/"                                 & del "package-lock.json"
cd "%package_dir%/../../libraries/utils/"                        & del "package-lock.json"
cd "%package_dir%/../../libraries/bio/"                          & del "package-lock.json"
cd "%package_dir%/../../packages/Bio/"                           & del "package-lock.json"
cd "%package_dir%/../../packages/MolecularLiabilityBrowserData"  & del "package-lock.json"
cd "%package_dir%/"                                              & del "package-lock.json"

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
