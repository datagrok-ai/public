set package_dir=%cd%

call npm uninstall -g datagrok-api @datagrok-libraries/utils @datagrok-libraries/bio @datagrok/bio @datagrok/molecular-liability-browser-data

cd "%package_dir%/../../js-api/"                                 & rmdir /s /q "node_modules"
cd "%package_dir%/../../libraries/utils/"                        & rmdir /s /q "node_modules"
cd "%package_dir%/../../libraries/ml/"                           & rmdir /s /q "node_modules"
cd "%package_dir%/../../libraries/bio/"                          & rmdir /s /q "node_modules"
cd "%package_dir%/../../packages/Bio/"                           & rmdir /s /q "node_modules"
cd "%package_dir%/../../packages/MolecularLiabilityBrowserData"  & rmdir /s /q "node_modules"
cd "%package_dir%/"                                              & rmdir /s /q "node_modules"

rem cd "%package_dir%/../../js-api/"                                 & del "package-lock.json"
rem cd "%package_dir%/../../libraries/utils/"                        & del "package-lock.json"
rem cd "%package_dir%/../../libraries/ml/"                           & del "package-lock.json"
rem cd "%package_dir%/../../libraries/bio/"                          & del "package-lock.json"
rem cd "%package_dir%/../../packages/Bio/"                           & del "package-lock.json"
rem cd "%package_dir%/../../packages/MolecularLiabilityBrowserData"  & del "package-lock.json"
rem cd "%package_dir%/"                                              & del "package-lock.json"
