set package_dir=%cd%

call setup-unlink-clean.cmd

cd "%package_dir%/../../js-api/"                                 & call npm install
cd "%package_dir%/../../libraries/utils/"                        & call npm install
cd "%package_dir%/../../libraries/ml/"                           & call npm install
cd "%package_dir%/../../libraries/bio/"                          & call npm install
cd "%package_dir%/../../packages/Bio/"                           & call npm install
cd "%package_dir%/../../packages/MolecularLiabilityBrowserData"  & call npm install
cd "%package_dir%/"                                              & call npm install

cd "%package_dir%/../../js-api/"                                 & call npm link
cd "%package_dir%/../../libraries/utils/"                        & call npm link
cd "%package_dir%/../../libraries/ml/"                           & call npm link
cd "%package_dir%/../../libraries/bio/"                          & call npm link
cd "%package_dir%/../../packages/Bio/"                           & call npm link
cd "%package_dir%/../../packages/MolecularLiabilityBrowserData"  & call npm link
cd "%package_dir%/"                                              & call npm link

cd "%package_dir%/../../js-api/"                                 & call npm run link-all
cd "%package_dir%/../../libraries/utils/"                        & call npm run link-all
cd "%package_dir%/../../libraries/ml/"                           & call npm run link-all
cd "%package_dir%/../../libraries/bio/"                          & call npm run link-all
cd "%package_dir%/../../packages/Bio/"                           & call npm run link-all
cd "%package_dir%/../../packages/MolecularLiabilityBrowserData"  & call npm run link-all
cd "%package_dir%/"                                              & call npm run link-all

cd "%package_dir%/../../js-api/"                                 & call npm run build
cd "%package_dir%/../../libraries/utils/"                        & call npm run build
cd "%package_dir%/../../libraries/ml/"                           & call npm run build
cd "%package_dir%/../../libraries/bio/"                          & call npm run build
cd "%package_dir%/../../packages/Bio/"                           & call npm run build
cd "%package_dir%/../../packages/MolecularLiabilityBrowserData"  & call npm run build
cd "%package_dir%/"                                              & call npm run build
