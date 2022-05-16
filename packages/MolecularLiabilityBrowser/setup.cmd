set package_dir=%cd%
call npm unlink datagrok-api @datagrok-libraries/utils @datagrok-libraries/bio

rem echo js-api ...
rem cd ../../js-api
rem call npm install
rem call npm run build
rem call npm link

echo Utils ...
cd %package_dir%/../../libraries/utils
rmdir /s /q node_modules
call npm install
call npm run build
call npm link

echo Bio ...
cd %package_dir%/../../libraries/bio
rmdir /s /q node_modules
call npm install
call npm run build
call npm link

echo MLB ...
cd %package_dir%
rmdir /s /q node_modules
call npm install
call npm run link-all
call npm run build
