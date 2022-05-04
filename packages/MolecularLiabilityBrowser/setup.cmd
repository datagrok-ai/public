set package_dir=%cd%

cd ../../js-api
call npm install
call npm link
cd %package_dir%
rmdir /s /q node_modules
call npm install
call npm link datagrok-api
call npm run build