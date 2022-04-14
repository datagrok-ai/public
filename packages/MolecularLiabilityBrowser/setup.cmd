set package_dir=%cd%

cd ../../js-api
call npm install
call npm link
cd %package_dir%
call npm install
call npm link datagrok-api
webpack