cd ../../js-api
call npm install
call npm link
cd ../libraries/utils
call npm install
call npm link
call npm link datagrok-api
cd ../../libraries/ml
call npm install
call npm link
call npm link datagrok-api @datagrok-libraries/utils
cd ../../packages/Chem
call npm install
call npm link datagrok-api @datagrok-libraries/utils @datagrok-libraries/ml
webpack
