call npm unlink datagrok-api
call npm unlink @datagrok-tools/utils
cd ../../js-api
call npm install
call npm link
cd ../libraries/utils
call npm install
call npm link
cd ../../packages/Chem
call npm install
call npm link datagrok-api @datagrok-libraries/utils
webpack