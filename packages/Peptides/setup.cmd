call npm unlink datagrok-api
call npm unlink @datagrok-libraries/utils
call npm unlink @datagrok-libraries/ml
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
call npm link @datagrok-libraries/utils
cd ../../packages/Peptides
call npm install
call npm link datagrok-api @datagrok-libraries/utils @datagrok-libraries/ml
webpack