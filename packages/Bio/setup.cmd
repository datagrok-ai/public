cd ../../js-api
call npm install
call npm link
cd ../libraries/utils
call npm install
call npm link
call npm link datagrok-api
cd ../libraries/ml
call npm install
call npm link
call npm link @datagrok-libraries/utils
cd ../libraries/bio
call npm install
call npm link
call npm link @datagrok-libraries/utils
cd ../../packages/Bio
call npm install
call npm link datagrok-api @datagrok-libraries/bio @datagrok-libraries/utils @datagrok-libraries/ml
webpack
