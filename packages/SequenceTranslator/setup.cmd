cd ../../js-api
call npm install
call npm link
cd ../libraries/utils
call npm install
call npm link
call npm link datagrok-api
cd ../../packages/SequenceTranslator
call npm install
call npm link datagrok-api @datagrok-libraries/utils
webpack