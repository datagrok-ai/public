cd ../../js-api
call npm install
call npm link
cd ../libraries/bio
call npm install
call npm link
cd ../../packages/Bio
call npm install
call npm link datagrok-api @datagrok-libraries/bio
webpack
