npm unlink datagrok-api
npm unlink @datagrok-tools/utils
cd ../../js-api
npm link
cd ../libraries/utils
npm link
cd ../../packages/Chem
npm link datagrok-api @datagrok-libraries/utils
npm install
webpack