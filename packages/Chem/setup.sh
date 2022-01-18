npm unlink datagrok-api
npm unlink @datagrok-libraries/utils
npm unlink @datagrok-libraries/ml
cd ../../js-api
npm install
npm link
cd ../libraries/utils
npm install
npm link
cd ../../libraries/ml
npm install
npm link
cd ../../packages/Chem
npm install
npm link datagrok-api @datagrok-libraries/utils @datagrok-libraries/ml 
