//name: globalValidationTest
//description: Global validators test
//language: javascript
//tags: simulation, test
//input: double a {caption: First value; validatorFunc: Compute:GlobalValidatorDemoFactory; validatorFuncOptions: { "max": 100 } }
//input: double b {caption: Second value; validatorFunc: Compute:GlobalValidatorDemoFactory; validatorFuncOptions: { "max": 100 } }
//input: double c {caption: Third value; validatorFunc: Compute:GlobalValidatorDemoFactory; validatorFuncOptions: { "max": 100 } }
//editor: Compute:RichFunctionViewEditor
//output: double s

s = a + b + c;
