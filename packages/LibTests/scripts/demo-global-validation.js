//name: globalValidationTest
//description: Global validators
//language: javascript
//tags: simulation
//input: double a {caption: First value; validatorFunc: LibTests:GlobalValidatorDemoFactory; validatorOptions: { "max": 100 } }
//input: double b {caption: Second value; validatorFunc: LibTests:GlobalValidatorDemoFactory; validatorOptions: { "max": 100 } }
//input: double c {caption: Third value; validatorFunc: LibTests:GlobalValidatorDemoFactory; validatorOptions: { "max": 100 } }
//editor: Compute:RichFunctionViewEditor
//output: double s

s = a + b + c;
