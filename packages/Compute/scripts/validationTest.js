//name: validationTest
//description: Custom validators test
//language: javascript
//tags: simulation, test
//input: double a {caption: First value; validatorFunc: Compute:RangeValidatorFactory; validatorFuncOptions: { "min": 1, "max": 10 } }
//input: double b {caption: Second value; validatorFunc: Compute:RangeValidatorFactory; validatorFuncOptions: { "min": 20, "max": 100 }  }
//input: double x=0 {caption: Async suggestions; validatorFunc: Compute:AsyncValidatorDemoFactory }
//editor: Compute:RichFunctionViewEditor
//output: double c

c = a + b;
