//name: validationTest
//description: Custom validators test
//language: javascript
//tags: simulation
//input: double a {caption: First value; validator: LibTests:RangeValidatorFactory; validatorOptions: { "min": 1, "max": 10 } }
//input: double b {caption: Second value; validator: LibTests:RangeValidatorFactory; validatorOptions: { "min": 20, "max": 100 }  }
//input: double x=0 {caption: Async suggestions; validator: LibTests:AsyncValidatorDemoFactory }
//editor: Compute:RichFunctionViewEditor
//output: double c

c = a + b;
