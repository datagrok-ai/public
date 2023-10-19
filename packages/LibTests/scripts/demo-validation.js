//name: validationDemo
//description: Custom validators demo
//language: javascript
//tags: simulation, demo
//input: double a {caption: First value; validator: LibTests:RangeValidatorFactory; validatorOptions: { "min": 1, "max": 10 } }
//input: double b {caption: Second value; validator: LibTests:RangeValidatorFactory; validatorOptions: { "min": 20, "max": 100 }  }
//input: double x=0 {caption: Async suggestions; validator: LibTests:AsyncValidatorDemoFactory }
//editor: Compute:RichFunctionViewEditor
//output: double c

c = a + b;
