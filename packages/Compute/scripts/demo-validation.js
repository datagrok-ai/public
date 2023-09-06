//name: Validation demo
//description: Custom validators demo
//language: javascript
//tags: simulation, demo
//input: double a {caption: First value; validator: Compute:RangeValidatorFactory; validatorOptions: { "min": 1, "max": 10 } }
//input: double b {caption: Second value; validator: Compute:RangeValidatorFactory; validatorOptions: { "min": 20, "max": 100 }  }
//editor: Compute:RichFunctionViewEditor
//output: double c

c = a + b
