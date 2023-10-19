//name: globalValidationDemo
//description: Global validators demo
//language: javascript
//tags: simulation, demo
//input: double a {caption: First value; validator: LibTests:GlobalValidatorDemoFactory; validatorOptions: { "max": 100 } }
//input: double b {caption: Second value; validator: LibTests:GlobalValidatorDemoFactory; validatorOptions: { "max": 100 } }
//input: double c {caption: Third value; validator: LibTests:GlobalValidatorDemoFactory; validatorOptions: { "max": 100 } }
//editor: Compute:RichFunctionViewEditor
//output: double s

s = a + b + c;
