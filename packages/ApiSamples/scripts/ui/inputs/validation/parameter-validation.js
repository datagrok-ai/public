//name: parameter-validation
//language: javascript
//input: int foo = 5 { validator: bar > 3 }
//input: double bar = 2 { min: 0; max: 10 }

grok.shell.info('foo: ' + foo ?? '' + '\n' + 'bar: ' + bar ?? '');