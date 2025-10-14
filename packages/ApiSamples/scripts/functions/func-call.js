let f = DG.Func.find({name: 'Sin'})[0];

// inspect function parameters
let inputs = f.inputs.map(input => input.propertyType + ' ' + input.name);
let outputs = f.outputs.map(output => output.propertyType + ' ' + output.name);
grok.shell.info(f.name +  '(' + inputs + ')' + ': ' + outputs);

// simple way: f.apply
// Any function can be called asynchronously
let result = await f.apply({x: 0.5});
grok.shell.info(result);
// If we know that this is sync function, AND it's either a core function or a plugin is loaded already,
// then we can call it synchronously using applySync:
grok.shell.info('Sync call: ' + f.applySync({x: 0.5}));

// more complex: preparing a FuncCall that lets you track inputs and outputs, intercept execution, etc
let fc = f.prepare({"x": 0.5});

grok.shell.info('Input name and type: ' + fc.inputParams['x'].name + ' ' + fc.inputParams['x'].property.propertyType);
grok.shell.info('Input value:' + fc.inputs['x']);

grok.shell.info(fc.inputs['x']);
let call = await fc.call();
grok.shell.info(call.outputs['result'])