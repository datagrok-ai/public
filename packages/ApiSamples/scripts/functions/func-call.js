let f = DG.Func.find({name: 'Sin'})[0];

// inspect function parameters
let s = f.name + '(' + f.inputs.map(input => input.propertyType + ' ' + input.name) + ')' + ': ' + f.outputs.map(output => output.propertyType + ' ' + output.name);

grok.shell.info(s);

// simple way: f.apply
f.apply({x: 0.5}).then((result) => grok.shell.info(result));

// more complex: preparing a FuncCall that lets you track inputs and outputs, intercept execution, etc
let fc = f.prepare({"x": 0.5});

grok.shell.info(fc.inputs['x']);
fc.call().then((call) => grok.shell.info(call.outputs['result']));