//tags: Script, Package
//help-url: https://datagrok.ai/help/compute/scripting
// An example of using scripting (R, Python etc.)

let result1 = await grok.functions.call('Samples:PythonDup', {'s': '1010'});
let result2 = await grok.functions.call('Samples:RDup', {'s': '0101'});
grok.shell.info(result1);
grok.shell.info(result2)
