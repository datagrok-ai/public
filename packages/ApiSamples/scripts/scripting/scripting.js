//tags: Script, Package
//help-url: https://datagrok.ai/help/compute/scripting
// An example of using scripting (R, Python etc.)

grok.functions.call('Samples:PythonDup', {'s': '1010'}).then(result => grok.shell.info(result));
grok.functions.call('Samples:RDup', {'s': '0101'}).then(result => grok.shell.info(result));
