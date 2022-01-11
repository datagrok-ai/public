//tags: Script, Package
//help-url: https://datagrok.ai/help/compute/scripting
// An example of using scripting (R, Python etc.)

grok.functions.call('DemoScripts:PythonDup', {'s': '1010'}).then(result => grok.shell.info(result));
grok.functions.call('DemoScripts:RDup', {'s': '0101'}).then(result => grok.shell.info(result));
