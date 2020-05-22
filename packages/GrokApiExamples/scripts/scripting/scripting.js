// An example of using scripting (R, Python etc.)

grok.functions.call('PythonScripts:PythonDup', {'s': '1010'}).then(result => grok.shell.info(result));
grok.functions.call('RScripts:RDup', {'s': '0101'}).then(result => grok.shell.info(result));
