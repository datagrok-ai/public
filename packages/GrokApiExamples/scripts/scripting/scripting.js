// An example of using scripting (R, Python etc.)

grok.functions.callFunc('PythonScripts:PythonDup', {'s': '1010'}).then(result => grok.shell.info(result));
grok.functions.callFunc('RScripts:RDup', {'s': '0101'}).then(result => grok.shell.info(result));
