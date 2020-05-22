// An example of using scripting (R, Python etc.)

grok.callFunc('PythonScripts:PythonDup', {'s': '1010'}).then(result => grok.shell.balloon.info(result));
grok.callFunc('RScripts:RDup', {'s': '0101'}).then(result => grok.shell.balloon.info(result));
