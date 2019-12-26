// An example of using scripting (R, Python etc.)

gr.callFunc('PythonScripts:PythonDup', {'s': '1010'}).then(result => gr.balloon.info(result));
gr.callFunc('RScripts:RDup', {'s': '0101'}).then(result => gr.balloon.info(result));
