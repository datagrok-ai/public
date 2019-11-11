// An example of using scripting (R, Python etc.)

gr.callFunc('Demo:PythonScripts:PythonDup', {'s': '1010'}).then(result => gr.balloon.info(result));
gr.callFunc('Demo:RScripts:RDup', {'s': '0101'}).then(result => gr.balloon.info(result));
