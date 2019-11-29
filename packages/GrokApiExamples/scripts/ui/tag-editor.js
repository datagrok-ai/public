let v = gr.newView('demo: tag editor');

let editor = TagEditor.create();
editor.addTag('demo');
editor.addTag('test');
editor.addTag(1234);

v.append(editor.root);
