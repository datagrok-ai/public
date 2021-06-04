let f = DG.Func.find({name: 'Sin'})[0];
f.prepare().getEditor().then((editor) => grok.shell.o = editor);