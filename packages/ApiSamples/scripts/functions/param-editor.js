let f = DG.Func.find({name: 'Sin'})[0];
grok.shell.o = await f.prepare().getEditor();