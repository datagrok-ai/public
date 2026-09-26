let view = grok.shell.newView('hello');
view.append(ui.divText('hi'));
let sp = DG.Viewer.fromType('Scatter plot', grok.data.demo.demog());
grok.shell.newView('again');
