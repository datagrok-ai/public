// a view that hosts multiple other views

let custom = DG.View.fromRoot(ui.div([]));
custom.toolbox = ui.divText('my toolbox');
custom.setRibbonPanels([[ui.button('hey')]]);

let multiView = new DG.MultiView({
  viewFactories: {
    'custom': () => custom,
    'demog': () => DG.TableView.create(grok.data.demo.demog(), false),
    'apps': () => DG.View.createByType(DG.View.APPS),
    'users': () => DG.View.createByType(DG.View.USERS),
  }
});

grok.shell.addView(multiView);